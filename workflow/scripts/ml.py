#!/usr/bin/env python3
"""Subject-grouped nested cross-validation: 16S vs shotgun vs combined feature sets (Figure 5).

Leakage control
  * Folds are grouped by subject (repeated visits never straddle train/test), outer AND inner.
  * Prevalence filtering, pseudocount choice and CLR are fitted inside every training fold
    (CompositionalBlocks) - nothing is learned from test data.
  * Hyper-parameters are tuned only in the inner loop.
  * Feature sets are compared on the same specimens when ml.paired_only is true.
Confidence intervals: subject-level (cluster) bootstrap of the repeat-averaged out-of-fold predictions.
SHAP is computed on a model refit on all data and is for interpretation only, not for performance.
"""
import argparse
import json
import sys
import warnings
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402
from sklearn.base import BaseEstimator, TransformerMixin  # noqa: E402
from sklearn.linear_model import LogisticRegression  # noqa: E402
from sklearn.ensemble import RandomForestClassifier  # noqa: E402
from sklearn.metrics import average_precision_score, roc_auc_score, roc_curve  # noqa: E402
from sklearn.model_selection import GridSearchCV, StratifiedGroupKFold  # noqa: E402
from sklearn.pipeline import Pipeline  # noqa: E402
from sklearn.preprocessing import StandardScaler  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402

warnings.filterwarnings("ignore")
FS_COLORS = {"16s": "#CC6677", "mgx": "#44AA99", "combined": "#332288"}
PLAT_LABEL = {"16s": "16S", "mgx": "MGX"}


class CompositionalBlocks(BaseEstimator, TransformerMixin):
    """Per-block prevalence filter -> closure -> pseudocount -> CLR, all fitted on training data only."""

    def __init__(self, blocks=None, feature_names=None, min_prev=0.1):
        self.blocks = blocks
        self.feature_names = feature_names
        self.min_prev = min_prev

    def fit(self, X, y=None):
        X = np.asarray(X, float)
        self.keep_, self.pseudo_ = {}, {}
        for b, idx in self.blocks.items():
            idx = np.asarray(idx)
            prev = (X[:, idx] > 0).mean(axis=0)
            keep = idx[prev >= self.min_prev]
            if len(keep) < 2:                                  # never return an empty block
                keep = idx[np.argsort(-prev)[: min(2, len(idx))]]
            sub = X[:, keep]
            sub = sub / np.where(sub.sum(1, keepdims=True) > 0, sub.sum(1, keepdims=True), 1)
            pos = sub[sub > 0]
            self.keep_[b], self.pseudo_[b] = keep, (pos.min() / 2 if pos.size else 1e-6)
        return self

    def transform(self, X):
        X = np.asarray(X, float)
        out = []
        for b in self.blocks:
            sub = X[:, self.keep_[b]]
            s = sub.sum(1, keepdims=True)
            sub = sub / np.where(s > 0, s, 1) + self.pseudo_[b]
            sub = sub / sub.sum(1, keepdims=True)
            lx = np.log(sub)
            out.append(lx - lx.mean(1, keepdims=True))
        return np.hstack(out)

    def get_feature_names_out(self, input_features=None):
        names = np.asarray(self.feature_names)
        return np.concatenate([names[self.keep_[b]] for b in self.blocks])


def make_model(name, prep, n_jobs, seed):
    steps = [("prep", prep)]
    if name == "rf":
        steps.append(("clf", RandomForestClassifier(n_estimators=300, class_weight="balanced_subsample",
                                                     n_jobs=n_jobs, random_state=seed)))
        grid = {"clf__max_features": ["sqrt", 0.2], "clf__min_samples_leaf": [1, 5]}
    elif name == "xgb":
        from xgboost import XGBClassifier
        steps.append(("clf", XGBClassifier(n_estimators=300, learning_rate=0.05, subsample=0.8, colsample_bytree=0.5,
                                           eval_metric="logloss", tree_method="hist", n_jobs=n_jobs,
                                           random_state=seed, verbosity=0)))
        grid = {"clf__max_depth": [2, 4], "clf__min_child_weight": [1, 5]}
    elif name == "logreg":
        steps += [("scale", StandardScaler()),
                  ("clf", LogisticRegression(penalty="elasticnet", solver="saga", max_iter=5000,
                                             class_weight="balanced", random_state=seed))]
        grid = {"clf__C": [0.01, 0.1, 1.0], "clf__l1_ratio": [0.2, 0.8]}
    else:
        raise ValueError(name)
    return Pipeline(steps), grid


def task_labels(spec: pd.DataFrame, task: str):
    d = spec.copy()
    if task == "IBD_vs_nonIBD":
        d["y"] = (d.diagnosis != "nonIBD").astype(float)
    elif task in ("CD_vs_nonIBD", "UC_vs_nonIBD"):
        pos = task.split("_")[0]
        d = d[d.diagnosis.isin([pos, "nonIBD"])]
        d["y"] = (d.diagnosis == pos).astype(float)
    elif task == "CD_vs_UC":
        d = d[d.diagnosis.isin(["CD", "UC"])]
        d["y"] = (d.diagnosis == "CD").astype(float)
    elif task == "preflare":
        d = d[d.diagnosis != "nonIBD"]
        d["y"] = d["flare_next"]
    else:
        raise ValueError(f"unknown task {task}")
    return d.dropna(subset=["y"])


def build_matrix(d, layers, layer_dir, fs_layers):
    """Rows = specimens in d; columns = proportions from the requested layers (block-prefixed names)."""
    cols, blocks, start = [], {}, 0
    for lname in fs_layers:
        L = mb.read_layer(Path(layer_dir) / f"{lname}.tsv")
        P = mb.to_proportions(L)
        plat = "sample_16s" if lname.startswith("16s") else "sample_mgx"
        sub = P.reindex(columns=d[plat].tolist()).T
        sub.index = d["specimen_id"].tolist()
        tag = "16S" if lname.startswith("16s") else "MGX"
        sub.columns = [f"[{tag}] {c}" for c in sub.columns]
        cols.append(sub)
        blocks[lname] = list(range(start, start + sub.shape[1]))
        start += sub.shape[1]
    X = pd.concat(cols, axis=1)
    return X, blocks


def permute_labels(y, groups, rng):
    """Shuffle complete subject-level label vectors without breaking within-subject structure.

    Constant-within-subject labels are permuted between subjects. For longitudinal labels
    that vary by visit (for example ``preflare``), complete label vectors are permuted
    between subjects with the same number of observations. This avoids the previous
    specimen-level fallback, which destroyed the within-subject dependence structure.
    """
    groups = np.asarray(groups)
    y = np.asarray(y)
    subjects = pd.unique(groups)
    idx_by_subject = {s: np.flatnonzero(groups == s) for s in subjects}
    labels_by_subject = {s: y[idx_by_subject[s]].copy() for s in subjects}

    if all(np.all(v == v[0]) for v in labels_by_subject.values()):
        shuffled = rng.permutation(subjects)
        out = np.empty_like(y)
        for target, source in zip(subjects, shuffled):
            out[idx_by_subject[target]] = labels_by_subject[source][0]
        return out

    # Longitudinal labels may vary within a subject. Only exchange complete
    # vectors among subjects with the same visit count so every row keeps its
    # original within-subject position and vector length.
    out = y.copy()
    by_n = {}
    for s, idx in idx_by_subject.items():
        by_n.setdefault(len(idx), []).append(s)
    for subject_group in by_n.values():
        shuffled = rng.permutation(subject_group)
        for target, source in zip(subject_group, shuffled):
            out[idx_by_subject[target]] = labels_by_subject[source]
    return out


def cluster_bootstrap_auc(y, p, groups, B, rng):
    subj = np.unique(groups)
    idx_by = {s: np.where(groups == s)[0] for s in subj}
    vals = []
    for _ in range(B):
        pick = rng.choice(subj, len(subj), replace=True)
        ii = np.concatenate([idx_by[s] for s in pick])
        if len(np.unique(y[ii])) == 2:
            vals.append(roc_auc_score(y[ii], p[ii]))
    return (np.percentile(vals, 2.5), np.percentile(vals, 97.5)) if len(vals) > 10 else (np.nan, np.nan)


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--specimens", required=True)
    ap.add_argument("--layer-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--figdir", required=True)
    ap.add_argument("--permute-labels", action="store_true",
                    help="negative control: shuffle labels at subject level (expected AUROC ~ 0.5)")
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    mc, seed = cfg["ml"], cfg["stats"]["seed"]
    out = Path(a.outdir); out.mkdir(parents=True, exist_ok=True)
    Path(a.figdir).mkdir(parents=True, exist_ok=True)
    spec = pd.read_csv(a.specimens, sep="\t", na_values=mb.NA_STRINGS, keep_default_na=False)
    rng = np.random.default_rng(seed)
    fsets = {k: v for k, v in mc["feature_sets"].items()}
    need_plat = lambda ls: ["sample_16s" if l.startswith("16s") else "sample_mgx" for l in ls]  # noqa: E731
    all_plat = sorted({p for v in fsets.values() for p in need_plat(v)})

    results, oof_rows, status, shap_rows = [], [], [], []
    roc_store, shap_store = {}, {}
    for task in mc["tasks"]:
        dt = task_labels(spec, task)
        for fs, fl in fsets.items():
            plats = all_plat if mc["paired_only"] else need_plat(fl)
            d = dt.dropna(subset=plats).reset_index(drop=True)
            n_subj = d.groupby("y")["subject_id"].nunique()
            if len(n_subj) < 2 or n_subj.min() < mc["outer_folds"]:
                status.append({"task": task, "feature_set": fs, "status": "skipped",
                               "reason": f"too few subjects per class: {n_subj.to_dict()}"})
                continue
            X, blocks = build_matrix(d, fl, a.layer_dir, fl)
            y, g = d["y"].to_numpy().astype(int), d["subject_id"].to_numpy()
            if a.permute_labels:
                y = permute_labels(y, g, np.random.default_rng(seed))
            status.append({"task": task, "feature_set": fs, "status": "ok", "n_specimens": len(d),
                           "n_subjects": int(d.subject_id.nunique()), "n_pos": int(y.sum()), "n_neg": int((1 - y).sum())})
            print(f"[ml] {task} / {fs}: n={len(d)} subjects={d.subject_id.nunique()} pos={y.sum()}", file=sys.stderr)
            for model in mc["models"]:
                oof = np.full((mc["repeats"], len(y)), np.nan)
                for r in range(mc["repeats"]):
                    outer = StratifiedGroupKFold(mc["outer_folds"], shuffle=True, random_state=seed + r)
                    for k, (tr, te) in enumerate(outer.split(X, y, g)):
                        prep = CompositionalBlocks(blocks, list(X.columns), mc["min_prevalence"])
                        pipe, grid = make_model(model, prep, mc["n_jobs"], seed + r)
                        inner_k = int(min(mc["inner_folds"], np.bincount(y[tr]).min(), len(np.unique(g[tr]))))
                        gs = GridSearchCV(pipe, grid, scoring="roc_auc", refit=True, n_jobs=1,
                                          cv=StratifiedGroupKFold(max(2, inner_k), shuffle=True, random_state=seed + r))
                        gs.fit(X.iloc[tr].to_numpy(), y[tr], groups=g[tr])
                        oof[r, te] = gs.predict_proba(X.iloc[te].to_numpy())[:, 1]
                for r in range(mc["repeats"]):
                    results.append({"task": task, "feature_set": fs, "model": model, "repeat": r,
                                    "auroc": roc_auc_score(y, oof[r]), "auprc": average_precision_score(y, oof[r])})
                pbar = oof.mean(axis=0)
                lo, hi = cluster_bootstrap_auc(y, pbar, g, mc["bootstrap"], rng)
                roc_store[(task, fs, model)] = (y, pbar, roc_auc_score(y, pbar), lo, hi)
                oof_rows.append(pd.DataFrame({"task": task, "feature_set": fs, "model": model,
                                              "specimen_id": d.specimen_id, "subject_id": g, "y": y, "p_oof_mean": pbar}))
    res = pd.DataFrame(results)
    if res.empty:
        pd.DataFrame(status).to_csv(out / "ml_status.tsv", sep="\t", index=False)
        raise SystemExit("no task had enough data; see ml_status.tsv")
    summ = (res.groupby(["task", "feature_set", "model"])
            .agg(auroc_mean=("auroc", "mean"), auroc_sd=("auroc", "std"), auprc_mean=("auprc", "mean"),
                 repeats=("repeat", "nunique")).reset_index())
    ci = pd.DataFrame([{"task": k[0], "feature_set": k[1], "model": k[2], "auroc_oof_mean_pred": v[2],
                        "ci_low": v[3], "ci_high": v[4]} for k, v in roc_store.items()])
    summ = summ.merge(ci, on=["task", "feature_set", "model"])
    res.to_csv(out / "ml_cv_results.tsv", sep="\t", index=False)
    summ.to_csv(out / "ml_summary.tsv", sep="\t", index=False)
    pd.concat(oof_rows).to_csv(out / "ml_oof_predictions.tsv", sep="\t", index=False)
    pd.DataFrame(status).to_csv(out / "ml_status.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ SHAP (final fit on all data)
    shap_model = mc["shap_model"]
    if shap_model in mc["models"] and not a.permute_labels:
        import shap
        for task in sorted(res.task.unique()):
            fs = "combined" if "combined" in fsets and ((summ.task == task) & (summ.feature_set == "combined")).any() else None
            if fs is None:
                continue
            plats = all_plat if mc["paired_only"] else need_plat(fsets[fs])
            d = task_labels(spec, task).dropna(subset=plats).reset_index(drop=True)
            X, blocks = build_matrix(d, fsets[fs], a.layer_dir, fsets[fs])
            y, g = d["y"].to_numpy().astype(int), d["subject_id"].to_numpy()
            prep = CompositionalBlocks(blocks, list(X.columns), mc["min_prevalence"])
            pipe, grid = make_model(shap_model, prep, mc["n_jobs"], seed)
            gs = GridSearchCV(pipe, grid, scoring="roc_auc", n_jobs=1,
                              cv=StratifiedGroupKFold(max(2, mc["inner_folds"]), shuffle=True, random_state=seed))
            gs.fit(X.to_numpy(), y, groups=g)
            best = gs.best_estimator_
            Xt = best.named_steps["prep"].transform(X.to_numpy())
            names = best.named_steps["prep"].get_feature_names_out()
            sv = shap.TreeExplainer(best.named_steps["clf"]).shap_values(Xt)
            if isinstance(sv, list):
                sv = sv[1]
            elif getattr(sv, "ndim", 2) == 3:
                sv = sv[:, :, 1]
            imp = pd.DataFrame({"task": task, "feature": names, "mean_abs_shap": np.abs(sv).mean(0),
                                "mean_shap_signed_corr": [np.corrcoef(Xt[:, i], sv[:, i])[0, 1] if sv[:, i].std() > 0 else 0.0
                                                          for i in range(len(names))]})
            imp = imp.sort_values("mean_abs_shap", ascending=False)
            shap_rows.append(imp)
            shap_store[task] = imp.head(15)
            plt.figure()
            shap.summary_plot(sv, Xt, feature_names=list(names), max_display=15, show=False)
            plt.title(f"SHAP: {task} ({shap_model}, combined features)", fontsize=9)
            plt.tight_layout()
            plt.savefig(Path(a.figdir) / f"Fig5_shap_beeswarm_{task}.png", dpi=250)
            plt.close()
        if shap_rows:
            pd.concat(shap_rows).to_csv(out / "ml_shap_importance.tsv", sep="\t", index=False)

    # ------------------------------------------------------------------ Figure 5
    tasks = sorted(res.task.unique())
    nrow = 2 if shap_store else 1
    fig, axes = plt.subplots(nrow, len(tasks), figsize=(4.6 * len(tasks), 4.4 * nrow), squeeze=False)
    for j, task in enumerate(tasks):
        ax = axes[0, j]
        ax.plot([0, 1], [0, 1], ls="--", c="grey", lw=0.8)
        for fs in fsets:
            sub = summ[(summ.task == task) & (summ.feature_set == fs)]
            if sub.empty:
                continue
            best = sub.sort_values("auroc_oof_mean_pred", ascending=False).iloc[0]
            y_, p_, auc_, lo, hi = roc_store[(task, fs, best.model)]
            fpr, tpr, _ = roc_curve(y_, p_)
            ax.plot(fpr, tpr, c=FS_COLORS.get(fs, "k"), lw=1.8,
                    label=f"{fs} ({best.model}) {auc_:.2f} [{lo:.2f}-{hi:.2f}]")
        ax.set_title(task.replace("_", " "), fontsize=10, fontweight="bold")
        ax.set_xlabel("1 - specificity"); ax.set_ylabel("sensitivity")
        ax.legend(loc="lower right", fontsize=7, frameon=False)
        if shap_store and task in shap_store:
            ax2 = axes[1, j]
            t = shap_store[task].iloc[::-1]
            colr = ["#CC6677" if f.startswith("[16S]") else "#44AA99" for f in t.feature]
            ax2.barh(t.feature, t.mean_abs_shap, color=colr)
            ax2.set_xlabel("mean |SHAP|"); ax2.tick_params(axis="y", labelsize=7)
            ax2.set_title("top features (combined model)", fontsize=9)
        elif shap_store:
            axes[1, j].axis("off")
    fig.suptitle("Subject-grouped nested CV, AUROC [95% cluster-bootstrap CI]; red = 16S, green = shotgun features",
                 fontsize=9)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(Path(a.figdir) / f"Fig5_ml_performance.{ext}", dpi=300)
    print(summ.sort_values(["task", "auroc_oof_mean_pred"], ascending=[True, False]).to_string(index=False), file=sys.stderr)


if __name__ == "__main__":
    main()