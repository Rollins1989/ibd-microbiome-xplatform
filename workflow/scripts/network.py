#!/usr/bin/env python3
"""SparCC co-occurrence networks per diagnosis group, with module detection and comparison.

Design choices
  * SparCC (compositionally aware) on raw counts, median of Dirichlet-resampled estimates.
  * Edges: |r| >= min_abs_r AND permutation-null FDR <= fdr (pooled null across all pairs).
  * Same feature set for every group (prevalence >= threshold in ALL groups) so the networks are comparable.
  * Default uses one baseline sample per subject; using all visits makes correlations pseudo-replicated
    and the FDR anti-conservative (the config option exists but is flagged in the output).
Exploratory by nature; interpret hubs/modules as hypotheses.
"""
import argparse
import sys
from itertools import combinations
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402
import networkx as nx  # noqa: E402
import numpy as np  # noqa: E402
import pandas as pd  # noqa: E402
import yaml  # noqa: E402

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402

COLORS = {"nonIBD": "#0072B2", "CD": "#D55E00", "UC": "#009E73"}


def build_graph(cor, q, feats, min_r, fdr):
    G = nx.Graph()
    G.add_nodes_from(feats)
    for i, j in combinations(range(len(feats)), 2):
        if abs(cor[i, j]) >= min_r and q[i, j] <= fdr:
            G.add_edge(feats[i], feats[j], r=float(cor[i, j]), fdr=float(q[i, j]))
    return G


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--config", required=True)
    ap.add_argument("--meta", required=True)
    ap.add_argument("--layer-dir", required=True)
    ap.add_argument("--outdir", required=True)
    ap.add_argument("--figdir", required=True)
    a = ap.parse_args()
    cfg = yaml.safe_load(open(a.config))
    nc, seed = cfg["networks"], cfg["stats"]["seed"]
    layer = nc["layer"]
    meta = pd.read_csv(a.meta, sep="\t", na_values=mb.NA_STRINGS, keep_default_na=False)
    counts = mb.read_layer(Path(a.layer_dir) / f"{layer}.tsv")
    plat = "16S" if layer.startswith("16s") else "MGX"
    m = meta[(meta.data_type == plat) & meta.sample_id.isin(counts.columns)]
    if nc["samples"] == "baseline":
        m = m[m.is_baseline.astype(str).str.lower() == "true"]
    groups = list(cfg["design"]["diagnosis_levels"])
    ids = {g: m.loc[m.diagnosis == g, "sample_id"].tolist() for g in groups}
    for g, v in ids.items():
        if len(v) < 20:
            print(f"[network] WARNING: group {g} has only {len(v)} samples; correlations will be unstable", file=sys.stderr)
    prev = pd.DataFrame({g: (counts[ids[g]] > 0).mean(axis=1) for g in groups})
    feats = prev.index[(prev >= nc["min_group_prevalence"]).all(axis=1)].tolist()
    if len(feats) < 5:
        raise SystemExit(f"only {len(feats)} features pass the prevalence filter in all groups")
    print(f"[network] layer={layer} features={len(feats)} samples={ {g: len(v) for g, v in ids.items()} }", file=sys.stderr)

    graphs, edge_rows, node_rows, metric_rows = {}, [], [], []
    for gi, g in enumerate(groups):
        X = counts.loc[feats, ids[g]].T.to_numpy()
        cor, q = mb.sparcc_with_fdr(X, nc["sparcc"]["iterations"], nc["sparcc"]["permutations"], seed + gi)
        G = build_graph(cor, q, feats, nc["min_abs_r"], nc["fdr"])
        graphs[g] = G
        mods = {}
        if G.number_of_edges():
            comms = nx.community.louvain_communities(G, seed=seed, weight=None)
            comms = sorted([c for c in comms if len(c) > 1], key=len, reverse=True)
            for k, c in enumerate(comms, 1):
                for n in c:
                    mods[n] = k
            modularity = nx.community.modularity(G, nx.community.louvain_communities(G, seed=seed))
        else:
            comms, modularity = [], np.nan
        deg = dict(G.degree())
        for n in feats:
            node_rows.append({"group": g, "feature": n, "degree": deg[n], "module": mods.get(n, 0)})
        for u, v, d in G.edges(data=True):
            edge_rows.append({"group": g, "feature_a": u, "feature_b": v, "r": d["r"], "fdr": d["fdr"]})
        hubs = sorted(deg, key=deg.get, reverse=True)[:5]
        pos = sum(1 for *_, d in G.edges(data=True) if d["r"] > 0)
        metric_rows.append({"group": g, "n_samples": len(ids[g]), "n_nodes": len(feats), "n_edges": G.number_of_edges(),
                            "density": nx.density(G), "pct_positive": 100 * pos / max(1, G.number_of_edges()),
                            "n_modules": len(comms), "modularity": modularity, "top_hubs": ";".join(hubs),
                            "sample_scheme": nc["samples"]})
    out = Path(a.outdir)
    out.mkdir(parents=True, exist_ok=True)
    pd.DataFrame(metric_rows).to_csv(out / "network_metrics.tsv", sep="\t", index=False)
    pd.DataFrame(edge_rows).to_csv(out / "network_edges.tsv", sep="\t", index=False)
    pd.DataFrame(node_rows).to_csv(out / "network_nodes.tsv", sep="\t", index=False)

    cmp_rows = []
    for g1, g2 in combinations(groups, 2):
        e1 = {frozenset((u, v)): d["r"] for u, v, d in graphs[g1].edges(data=True)}
        e2 = {frozenset((u, v)): d["r"] for u, v, d in graphs[g2].edges(data=True)}
        shared = set(e1) & set(e2)
        union = set(e1) | set(e2)
        flips = sum(1 for e in shared if np.sign(e1[e]) != np.sign(e2[e]))
        cmp_rows.append({"group_1": g1, "group_2": g2, "edges_1": len(e1), "edges_2": len(e2), "shared": len(shared),
                         "unique_1": len(set(e1) - set(e2)), "unique_2": len(set(e2) - set(e1)),
                         "jaccard": len(shared) / len(union) if union else np.nan, "sign_flips": flips})
    pd.DataFrame(cmp_rows).to_csv(out / "network_comparison.tsv", sep="\t", index=False)

    # ----- figure
    U = nx.Graph()
    U.add_nodes_from(feats)
    for G in graphs.values():
        U.add_edges_from(G.edges())
    pos = nx.spring_layout(U, seed=seed, k=1.5 / np.sqrt(len(feats)))
    fig, axes = plt.subplots(1, len(groups), figsize=(5.2 * len(groups), 5.4))
    pal = plt.get_cmap("tab10")
    nodes = pd.DataFrame(node_rows)
    for ax, g in zip(np.atleast_1d(axes), groups):
        G = graphs[g]
        nd = nodes[nodes.group == g].set_index("feature")
        col = [pal((nd.loc[n, "module"] - 1) % 10) if nd.loc[n, "module"] > 0 else "#cccccc" for n in G.nodes()]
        size = [40 + 60 * G.degree(n) for n in G.nodes()]
        ec = ["#B2182B" if d["r"] > 0 else "#2166AC" for *_, d in G.edges(data=True)]
        nx.draw_networkx_edges(G, pos, ax=ax, edge_color=ec, alpha=0.6, width=1.0)
        nx.draw_networkx_nodes(G, pos, ax=ax, node_color=col, node_size=size, linewidths=0.4, edgecolors="black")
        hubs = sorted(G.degree, key=lambda x: -x[1])[:5]
        nx.draw_networkx_labels(G, pos, labels={n: n for n, d in hubs if d > 0}, font_size=7, ax=ax)
        row = [r for r in metric_rows if r["group"] == g][0]
        ax.set_title(f"{g}  (n={row['n_samples']}, edges={row['n_edges']}, modules={row['n_modules']})",
                     color=COLORS.get(g, "black"), fontsize=10, fontweight="bold")
        ax.axis("off")
    fig.suptitle(f"SparCC networks ({layer}); red = positive, blue = negative; node colour = module; "
                 f"|r|>={nc['min_abs_r']}, FDR<={nc['fdr']}", fontsize=9)
    fig.tight_layout()
    Path(a.figdir).mkdir(parents=True, exist_ok=True)
    for ext in ("pdf", "png"):
        fig.savefig(Path(a.figdir) / f"Fig8_cooccurrence_networks.{ext}", dpi=300)
    print("[network] done", file=sys.stderr)


if __name__ == "__main__":
    main()
