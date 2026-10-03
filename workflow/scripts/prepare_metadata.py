#!/usr/bin/env python3
"""Validate the sample sheet and derive clinical variables (activity, pre-flare labels, flare events)."""
import argparse
import sys
from pathlib import Path

import yaml

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "lib"))
import mbstats as mb  # noqa: E402


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--config", required=True, help="YAML containing the 'design' block (resolved config dump)")
    ap.add_argument("--out-samples", required=True)
    ap.add_argument("--out-specimens", required=True)
    ap.add_argument("--out-events", required=True)
    a = ap.parse_args()

    design = yaml.safe_load(open(a.config))["design"]
    df = mb.load_samples(a.samples, design["diagnosis_levels"])
    meta, spec, ev = mb.derive_clinical(df, design)
    for df_, p in ((meta, a.out_samples), (spec, a.out_specimens), (ev, a.out_events)):
        Path(p).parent.mkdir(parents=True, exist_ok=True)
        df_.to_csv(p, sep="\t", index=False, na_rep="NA")
    n_pair = spec[["sample_16s", "sample_mgx"]].notna().all(axis=1).sum() if {"sample_16s", "sample_mgx"} <= set(spec) else 0
    print(f"[prepare_metadata] samples={len(meta)} specimens={len(spec)} subjects={spec.subject_id.nunique()} "
          f"paired_specimens={n_pair} flare_events={ev.event_id.nunique()}", file=sys.stderr)


if __name__ == "__main__":
    main()
