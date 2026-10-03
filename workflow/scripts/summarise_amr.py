#!/usr/bin/env python3
"""Merge RGI (CARD) and ABRicate hits per MAG into one long table: mag, source, gene, drug_class, identity.

RGI: only 'Perfect' and 'Strict' hits are kept (Loose hits are excluded by default; --include-loose to keep).
ABRicate: identity/coverage thresholds were applied at run time.
"""
import argparse
import sys
from pathlib import Path

import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--rgi-dir", required=True)
    ap.add_argument("--abricate-dir", required=True)
    ap.add_argument("--include-loose", action="store_true")
    ap.add_argument("--out", required=True)
    a = ap.parse_args()
    rows = []
    for f in sorted(Path(a.rgi_dir).glob("*.rgi.txt")):
        mag = f.name[: -len(".rgi.txt")]
        d = pd.read_csv(f, sep="\t")
        if d.empty:
            continue
        if not a.include_loose:
            d = d[d["Cut_Off"].isin(["Perfect", "Strict"])]
        for _, r in d.iterrows():
            rows.append({"mag": mag, "source": "rgi", "gene": r["Best_Hit_ARO"], "drug_class": r.get("Drug Class"),
                         "identity": r.get("Best_Identities")})
    for f in sorted(Path(a.abricate_dir).glob("*.tsv")):
        mag, _, db = f.name[: -len(".tsv")].rpartition(".")
        d = pd.read_csv(f, sep="\t")
        for _, r in d.iterrows():
            rows.append({"mag": mag, "source": f"abricate_{db}", "gene": r["GENE"], "drug_class": r.get("RESISTANCE"),
                         "identity": r.get("%IDENTITY")})
    pd.DataFrame(rows, columns=["mag", "source", "gene", "drug_class", "identity"]).to_csv(a.out, sep="\t", index=False, na_rep="NA")
    print(f"[summarise_amr] {len(rows)} hits", file=sys.stderr)


if __name__ == "__main__":
    main()
