"""Riaz counts (Entrez rows) -> log2 CPM keyed by Ensembl gene id.

Writes a genes x samples matrix restricted to Human-GEM genes, plus a tidy
sample table carrying treatment timepoint and best overall response.
"""
import argparse, sys
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--counts", required=True)
ap.add_argument("--clinical", required=True)
ap.add_argument("--gem-genes", required=True)
ap.add_argument("--out-expr", required=True)
ap.add_argument("--out-samples", required=True)
a = ap.parse_args()

counts = pd.read_csv(a.counts, sep="\t", index_col=0)
counts.index = counts.index.astype(str)        # Entrez ids
symbols = counts.pop("HUGO") if "HUGO" in counts.columns else None
counts = counts.apply(pd.to_numeric, errors="coerce").dropna(how="all")
print(f"counts: {counts.shape[0]} genes x {counts.shape[1]} samples", file=sys.stderr)

# Entrez -> Ensembl, from the Human-GEM gene annotation table.
gem = pd.read_csv(a.gem_genes, sep="\t", dtype=str)
gem = gem[["genes", "geneEntrezID", "geneSymbols"]].dropna(subset=["geneEntrezID"])
gem["geneEntrezID"] = gem["geneEntrezID"].str.split(";")
gem = gem.explode("geneEntrezID")
gem["geneEntrezID"] = gem["geneEntrezID"].str.strip()
e2g = gem.drop_duplicates("geneEntrezID").set_index("geneEntrezID")["genes"]

# CPM then log2: E-Flux needs a within-sample relative scale, so library-size
# normalisation is the only correction required (no cross-sample quantile step,
# which would blunt the between-patient contrast we are testing).
cpm = counts.div(counts.sum(axis=0), axis=1) * 1e6
logcpm = np.log2(cpm + 1.0)

shared = logcpm.index.intersection(e2g.index)
expr = logcpm.loc[shared].copy()
expr.index = e2g.loc[shared].values
expr = expr.groupby(level=0).max()          # collapse duplicate Ensembl ids
print(f"mapped to {expr.shape[0]} Human-GEM genes "
      f"({expr.shape[0]}/{gem['genes'].nunique()} of the model's genes)", file=sys.stderr)

# Sample table: column names are <PatientID>_<Pre|On>.
clin = pd.read_csv(a.clinical)
clin["PatientID"] = clin["PatientID"].astype(str)
resp = clin.drop_duplicates("PatientID").set_index("PatientID")

rows = []
for col in expr.columns:
    stem = col.split(".")[0]                 # Pt109_On.2 -> Pt109_On
    pid, _, tp = stem.rpartition("_")
    c = resp.loc[pid] if pid in resp.index else None
    rows.append({
        "sample": col,
        "patient": pid,
        "timepoint": tp,
        "replicate": "." in col,
        "BOR": None if c is None else c.get("BOR"),
        "response3": None if c is None else c.get("myBOR3"),
        "cohort": None if c is None else c.get("Cohort"),
        "PFS": None if c is None else c.get("PFS"),
        "PFS_event": None if c is None else c.get("PFS_SOR"),
        "OS": None if c is None else c.get("OS"),
        "OS_event": None if c is None else c.get("OS_SOR"),
    })
samples = pd.DataFrame(rows)
# Riaz codes PFS_SOR/OS_SOR as 1 = censored, so events are the zeros.
for col in ("PFS_event", "OS_event"):
    samples[col] = pd.to_numeric(samples[col], errors="coerce")

expr.to_csv(a.out_expr, sep="\t")
samples.to_csv(a.out_samples, sep="\t", index=False)
print(samples["timepoint"].value_counts().to_string(), file=sys.stderr)
print(samples["response3"].value_counts(dropna=False).to_string(), file=sys.stderr)
print(f"wrote {a.out_expr} and {a.out_samples}", file=sys.stderr)
