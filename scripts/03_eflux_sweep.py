"""Glutamine dose sweep: per-sample E-Flux model -> metabolite production capacity.

For every tumour sample the transcriptome caps reaction fluxes (E-Flux); the
glutamine exchange bound is then stepped across a dose range and the maximum
net production of each panel metabolite is recorded while biomass is held near
its own optimum.
"""
import argparse, os, pickle, sys, time
from multiprocessing import Pool

import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gem_context as gc

MODEL = None
PANEL = None
DOSES = None
VMAX = 0.3
FLOOR = 0.005


def _init(model_path, met_ids, doses, vmax, floor):
    global MODEL, PANEL, DOSES, VMAX, FLOOR
    with open(model_path, "rb") as fh:
        MODEL = pickle.load(fh)
    MODEL.solver = "glpk"
    PANEL, DOSES, VMAX, FLOOR = met_ids, doses, vmax, floor


def _one_sample(args):
    sample, expr = args
    t0 = time.time()
    scores = gc.reaction_scores(MODEL, expr)
    rows = []
    for dose in DOSES:
        with MODEL:
            gc.apply_medium(MODEL, gln_dose=dose)
            gc.apply_eflux(MODEL, scores, vmax=VMAX, floor_frac=FLOOR)
            fluxes, bio = gc.metabolite_turnover(MODEL, PANEL)
        for met, val in fluxes.items():
            rows.append({"sample": sample, "dose": dose, "met_id": met,
                         "turnover": val, "biomass": bio})
    print(f"  {sample} done in {time.time()-t0:.0f}s", file=sys.stderr, flush=True)
    return rows


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--model", required=True)
    ap.add_argument("--expr", required=True)
    ap.add_argument("--samples", required=True)
    ap.add_argument("--panel", required=True)
    ap.add_argument("--out", required=True)
    ap.add_argument("--timepoint", default="Pre",
                    help="Pre, On, or all")
    ap.add_argument("--doses", default="1.0,1.25,1.5,1.75,2.0,2.5")
    ap.add_argument("--vmax", type=float, default=0.3,
                    help="E-Flux bound scale. Must be of the same order as the "
                         "realised fluxes, not the uptake bounds: at 0.3 the "
                         "bounds bind and biomass itself becomes "
                         "transcriptome-limited; at 1 or above nothing binds "
                         "and the readout collapses to stoichiometry; at 0.1 "
                         "the model is infeasible.")
    ap.add_argument("--floor-frac", type=float, default=0.005,
                    help="minimum bound as a fraction of vmax")
    ap.add_argument("--procs", type=int, default=4)
    ap.add_argument("--limit", type=int, default=0, help="smoke-test subset")
    a = ap.parse_args()

    doses = [float(x) for x in a.doses.split(",")]
    panel = pd.read_csv(a.panel, sep="\t")
    met_ids = panel["met_id"].tolist()

    expr = pd.read_csv(a.expr, sep="\t", index_col=0)
    samples = pd.read_csv(a.samples, sep="\t")
    keep = samples if a.timepoint == "all" else samples[samples["timepoint"] == a.timepoint]
    keep = keep[~keep["replicate"].astype(bool)]
    cols = [c for c in keep["sample"] if c in expr.columns]
    if a.limit:
        cols = cols[:a.limit]
    print(f"{len(cols)} samples x {len(doses)} doses x {len(met_ids)} metabolites",
          file=sys.stderr)

    jobs = [(c, expr[c].to_dict()) for c in cols]
    with Pool(a.procs, initializer=_init,
              initargs=(a.model, met_ids, doses, a.vmax, a.floor_frac)) as pool:
        out = []
        for rows in pool.imap_unordered(_one_sample, jobs):
            out.extend(rows)
    df = pd.DataFrame(out)
    df.to_csv(a.out, sep="\t", index=False)
    print(f"wrote {a.out}  ({len(df)} rows)", file=sys.stderr)


if __name__ == "__main__":
    main()
