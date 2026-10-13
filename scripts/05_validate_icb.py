"""Falsification step: does predicted production capacity track ICB outcome?

Any ranking built partly on literature priors can be self-fulfilling, so the
ranked metabolites are tested against the one thing the priors cannot touch:
response and survival on anti-PD-1 in the Riaz cohort. Every panel metabolite
is tested, not only the top-ranked one, with FDR correction across the panel.
"""
import argparse
import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests
from statsmodels.duration.hazard_regression import PHReg

ap = argparse.ArgumentParser()
ap.add_argument("--sweep", required=True)
ap.add_argument("--samples", required=True)
ap.add_argument("--panel", required=True)
ap.add_argument("--ranking", required=True)
ap.add_argument("--out", required=True)
ap.add_argument("--dose", type=float, default=1.0,
                help="dose at which to take the baseline capacity")
a = ap.parse_args()

sweep = pd.read_csv(a.sweep, sep="\t")
samples = pd.read_csv(a.samples, sep="\t")
panel = pd.read_csv(a.panel, sep="\t")
lab = dict(zip(panel["met_id"], panel["label"]))

base = sweep[np.isclose(sweep["dose"], a.dose)]
wide = base.pivot_table(index="sample", columns="met_id", values="turnover")
meta = samples.set_index("sample")
dat = wide.join(meta[["response3", "PFS", "PFS_event", "OS", "OS_event",
                      "cohort", "patient"]], how="inner")

# Riaz codes *_SOR as 1 = censored, so the event indicator is the complement.
dat["pfs_evt"] = 1 - dat["PFS_event"]
dat["os_evt"] = 1 - dat["OS_event"]

rows = []
for met in wide.columns:
    x = dat[met]
    if x.notna().sum() < 10 or x.std(ddof=1) < 1e-12:
        rows.append({"met_id": met, "label": lab.get(met, met),
                     "note": "no variance across samples"})
        continue
    sub = dat[dat["response3"].isin(["PRCR", "SDPD"])]
    r = sub[sub["response3"] == "PRCR"][met].dropna()
    nr = sub[sub["response3"] == "SDPD"][met].dropna()
    rec = {"met_id": met, "label": lab.get(met, met),
           "n_resp": len(r), "n_nonresp": len(nr),
           "median_resp": r.median(), "median_nonresp": nr.median()}
    if len(r) >= 3 and len(nr) >= 3:
        u, p = stats.mannwhitneyu(r, nr, alternative="two-sided")
        rec["auc"] = u / (len(r) * len(nr))
        rec["p_response"] = p
    # Cox on z-scored capacity, so the HR is per SD
    for tag, tcol, ecol in (("pfs", "PFS", "pfs_evt"), ("os", "OS", "os_evt")):
        cox = dat[[met, tcol, ecol]].dropna()
        cox = cox[cox[tcol] > 0]
        # A predictor that is near-constant makes the partial likelihood flat:
        # the fit "succeeds" with an absurd HR and p = 1. Require real spread.
        if (len(cox) >= 15 and cox[ecol].sum() >= 5
                and cox[met].nunique() >= 5):
            cap = cox[met].to_numpy(dtype=float)
            cap = (cap - cap.mean()) / cap.std(ddof=1)
            try:
                fit = PHReg(cox[tcol].to_numpy(dtype=float),
                            cap.reshape(-1, 1),
                            status=cox[ecol].to_numpy(dtype=float)).fit()
                if not np.isfinite(fit.params[0]) or abs(fit.params[0]) > 10:
                    raise ValueError("Cox fit did not converge to a sane HR")
                rec[f"hr_{tag}"] = float(np.exp(fit.params[0]))
                rec[f"p_{tag}"] = float(fit.pvalues[0])
            except Exception as exc:
                rec[f"p_{tag}"] = np.nan
                rec[f"note_{tag}"] = str(exc)[:60]
        else:
            rec[f"note_{tag}"] = "insufficient events or near-constant predictor"
    rows.append(rec)

res = pd.DataFrame(rows)
for col in ("p_response", "p_pfs", "p_os"):
    if col in res and res[col].notna().any():
        ok = res[col].notna()
        res.loc[ok, col.replace("p_", "q_")] = multipletests(
            res.loc[ok, col], method="fdr_bh")[1]

rank = pd.read_csv(a.ranking, sep="\t")[["met_id", "rank_weighted", "rank_flux_only"]]
res = res.merge(rank, on="met_id", how="left").sort_values("rank_weighted")
res.to_csv(a.out, sep="\t", index=False)

pd.set_option("display.width", 200)
show = [c for c in ["rank_weighted", "label", "n_resp", "n_nonresp", "auc",
                    "p_response", "q_response", "hr_pfs", "p_pfs", "q_pfs",
                    "hr_os", "p_os", "q_os"] if c in res]
print(res[show].round(4).to_string(index=False))
