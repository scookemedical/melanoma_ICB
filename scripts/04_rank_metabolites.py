"""Rank panel metabolites by predicted impact of glutamine supplementation.

Impact score has three interpretable, separable terms:

  1. response   - fractional change in production capacity across the dose
                  sweep, estimated per sample from the GEM (data-driven)
  2. sensitivity- kinetic responsiveness of the metabolite's effector,
                  Km / (Km + C): the elasticity of fractional effector
                  occupancy with respect to a fractional concentration change.
                  Near 1 when the effector is far from saturation, near 0 when
                  it is saturated and a flux change therefore buys nothing.
  3. coupling   - flux-to-concentration coupling x literature evidence weight
                  (hand-curated; see config/metabolite_panel.tsv)

Terms 2 and 3 are priors, so the flux-only ranking is reported alongside the
weighted one and a weight-perturbation analysis reports how stable rank 1 is.
"""
import argparse
import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--sweep", required=True)
ap.add_argument("--panel", required=True)
ap.add_argument("--out-ranking", required=True)
ap.add_argument("--out-stability", required=True)
ap.add_argument("--dose-lo", type=float, default=1.0)
ap.add_argument("--dose-hi", type=float, default=2.0)
ap.add_argument("--boot", type=int, default=2000)
ap.add_argument("--seed", type=int, default=0)
a = ap.parse_args()

rng = np.random.default_rng(a.seed)
sweep = pd.read_csv(a.sweep, sep="\t")
panel = pd.read_csv(a.panel, sep="\t")

# ---- term 1: per-sample fractional response over the dose range ----------
lo = sweep[np.isclose(sweep["dose"], a.dose_lo)].set_index(["sample", "met_id"])["flux"]
hi = sweep[np.isclose(sweep["dose"], a.dose_hi)].set_index(["sample", "met_id"])["flux"]
pair = pd.concat([lo.rename("lo"), hi.rename("hi")], axis=1).dropna()
pair = pair[pair["lo"] > 1e-9]
pair["frac"] = (pair["hi"] - pair["lo"]) / pair["lo"]

# dose-response slope across every dose, as a shape check on linearity
slope = (sweep.dropna(subset=["flux"])
         .groupby(["sample", "met_id"])
         .apply(lambda g: np.polyfit(g["dose"], g["flux"], 1)[0]
                if g["dose"].nunique() > 1 else np.nan,
                include_groups=False)
         .rename("slope"))

rows = []
for met_id, grp in pair.groupby(level="met_id"):
    v = grp["frac"].to_numpy()
    bs = rng.choice(v, size=(a.boot, v.size), replace=True).mean(axis=1)
    rows.append({
        "met_id": met_id,
        "n_samples": v.size,
        "response_mean": v.mean(),
        "response_sd": v.std(ddof=1) if v.size > 1 else np.nan,
        "response_lo95": np.percentile(bs, 2.5),
        "response_hi95": np.percentile(bs, 97.5),
        "baseline_flux": grp["lo"].mean(),
    })
resp = pd.DataFrame(rows)

df = panel.merge(resp, on="met_id", how="inner")

# ---- term 2: kinetic sensitivity -----------------------------------------
df["sensitivity"] = df["km_uM"] / (df["km_uM"] + df["conc_uM"])

# ---- composite -----------------------------------------------------------
df["prior"] = df["pool_coupling"] * df["evidence_weight"]
df["impact_score"] = df["response_mean"].abs() * df["sensitivity"] * df["prior"]
df["signed_impact"] = df["impact_score"] * df["direction"]

df["rank_flux_only"] = df["response_mean"].abs().rank(ascending=False).astype(int)
df["rank_weighted"] = df["impact_score"].rank(ascending=False).astype(int)
df = df.sort_values("impact_score", ascending=False)

cols = ["rank_weighted", "rank_flux_only", "label", "met_id", "pathway",
        "baseline_flux", "response_mean", "response_lo95", "response_hi95",
        "sensitivity", "pool_coupling", "evidence_weight", "impact_score",
        "direction", "signed_impact", "n_samples", "effector", "km_uM",
        "conc_uM", "citation"]
df[cols].to_csv(a.out_ranking, sep="\t", index=False)

# ---- weight-perturbation stability --------------------------------------
# Re-rank under multiplicative jitter on the two prior terms (log-uniform,
# +/- 50%) and on the response estimate (its own bootstrap SE).
n = len(df)
top_counts = pd.Series(0, index=df["met_id"])
top3_counts = pd.Series(0, index=df["met_id"])
se = (df["response_hi95"] - df["response_lo95"]) / 3.92
for _ in range(a.boot):
    jitter_sens = np.exp(rng.uniform(np.log(0.5), np.log(1.5), n))
    jitter_prior = np.exp(rng.uniform(np.log(0.5), np.log(1.5), n))
    resp_draw = rng.normal(df["response_mean"].to_numpy(),
                           np.nan_to_num(se.to_numpy(), nan=0.0))
    score = (np.abs(resp_draw)
             * (df["sensitivity"].to_numpy() * jitter_sens)
             * (df["prior"].to_numpy() * jitter_prior))
    order = np.argsort(-score)
    top_counts.iloc[order[0]] += 1
    for i in order[:3]:
        top3_counts.iloc[i] += 1

stab = pd.DataFrame({
    "met_id": df["met_id"].to_numpy(),
    "label": df["label"].to_numpy(),
    "p_rank1": (top_counts.to_numpy() / a.boot),
    "p_top3": (top3_counts.to_numpy() / a.boot),
}).sort_values("p_rank1", ascending=False)
stab.to_csv(a.out_stability, sep="\t", index=False)

pd.set_option("display.width", 200)
print(df[["rank_weighted", "label", "response_mean", "sensitivity",
          "prior", "impact_score"]].round(4).to_string(index=False))
print()
print(stab.round(3).to_string(index=False))
