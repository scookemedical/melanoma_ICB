# Results

Riaz anti-PD-1 melanoma cohort, 51 pre-treatment tumours, Human-GEM v2.1.1,
glutamine uptake swept 1.0-2.5x the plasma reference (0.6 mM), `vmax = 0.3`.
Full tables in `results/`.

## 1. Ranking

`response` is **yield**: mol of extra metabolite turnover per mol of extra
glutamine taken up, averaged over samples. Negative means glutamine load
*reduces* that metabolite's synthesis flux.

| Rank | Metabolite | Yield (mol/mol) | 95% CI | Sensitivity | Prior | Impact | P(rank 1) |
|---|---|---|---|---|---|---|---|
| 1 | ornithine | +0.337 | 0.30 to 0.38 | 0.77 | 0.24 | 0.062 | 0.55 |
| 2 | ammonia | +0.256 | 0.19 to 0.33 | 0.71 | 0.30 | 0.055 | 0.40 |
| 3 | alpha-ketoglutarate | +0.103 | 0.01 to 0.19 | 0.27 | 0.86 | 0.024 | 0.04 |
| 4 | fumarate | −0.098 | −0.13 to −0.06 | 0.55 | 0.30 | 0.016 | 0.005 |
| 5 | succinate | +0.060 | 0.02 to 0.11 | 0.33 | 0.30 | 0.006 | 0 |
| 6 | proline | −0.125 | −0.16 to −0.09 | 0.23 | 0.20 | 0.006 | 0 |
| 7 | citrate | −0.086 | **−0.19 to +0.01** | 0.14 | 0.36 | 0.004 | 0 |
| 8 | glutathione | −0.023 | −0.05 to −0.002 | 0.38 | 0.24 | 0.002 | 0 |
| 9-13 | aspartate, glutamate, asparagine, UDP-GlcNAc, spermidine | <0.05 | | | | ~0 | 0 |

(S)-2-hydroxyglutarate and putrescine carry **no flux at any dose** in the
parsimonious solution and cannot be ranked.

**Rank 1 is not robust, and rank 1 is not alpha-ketoglutarate.** Under +/-50%
jitter on the prior terms, ornithine holds rank 1 in 55% of draws and ammonia in
40% — the two are not separable. Alpha-ketoglutarate, the a priori favourite and
the only panel member with direct in vivo melanoma evidence, lands third
(P(rank 1) = 0.04) because its effector dioxygenases are closer to saturation
(sensitivity 0.27) than ornithine's or ammonia's.

**Ornithine's rank 1 rests on a threshold, not a gradient.** Its baseline
turnover is zero in 38 of 41 evaluable tumours: ornithine synthesis is *off* at
reference glutamine and *switches on* under load. That is an interesting
qualitative prediction — a pathway recruited only by glutamine excess, feeding
the arginase/polyamine immunosuppressive axis — but the yield estimate is
dominated by an on/off transition, its Cox model is undefined for lack of
variance, and it should be treated as the weakest-supported entry in the table
despite its position.

The biological shape of the result is coherent: glutamine load is routed to
nitrogen disposal (ornithine, ammonia) and to alpha-ketoglutarate, while
*reducing* flux to the carbon sinks fed by anaplerosis (citrate, fumarate,
proline). That matches the in vivo observation that dietary glutamine raised
intratumoural alpha-ketoglutarate without raising biosynthetic intermediates.

## 2. Falsification against ICB outcome — the ranking does not survive it

Baseline (dose 1.0) turnover vs outcome in the same tumours; 9 responders
(PR/CR) vs 34 non-responders (SD/PD), FDR-corrected across the panel.

| Metabolite | Impact rank | AUC (response) | q | PFS HR per SD | q |
|---|---|---|---|---|---|
| **citrate** | 7 | 0.26 | 0.43 | **1.80** | **0.0024** |
| **glutathione** | 8 | 0.46 | 0.89 | **1.57** | **0.038** |
| succinate | 5 | 0.38 | 0.89 | 1.39 | 0.11 |
| UDP-GlcNAc | 12 | 0.52 | 0.89 | 1.44 | 0.11 |
| ornithine | 1 | 0.53 | 0.89 | not estimable | — |
| ammonia | 2 | 0.53 | 0.89 | 1.00 | 0.99 |
| alpha-ketoglutarate | 3 | 0.46 | 0.89 | 1.13 | 0.53 |

- **No metabolite's predicted turnover separates responders from
  non-responders** (every AUC 0.26-0.60, every q > 0.4).
- **Two associate with PFS after FDR correction, and both are ranked 7th and
  8th:** higher predicted citrate synthesis flux (HR 1.80 per SD, q = 0.002) and
  higher predicted glutathione synthesis flux (HR 1.57, q = 0.038) each predict
  *shorter* progression-free survival on nivolumab.
- **Nothing associates with overall survival.**
- The top three ranked metabolites show no outcome association at all.

So the impact ranking and the only outcome-linked signals disagree. The ranking
should not be used to pick an experiment on its own.

## 3. What is actually worth testing

The one hypothesis that came from the data rather than from the priors:
**tumours with a high predicted citrate and glutathione synthesis flux — a
lipogenic, redox-buffered phenotype — progress faster on anti-PD-1.** Citrate
and glutathione are the only two metabolites that are both outcome-associated
and predicted to move in the protective direction under glutamine load.

Two caveats specific to this hypothesis, and they are not small:

- **Citrate's predicted yield is not significantly different from zero**
  (−0.086, 95% CI −0.185 to +0.007). Its outcome association is solid; the claim
  that glutamine *lowers* it is not. Glutathione's yield does exclude zero
  (−0.023, CI −0.052 to −0.002) but is an order of magnitude smaller.
- A single cohort, 9 responders, two surviving hits out of 39 tests, and bulk
  tumour rather than compartment-resolved data.

So the defensible version is the prognostic half — high citrate/glutathione
synthesis flux marks faster progression — and the supplementation half is
unproven. This is a hypothesis to test, not a biomarker.

## 4. Numerical caveats that bound every number above

- **9.4% of sample-dose pairs (24/255) are infeasible** and drop out. Pinning
  glutamine uptake exactly gave 16% dropout rising monotonically with dose,
  which would have biased yield toward tumours that absorb glutamine well; a 10%
  uptake slack removed the dose trend (8/7/3/2/4 failures across the five
  doses). The residual failures are per-sample, not dose-driven.
- **32 of 51 tumours solve at all five doses**; the yield estimate uses the 41
  with both endpoint doses.
- **Dose is a multiple of a plasma reference uptake bound, not a supplement
  dose.** Nothing here converts to g/day.
- **Evidence weights and pool-coupling terms are provisional** and should be
  re-curated before use; the flux-only ranking (`rank_flux_only`) is reported so
  they can be bypassed entirely.
