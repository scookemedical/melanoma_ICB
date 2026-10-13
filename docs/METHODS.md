# Ranking glutamine-derived metabolites by predicted pathway impact

## What this answers, and what it does not

The question posed was: given a set amount of glutamine supplementation, rank
glutamine-derived metabolites by the impact of their stoichiometrically produced
amounts on cancer-promoting pathways, and name the most impactful one.

**The stoichiometric form of that question is not answerable, and this pipeline
deliberately does not answer it.** Glutamine catabolism is a branch point, not a
linear pathway: the same carbon and nitrogen partition between alpha-ketoglutarate
and the TCA cycle, glutathione, the hexosamine pathway, aspartate and nucleotides,
proline, asparagine, polyamines, and free ammonia. The split ratios are set by
enzyme abundance, redox state and nutrient availability, not by fixed coefficients.
Three specific failure points:

1. **Dose does not map onto intratumoural concentration.** Splanchnic first-pass
   extraction consumes most of an oral glutamine load, ASCT2/SLC1A5 transport
   saturates, and glutamine is steeply graded between tumour core and periphery
   (Pan et al., Nat Cell Biol 2016).
2. **The direct evidence in this disease contradicts fixed stoichiometry.** In
   transgenic melanoma, dietary glutamine raised intratumoural glutamine *and*
   alpha-ketoglutarate **without** raising biosynthetic intermediates
   (Ishak Gabra et al., Nat Commun 2020, PMID 32620791). Flux partitioned
   selectively; a stoichiometric model predicts proportional increases across
   every branch and would have been wrong.
3. **Amount produced is not impact.** Alpha-ketoglutarate matters because the Km
   of its effector dioxygenases sits near its physiological concentration.
   Glutathione is buffered in the millimolar range, so a large fractional flux
   change moves effector occupancy very little.

So the quantity swept here is a **glutamine uptake bound**, and the quantity
ranked is **predicted production capacity under a transcriptome-constrained
genome-scale model**, weighted by effector kinetics and literature evidence. The
output is a prioritised hypothesis list for experimental testing, not a
quantitative prediction of pathway magnitude.

## Inputs

| Input | Source | Notes |
|---|---|---|
| Human-GEM v2.1.1 (12,748 reactions, 8,374 metabolites, 2,870 genes) | `SysBioChalmers/Human-GEM` | genome-scale human reconstruction |
| Riaz et al., Cell 2017 (GSE91061) bulk RNA-seq, 118 samples / 56 pre- and 62 on-treatment | `riazn/bms038_analysis` (the authors' own analysis repo) | nivolumab-treated metastatic melanoma, with BOR, PFS and OS |

2,825 of the model's 2,870 genes (98.4%) are measured in the Riaz matrix.

### Why these sources and not GEO

This session's egress policy blocks `ftp.ncbi.nlm.nih.gov`, `www.ncbi.nlm.nih.gov`,
`ebi.ac.uk`, `zenodo.org`, `cbioportal.org`, `xenabrowser.net` and
`api.gdc.cancer.gov`; only GitHub and PyPI are reachable. The Riaz cohort was
obtained from the study authors' own GitHub analysis repository, which carries the
same count matrix and clinical table as the GEO record. Adding further cohorts
(Hugo GSE78220, Gide PRJEB23709, Liu, Van Allen) requires widening the network
policy for the environment.

## Pipeline

```
01_fetch_data.sh        acquire Human-GEM + Riaz
02_prepare_expression.py  counts -> log2 CPM, Entrez -> Ensembl, Human-GEM genes
03_eflux_sweep.py       per-sample E-Flux model; glutamine dose sweep;
                        max production capacity per panel metabolite
04_rank_metabolites.py  impact score, bootstrap CIs, weight-stability analysis
05_validate_icb.py      association with response / PFS / OS in Riaz
```

### Defined medium

Human-GEM ships with all 1,654 boundary reactions fully open, which supports an
unconstrained biomass flux of ~138 and is biologically meaningless. Uptake is
therefore closed and reopened as a plasma-referenced medium: uptake bound =
plasma concentration (mM) x 0.2 mmol/gDW/h per mM, with water, O2, H+, phosphate,
sulfate and ions left unlimited and secretion free. Vitamin B12
(aquacob(III)alamin), lipoic acid and gamma-tocopherol must be included or the
cofactor biomass pool cannot be produced and the model cannot grow at all.
Glutamine's reference is 0.6 mM; the sweep multiplies that bound only.

### Transcriptome constraint (E-Flux)

E-Flux (Colijn et al., PLoS Comput Biol 2009) scales reaction flux bounds by
GPR-resolved transcript abundance. It needs one LP per condition, which is what
makes a per-sample x per-dose sweep tractable; iMAT or tINIT would need a MILP
per sample. Three implementation choices matter:

- **OR = max, AND = min.** Classic E-Flux sums isozymes, which inflates
  many-isozyme reactions by an order of magnitude (observed 99th percentile 8.5
  under max vs several hundred under sum) and then compresses every other bound
  during normalisation.
- **Normalisation to the 99th percentile, clipped at 1,** for the same reason.
- **`vmax` on the order of the uptake fluxes, and a 1% bound floor.** If `vmax`
  is far above the uptake scale the bounds never bind and the result collapses
  back to a stoichiometric yield calculation — exactly the artefact this analysis
  exists to avoid. The floor is needed because one unexpressed subunit in an AND
  clause otherwise drives a bound to exactly zero and can make the model unable
  to grow; a transcript below detection means low, not absent.

### Production capacity

For each metabolite a temporary demand reaction (`met ->`) is added and maximised
with biomass floored at 90% of its own optimum, so the number is a capacity under
a still-proliferating phenotype rather than a free-for-all conversion yield.

### Impact score

```
impact = |fractional response| x sensitivity x pool_coupling x evidence_weight
sensitivity = Km / (Km + C)
```

`sensitivity` is the elasticity of fractional effector occupancy with respect to a
fractional concentration change: near 1 when the effector is far from saturation,
near 0 when it is saturated and extra flux buys nothing. This is the term that
separates "more is made" from "more matters". `pool_coupling` and
`evidence_weight` are hand-curated priors in `config/metabolite_panel.tsv`.

Because two of the four terms are priors, the ranking could be self-fulfilling.
Two guards are built in:

1. The **flux-only ranking** is reported alongside the weighted one.
2. A **weight-perturbation analysis** re-ranks under +/-50% log-uniform jitter on
   both prior terms and a normal draw on the response estimate, reporting the
   probability that each metabolite holds rank 1.

### Falsification

`05_validate_icb.py` tests every panel metabolite's baseline production capacity
against response (PRCR vs SDPD), PFS and OS in the Riaz cohort, FDR-corrected
across the panel. The priors cannot influence this step, so it is the one place
the ranking can be shown to be wrong.

## Known limitations

- **One cohort.** No external replication until the network policy allows a
  second dataset.
- **Bulk, not compartment-resolved.** Glutamine supplementation suppresses
  melanoma cell-autonomously, while glutamine *antagonism* enhances CD8 T cell
  antitumour immunity (Leone et al., Science 2019). Bulk tumour averages these
  opposing compartments. Deconvolution or scRNA-seq is required before any
  claim about net clinical direction, and three panel entries (asparagine,
  ammonia, spermidine) have compartment-dependent direction for this reason.
- **E-Flux gives relative, not absolute, flux.** Dose is a multiple of the plasma
  reference uptake bound; it is not convertible to a supplement dose in g/day.
- **Evidence weights are provisional** and should be re-curated by the user
  before the ranking is used to choose an experiment.
- **No TCGA-SKCM or healthy-tissue baseline**, so the medium is a plasma
  reference rather than a measured tumour interstitial composition.
