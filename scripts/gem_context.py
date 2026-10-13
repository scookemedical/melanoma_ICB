"""Shared helpers: defined medium, E-Flux contextualisation, production capacity.

E-Flux (Colijn et al., PLoS Comput Biol 2009) scales reaction flux bounds by
GPR-resolved transcript abundance. It needs one LP per condition, which is what
makes a per-sample x per-dose sweep tractable; iMAT/tINIT would require a MILP
per sample.
"""
import ast
import cobra

# ---------------------------------------------------------------- medium ----
# Uptake bound = plasma concentration (mM) x PLASMA_SCALE (mmol/gDW/h per mM).
# Concentrations are adult venous plasma reference values.
PLASMA_SCALE = 0.2

PLASMA_mM = {
    "MAR09034": 5.0,    # glucose
    "MAR09063": 0.60,   # glutamine   <-- the swept input
    "MAR09061": 0.33,   # alanine
    "MAR09066": 0.08,   # arginine
    "MAR09062": 0.05,   # asparagine
    "MAR09070": 0.01,   # aspartate
    "MAR09065": 0.05,   # cysteine
    "MAR09071": 0.06,   # glutamate
    "MAR09067": 0.25,   # glycine
    "MAR09038": 0.08,   # histidine
    "MAR09039": 0.06,   # isoleucine
    "MAR09040": 0.12,   # leucine
    "MAR09041": 0.18,   # lysine
    "MAR09042": 0.03,   # methionine
    "MAR09043": 0.06,   # phenylalanine
    "MAR09068": 0.17,   # proline
    "MAR09069": 0.12,   # serine
    "MAR09044": 0.15,   # threonine
    "MAR09045": 0.06,   # tryptophan
    "MAR09064": 0.06,   # tyrosine
    "MAR09046": 0.22,   # valine
    "MAR09035": 0.10,   # linoleate
    "MAR09036": 0.05,   # linolenate
    "MAR09083": 0.01,   # choline
    "MAR09146": 0.00002,  # folate
    "MAR09159": 0.0001,   # thiamin
    "MAR09143": 0.0003,   # riboflavin
    "MAR09378": 0.0005,   # nicotinamide
    "MAR09144": 0.00005,  # pyridoxine
    "MAR09145": 0.001,    # pantothenate
    "MAR09109": 0.000002, # biotin
    "MAR09147": 0.002,    # retinol
    "MAR09151": 0.025,    # alpha-tocopherol
    "MAR09169": 0.001,    # selenate
    "MAR09358": 0.01,     # hypoxanthine
    "MAR09362": 0.01,     # inosine
    "MAR09153": 0.005,    # gamma-tocopherol
    "MAR09167": 0.00005,  # lipoic acid
    "MAR09269": 0.0003,   # aquacob(III)alamin (vitamin B12)
}

# Freely available (ions, gases, water): uptake not limiting in vivo.
UNLIMITED_UPTAKE = [
    "MAR09047",  # H2O
    "MAR09048",  # O2
    "MAR09079",  # H+
    "MAR09078",  # HCO3-
    "MAR09072",  # Pi
    "MAR09074",  # sulfate
    "MAR09077",  # Na+
    "MAR09081",  # K+
    "MAR09082",  # Ca2+
    "MAR13072",  # Mg2+
    "MAR09076",  # Fe2+
    "MAR13073",  # Cu2+
]

GLN_EXCHANGE = "MAR09063"
BIOMASS = "MAR13082"


def apply_medium(model, gln_dose=1.0, scale=PLASMA_SCALE, fix_gln=True,
                 gln_slack=0.10):
    """Close all uptake, then open the defined medium.

    gln_dose multiplies the glutamine uptake bound only; 1.0 = plasma reference
    (0.6 mM). Secretion stays free (ub = 1000) for every boundary reaction.

    With `fix_gln` the glutamine exchange is *fixed* at the dose rather than
    merely capped by it. This matters more than it looks. Glutamine uptake is
    concentration-driven, so raising plasma glutamine raises uptake whether or
    not the cell needs the carbon; left as an upper bound, a growth-maximising
    model simply declines the extra glutamine and the dose response is zero.
    Fixing the flux forces the model to dispose of the supplemented load and
    asks the question actually posed -- given this much extra glutamine, where
    does it go -- which is also what was observed in vivo: dietary glutamine
    raised intratumoural alpha-ketoglutarate without raising biosynthetic
    intermediates (Ishak Gabra et al., Nat Commun 2020).
    """
    for rxn in model.reactions:
        if rxn.boundary:
            rxn.lower_bound = 0.0
            rxn.upper_bound = 1000.0
    for rid in UNLIMITED_UPTAKE:
        model.reactions.get_by_id(rid).lower_bound = -1000.0
    for rid, conc in PLASMA_mM.items():
        bound = conc * scale
        if rid == GLN_EXCHANGE:
            bound *= gln_dose
        model.reactions.get_by_id(rid).lower_bound = -bound
    if fix_gln:
        # Pinned exactly, a high dose plus tight expression caps leaves some
        # samples with no feasible way to dispose of the load: 16% of
        # sample-dose pairs went infeasible, and the rate rose with dose, which
        # would have biased the yield estimate toward tumours that happen to
        # absorb glutamine well. `gln_slack` gives the model a little room
        # (uptake between (1 - slack) and 1 x the dose), which is also the more
        # honest statement: uptake is concentration-driven, not clamped.
        gln = model.reactions.get_by_id(GLN_EXCHANGE)
        gln.upper_bound = gln.lower_bound * (1.0 - gln_slack)
    return model


# ------------------------------------------------------------ E-Flux ----
def _eval_gpr(node, expr, default, or_op=max):
    if isinstance(node, ast.Expression):
        return _eval_gpr(node.body, expr, default, or_op)
    if isinstance(node, ast.Name):
        return expr.get(node.id, default)
    if isinstance(node, ast.BoolOp):
        vals = [_eval_gpr(v, expr, default, or_op) for v in node.values]
        if not vals:
            return default
        # AND = min: a complex is limited by its scarcest subunit.
        # OR = max by default: isozymes are alternatives, and summing them
        # inflates many-isozyme reactions by an order of magnitude, which then
        # dominates the bound normalisation. Pass or_op=sum for classic E-Flux.
        return min(vals) if isinstance(node.op, ast.And) else or_op(vals)
    raise ValueError(f"unexpected GPR node: {ast.dump(node)}")


def reaction_scores(model, gene_expr, default=None, or_op=max):
    """GPR-resolve gene expression to a per-reaction activity score.

    Genes absent from the matrix fall back to `default` (median of measured
    model genes) rather than 0, so missing measurement is not read as absence.
    """
    if default is None:
        vals = sorted(gene_expr.values())
        default = vals[len(vals) // 2] if vals else 0.0
    scores = {}
    for rxn in model.reactions:
        rule = rxn.gene_reaction_rule.strip()
        if not rule:
            scores[rxn.id] = None          # no GPR -> leave bounds untouched
            continue
        try:
            tree = ast.parse(rule.replace(" and ", " and ").replace(" or ", " or "),
                             mode="eval")
            scores[rxn.id] = float(_eval_gpr(tree, gene_expr, default, or_op))
        except Exception:
            scores[rxn.id] = default
    return scores


def apply_eflux(model, scores, vmax=10.0, norm_pct=99.0, floor_frac=0.01):
    """Scale each GPR-backed reaction's bounds to its expression score.

    Bounds are normalised to the `norm_pct` percentile rather than the maximum
    so that one outlier reaction cannot compress every other bound toward zero;
    ratios above 1 are clipped. `vmax` must be of the same order as the medium
    uptake fluxes (~1 mmol/gDW/h here) or the bounds never bind and the result
    collapses back to a stoichiometric yield calculation.

    `floor_frac` keeps every bound at >= floor_frac * vmax. Without it a single
    unexpressed subunit inside an AND clause drives the cap to exactly zero and
    can knock out a reaction the biomass equation needs, leaving the model
    unable to grow at all: a transcript below detection means low, not absent.
    """
    import numpy as np
    finite = [s for s in scores.values() if s is not None]
    top = float(np.percentile(finite, norm_pct)) if finite else 1.0
    if top <= 0:
        top = 1.0
    for rxn in model.reactions:
        s = scores.get(rxn.id)
        if s is None or rxn.boundary:
            continue
        cap = vmax * max(floor_frac, min(1.0, s / top))
        if rxn.upper_bound > 0:
            rxn.upper_bound = min(rxn.upper_bound, cap)
        if rxn.lower_bound < 0:
            rxn.lower_bound = max(rxn.lower_bound, -cap)
    return model


# ------------------------------------------------------ metabolite turnover ----
def metabolite_turnover(model, met_ids, fraction_of_optimum=1.0):
    """Gross synthesis flux of each metabolite in the pFBA flux distribution.

    Returns ({met_id: turnover}, biomass).

    Why turnover and not a maximised demand reaction: maximising a demand
    reaction one metabolite at a time returns a *capacity envelope*, and in a
    network with ~12,700 reactions the optimum simply reroutes around any
    expression cap. Tested directly on this model, per-sample capacities varied
    by <1% while sitting at exact stoichiometric ratios of each other -- i.e.
    the readout had collapsed back onto stoichiometry and uptake bounds, which
    is precisely the artefact this analysis exists to avoid.

    Net production is identically zero for any internal metabolite at steady
    state, so the informative quantity is gross synthesis: the sum of positive
    production terms across the metabolite's reactions,

        turnover_i = sum_r max(0, S_ir * v_r)

    read off a parsimonious FBA solution (maximise biomass, then minimise total
    absolute flux). pFBA gives one well-defined flux distribution that has to
    respect every expression cap, so the result is transcriptome-dependent --
    the same readout style used by METAFlux and scFEA.
    """
    from cobra.flux_analysis import pfba
    bio = model.reactions.get_by_id(BIOMASS)
    sol = model.optimize()
    if sol.status != "optimal" or sol.objective_value is None:
        return {m: float("nan") for m in met_ids}, float("nan")
    bio_opt = sol.objective_value
    try:
        flux = pfba(model, fraction_of_optimum=fraction_of_optimum).fluxes
    except Exception:
        return {m: float("nan") for m in met_ids}, bio_opt
    out = {}
    for met_id in met_ids:
        met = model.metabolites.get_by_id(met_id)
        total = 0.0
        for rxn in met.reactions:
            made = rxn.metabolites[met] * flux.get(rxn.id, 0.0)
            if made > 0:
                total += made
        out[met_id] = total
    return out, bio_opt
