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


def apply_medium(model, gln_dose=1.0, scale=PLASMA_SCALE):
    """Close all uptake, then open the defined medium.

    gln_dose multiplies the glutamine uptake bound only; 1.0 = plasma reference
    (0.6 mM). Secretion stays free (ub = 1000) for every boundary reaction.
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


# -------------------------------------------------- production capacity ----
def production_capacity(model, met_ids, biomass_frac=0.9):
    """Max net production of each metabolite while biomass stays near optimum.

    Returns {met_id: flux}. A temporary demand reaction (met -> sink) is added
    and maximised; biomass is floored at `biomass_frac` x its own optimum so the
    cell is still required to grow, which is what makes the number a capacity
    under a proliferating phenotype rather than a free-for-all.
    """
    out = {}
    bio = model.reactions.get_by_id(BIOMASS)
    sol = model.optimize()
    if sol.status != "optimal" or sol.objective_value is None:
        return {m: float("nan") for m in met_ids}, float("nan")
    bio_opt = sol.objective_value
    floor = bio_opt * biomass_frac
    for met_id in met_ids:
        with model:
            bio.lower_bound = max(0.0, floor)
            met = model.metabolites.get_by_id(met_id)
            dm = cobra.Reaction(f"DM_{met_id}")
            dm.lower_bound, dm.upper_bound = 0.0, 1000.0
            model.add_reactions([dm])
            dm.add_metabolites({met: -1.0})
            model.objective = dm
            s = model.optimize()
            out[met_id] = (s.objective_value
                           if s.status == "optimal" and s.objective_value is not None
                           else float("nan"))
    return out, bio_opt
