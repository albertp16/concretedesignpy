# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
RC Column Concrete Jacketing — design boundary
==============================================

Orchestrates :mod:`concretedesignpy.calculators.column_jacket` (the engine)
into a single design report, and is the ONLY place units are converted.

    engine  : strictly mm, MPa, N, N*mm
    this API: mm, MPa, kN, kN*m

Nothing above this module should ever see newtons, and nothing below it should
ever see kilonewtons.  ``N_PER_KN`` and ``NMM_PER_KNM`` exist so the two
conversions are greppable.  The failure is silent otherwise: a report saying
``9737`` is right and ``9737000`` is also plausible-looking, and nothing
downstream would catch it.

Three properties of this module are load-bearing, not stylistic:

1.  **Validation is a safety layer.**  The engine will integrate a 5 mm jacket
    with 900 MPa reinforcement and return a confident number.
    :func:`validate_jacket_request` is where that gets refused.  Widening a
    bound is a decision, not a fix.
2.  **Advisories are part of the answer.**  Every result carries an
    ``advisories`` list of engineering traps the arithmetic cannot rule out.
    ``adequate`` is a statement about the COMPUTED CHECKS ONLY — it is
    deliberately not called ``safe`` or ``passes``.  A caller rendering only
    the DCRs is using this module wrong.
3.  **A refused check is recorded, not faked.**  The engine raises
    ``NotImplementedError`` for checks it has no model for (a partial jacket
    has no confinement, shear, interface, stiffness or detailing model).  Those
    land in ``unavailable`` and hold ``adequate`` false, because a section
    whose shear was never checked has not been shown adequate.

Reference: TN-RET-001 (APEC technical note, RC Column Concrete Jacketing),
ACI 318-19, ACI 562-16, ASCE 41-17, Mander et al. (1988) via fib Bulletin 14,
Lampropoulos, Tsioulou & Dritsos (2012).
"""

from __future__ import annotations

import hashlib
import math
from pathlib import Path

from concretedesignpy.calculators.column_jacket import (
    Concrete,
    JacketedColumn,
    Steel,
    TieSet,
    capacity_at_P,
    rect_perimeter_bars,
    utilisation,
)

__all__ = [
    "column_jacket_design",
    "build_jacketed_column",
    "validate_jacket_request",
    "provenance",
    "ENGINE_SHA256",
    "N_PER_KN",
    "NMM_PER_KNM",
]

# ──────────────────────────────────────────────
# Unit boundary — the only two conversions in the module
# ──────────────────────────────────────────────
N_PER_KN = 1e3
NMM_PER_KNM = 1e6


# ──────────────────────────────────────────────
# Provenance
# ──────────────────────────────────────────────
# The engine is a VENDORED copy.  That breaks single-source-of-truth by
# construction, so the hash is stamped into every result and drift becomes
# visible to every caller.  Do not edit column_jacket.py to fix a bug: fix it
# upstream, re-vendor the whole file, and bump ENGINE_VERSION here.
ENGINE_VERSION = "TN-RET-001-Rev01 + TN-RET-002 partial jackets (phases 1-3, 8)"
ENGINE_SOURCE_REF = (
    "Vendored verbatim from the APEC RC Column Concrete Jacketing service "
    "(app/engine/column_jacket.py), which itself vendors "
    "APEC knowledge/05_WORKING_FILES/Drafts/2026-08-07_RC_Column_Jacketing_Calc/"
    "column_jacket.py @ 2026-08-08. The file is byte-identical to that copy, so "
    "engine_sha256 compares directly against it."
)

_ENGINE_PATH = Path(__file__).with_name("column_jacket.py")

try:
    ENGINE_SHA256 = hashlib.sha256(_ENGINE_PATH.read_bytes()).hexdigest()
except OSError:  # pragma: no cover - only when the source is not on disk
    ENGINE_SHA256 = ""


def provenance():
    """Where every number in a result came from."""
    return {
        "engine_version": ENGINE_VERSION,
        "engine_sha256": ENGINE_SHA256,
        "engine_source": ENGINE_SOURCE_REF,
        "units": "mm, MPa, N internally; kN and kN.m in results",
        "code_basis": (
            "ACI 318-19 (SI, App. C, cross-checked vs ACI 318-25M); "
            "ACI 562-16 5.4/7.3/7.4; ASCE 41-17 material and classification "
            "provisions; Mander et al. (1988) via fib Bulletin 14"
        ),
        "disclaimer": (
            "Implements code equations only. Does not exercise engineering "
            "judgement and does not replace the Engineer of Record. Read "
            "TN-RET-001 Limitations and Common Errors before use."
        ),
    }


# ──────────────────────────────────────────────
# Input specification — bounds are the safety layer
# ──────────────────────────────────────────────
# (default, gt, ge, lt, le, integer)  — default None means the field is required.
_REQ = object()

_EXISTING_SPEC = {
    "width":                dict(default=_REQ, gt=100, le=3000),
    "depth":                dict(default=_REQ, gt=100, le=3000),
    "fc":                   dict(default=_REQ, gt=5, le=90),
    "fy":                   dict(default=_REQ, gt=200, le=600),
    "bar_dia":              dict(default=_REQ, gt=5, le=60),
    "bars_per_face_width":  dict(default=_REQ, ge=2, le=12, integer=True),
    "bars_per_face_depth":  dict(default=_REQ, ge=2, le=12, integer=True),
    "cover_to_bar_centre":  dict(default=_REQ, gt=15, le=200),
}

_JACKET_SPEC = {
    "thickness":                dict(default=_REQ, ge=50, le=500),
    "fc":                       dict(default=_REQ, gt=5, le=90),
    "fy":                       dict(default=_REQ, gt=200, le=600),
    "bar_dia":                  dict(default=_REQ, gt=5, le=60),
    "bars_per_face_width":      dict(default=_REQ, ge=2, le=12, integer=True),
    "bars_per_face_depth":      dict(default=_REQ, ge=2, le=12, integer=True),
    "clear_cover_to_tie":       dict(default=40.0, ge=15, le=120),
    "tie_dia":                  dict(default=_REQ, ge=6, le=32),
    "tie_spacing":              dict(default=_REQ, gt=25, le=600),
    "tie_fy":                   dict(default=_REQ, gt=200, le=600),
    "tie_legs_each_way":        dict(default=2, ge=2, le=8, integer=True),
    "bars_restrained_per_face": dict(default=2, ge=2, le=12, integer=True),
    "sides":                    dict(default="four",
                                     choices=("four", "three", "two", "one",
                                              "two_adjacent")),
    "spiral":                   dict(default=False, boolean=True),
}

_DEMAND_SPEC = {
    "Pu":           dict(default=_REQ, gt=-1e6, lt=1e6),
    "Mu":           dict(default=_REQ, ge=0, lt=1e6),
    "Vu":           dict(default=_REQ, ge=0, lt=1e6),
    "clear_height": dict(default=_REQ, gt=300, le=20000),
}

_CONSTRUCTION_SPEC = {
    "P0_at_casting":     dict(default=0.0, ge=0),
    "P_service":         dict(default=0.0, ge=0),
    "creep_coefficient": dict(default=2.0, ge=0, le=5),
    "continuity":        dict(default="discontinuous",
                              choices=("continuous", "discontinuous")),
    "interface_surface": dict(default="roughened_6mm",
                              choices=("roughened_6mm", "not_roughened",
                                       "monolithic", "steel")),
    "dowel_dia":         dict(default=16.0, ge=8, le=40),
    "dowels_per_row":    dict(default=2, ge=1, le=12, integer=True),
    "dowel_row_spacing": dict(default=250.0, gt=25, le=1200),
    "dowel_fy":          dict(default=415.0, gt=200, le=600),
}


def _validate_group(name, data, spec):
    """Apply defaults and bounds to one input group, or raise ValueError."""
    data = {} if data is None else dict(data)

    unknown = sorted(set(data) - set(spec))
    if unknown:
        raise ValueError(
            f"{name}: unknown field(s) {unknown}; expected one of "
            f"{sorted(spec)}"
        )

    out = {}
    for field, rule in spec.items():
        if field not in data:
            if rule["default"] is _REQ:
                raise ValueError(f"{name}.{field} is required")
            out[field] = rule["default"]
            continue

        value = data[field]

        if rule.get("choices") is not None:
            if value not in rule["choices"]:
                raise ValueError(
                    f"{name}.{field} must be one of {list(rule['choices'])}, "
                    f"got {value!r}"
                )
            out[field] = value
            continue

        if rule.get("boolean"):
            if not isinstance(value, bool):
                raise ValueError(f"{name}.{field} must be true or false, "
                                 f"got {value!r}")
            out[field] = value
            continue

        if isinstance(value, bool) or not isinstance(value, (int, float)):
            raise ValueError(f"{name}.{field} must be a number, got {value!r}")
        value = float(value)
        # A NaN or inf demand silently flips comparisons -- including the two
        # mandatory ACI 562 interface requirements, which are `>` tests.
        if math.isnan(value) or math.isinf(value):
            raise ValueError(f"{name}.{field} must be finite, got {value!r}")
        if rule.get("integer"):
            if value != int(value):
                raise ValueError(f"{name}.{field} must be a whole number, "
                                 f"got {value!r}")
            value = int(value)

        for key, op, text in (("gt", float.__le__, "greater than"),
                              ("ge", float.__lt__, "at least"),
                              ("lt", float.__ge__, "less than"),
                              ("le", float.__gt__, "at most")):
            bound = rule.get(key)
            if bound is not None and op(float(value), float(bound)):
                raise ValueError(
                    f"{name}.{field} must be {text} {bound:g}, got {value:g}"
                )
        out[field] = value

    return out


def validate_jacket_request(existing, jacket, demand, construction=None):
    """Validate and normalise the four input groups.

    Returns a dict of four dicts with defaults applied.  Raises ``ValueError``
    with an explanatory message for anything the model must not be pushed to.

    Two of the cross-field checks do real work and must not be relaxed
    casually:

    * the existing bar cage must retain a real lever arm.  Checking only
      ``2*cover < dimension`` is not enough — a 400 mm column with 199 mm cover
      passes that while placing every bar within 2 mm of mid-depth, and the
      engine would integrate it happily and return a confident, meaningless
      flexural capacity.
    * a jacket must be thick enough to physically house its own bars inside
      its own ties.
    """
    e = _validate_group("existing", existing, _EXISTING_SPEC)
    j = _validate_group("jacket", jacket, _JACKET_SPEC)
    d = _validate_group("demand", demand, _DEMAND_SPEC)
    c = _validate_group("construction", construction, _CONSTRUCTION_SPEC)

    clear = min(e["width"], e["depth"]) - 2 * e["cover_to_bar_centre"]
    if clear < 5 * e["bar_dia"]:
        raise ValueError(
            f"existing.cover_to_bar_centre {e['cover_to_bar_centre']:.0f} mm "
            f"leaves only {clear:.0f} mm between the outer bar rows of a "
            f"{min(e['width'], e['depth']):.0f} mm section; at least "
            f"{5 * e['bar_dia']:.0f} mm (5 bar diameters) is required for the "
            "section to have a meaningful lever arm"
        )

    # A jacket on two ADJACENT faces is unsymmetric about both principal axes,
    # so its capacity is a biaxial SURFACE rather than a curve.  Refused by
    # name rather than approximated.
    if j["sides"] == "two_adjacent":
        raise ValueError(
            "A jacket on two ADJACENT faces is unsymmetric about both "
            "principal axes, so its capacity is a biaxial SURFACE, not a "
            "uniaxial interaction curve. How that surface is constructed is "
            "unresolved -- TN-RET-002 TODO C2, BLOCKING. This module refuses "
            "rather than returning a curve that silently ignores M_z. One, "
            "two opposite, three and four faces are computed."
        )

    if c["P0_at_casting"] > d["Pu"]:
        raise ValueError(
            "P0_at_casting exceeds Pu; the load present at casting cannot be "
            "larger than the factored demand"
        )

    if j["thickness"] < j["bar_dia"] + 2 * j["tie_dia"] + 20:
        raise ValueError(
            f"jacket thickness {j['thickness']:.0f} mm cannot house a "
            f"{j['bar_dia']:.0f} mm bar inside {j['tie_dia']:.0f} mm ties with "
            "cover; increase thickness or reduce bar sizes"
        )

    return {"existing": e, "jacket": j, "demand": d, "construction": c}


# ──────────────────────────────────────────────
# Section assembly
# ──────────────────────────────────────────────
def build_jacketed_column(existing, jacket, demand=None, construction=None):
    """Assemble the engine's :class:`JacketedColumn` from validated inputs.

    Accepts the same input groups as :func:`column_jacket_design`.  ``demand``
    and ``construction`` are not used for the geometry; they are accepted so a
    caller can hand the same request to both functions.
    """
    req = validate_jacket_request(
        existing, jacket,
        demand if demand is not None else dict(Pu=0.0, Mu=0.0, Vu=0.0,
                                               clear_height=3000.0),
        construction,
    )
    return _build_column(req)


def _build_column(req):
    e, j = req["existing"], req["jacket"]

    conc_e = Concrete(fc=e["fc"], label="existing")
    conc_j = Concrete(fc=j["fc"], label="jacket")
    steel_e = Steel(fy=e["fy"])
    steel_j = Steel(fy=j["fy"])

    # Which faces are enlarged.  The per-face thicknesses default to
    # `thickness`, so a four-sided jacket builds exactly the section it always
    # did (TN-RET-002 / CODE_IMPLEMENTATION_PARTIAL_JACKET s2.1).
    t = j["thickness"]
    t_top, t_bot, t_left, t_right = {
        "four":  (t, t, t, t),
        "three": (t, t, t, 0.0),      # J3-U: top, bottom and one long face
        "two":   (t, t, 0.0, 0.0),    # J2-o: opposite faces
        "one":   (t, 0.0, 0.0, 0.0),  # J1-c: one face
    }[j["sides"]]

    B = e["width"] + t_left + t_right
    H = e["depth"] + t_top + t_bot
    z_core_left = -B / 2 + t_left
    zc_e = z_core_left + e["width"] / 2     # core centroid about the composite
    faces = (t_top > 0, t_bot > 0, t_left > 0, t_right > 0)

    bars_e = rect_perimeter_bars(
        e["width"], e["depth"], e["cover_to_bar_centre"],
        e["bars_per_face_width"], e["bars_per_face_depth"], e["bar_dia"],
        steel_e, "existing", t_top, zc_e,
    )
    # Cover to the jacket BAR CENTRE = clear cover + tie diameter + half bar.
    # Deriving it rather than accepting it prevents the geometrically
    # impossible layouts that arise when a caller passes clear cover as if it
    # were to centre.
    cover_to_jacket_bar = j["clear_cover_to_tie"] + j["tie_dia"] + j["bar_dia"] / 2.0
    bars_j = rect_perimeter_bars(
        B, H, cover_to_jacket_bar,
        j["bars_per_face_width"], j["bars_per_face_depth"], j["bar_dia"],
        steel_j, "jacket", 0.0, 0.0, faces,
    )

    leg_area = math.pi * j["tie_dia"] ** 2 / 4.0 * j["tie_legs_each_way"]
    ties = TieSet(Av_x=leg_area, Av_z=leg_area, s=j["tie_spacing"],
                  fyt=j["tie_fy"], db=j["tie_dia"])

    return JacketedColumn(
        b_e=e["width"], h_e=e["depth"], t=t,
        t_top=t_top, t_bot=t_bot, t_left=t_left, t_right=t_right,
        conc_e=conc_e, conc_j=conc_j,
        bars_e=bars_e, bars_j=bars_j, ties_j=ties,
        cover_j=j["clear_cover_to_tie"], spiral=j["spiral"],
        n_sup_b=j["bars_restrained_per_face"],
        n_sup_h=j["bars_restrained_per_face"],
        db_long=j["bar_dia"],
    )


# ──────────────────────────────────────────────
# Private helpers
# ──────────────────────────────────────────────
def _finite(x):
    """None for anything that cannot be reported as a number."""
    if x is None:
        return None
    x = float(x)
    return None if (math.isnan(x) or math.isinf(x)) else x


def _try(unavailable, what, fn):
    """Run an engine check, or record why it refused.

    The engine raises ``NotImplementedError`` for a check whose partial-jacket
    model does not exist, rather than answering with four-sided behaviour.
    Swallowing that into a number is precisely the failure this boundary exists
    to prevent, so the absence is carried into the result and ``adequate``
    stays false while any remain.
    """
    try:
        return fn()
    except NotImplementedError as exc:
        unavailable.append("{}: {}".format(what, exc))
        return None


def _advisory(severity, code, message, clause=None):
    return {"severity": severity, "code": code, "message": message,
            "clause": clause}


# ──────────────────────────────────────────────
# Main public function
# ──────────────────────────────────────────────
def column_jacket_design(existing, jacket, demand, construction=None,
                         n_points=90, n_layers=600):
    """Full design report for an RC column strengthened by a concrete jacket.

    Parameters
    ----------
    existing : dict
        The column AS MEASURED, not as drawn — feed it from a condition
        assessment.  ``width``, ``depth`` (mm), ``fc``, ``fy`` (MPa),
        ``bar_dia`` (mm), ``bars_per_face_width``, ``bars_per_face_depth``,
        ``cover_to_bar_centre`` (mm to bar centre).
    jacket : dict
        ``thickness`` (mm), ``fc``, ``fy`` (MPa), ``bar_dia`` (mm),
        ``bars_per_face_width``, ``bars_per_face_depth``,
        ``clear_cover_to_tie`` (mm, default 40), ``tie_dia`` (mm),
        ``tie_spacing`` (mm), ``tie_fy`` (MPa), ``tie_legs_each_way``
        (default 2), ``bars_restrained_per_face`` (default 2),
        ``sides`` (``"four"`` default / ``"three"`` / ``"two"`` opposite /
        ``"one"``), ``spiral`` (default False).

        ``bars_restrained_per_face`` drives the Mander k_e and is the input
        most often set wrong.  A plain perimeter hoop laterally restrains only
        the four corner bars — that is 2 per face.  Restraining every bar
        requires crossties (ACI 318-19 25.7.2.3); assuming full restraint when
        the detail is a plain hoop overstates f'cc by roughly 12%.
    demand : dict
        Factored demand on the JACKETED column: ``Pu`` (kN, compression
        positive), ``Mu`` (kN-m), ``Vu`` (kN), ``clear_height`` (mm).

        The trap this cannot check for you: after jacketing the column is
        several times stiffer and attracts more force.  If these values came
        from a model of the *un*-jacketed structure, every result below is
        optimistic.
    construction : dict, optional
        ``P0_at_casting`` (kN, sustained load present when the jacket is cast;
        0 = shored), ``P_service`` (kN), ``creep_coefficient`` (default 2.0),
        ``continuity`` (``"discontinuous"`` default / ``"continuous"``),
        ``interface_surface``, ``dowel_dia``, ``dowels_per_row``,
        ``dowel_row_spacing``, ``dowel_fy``.

        ``continuity`` is the single most consequential field.  A jacket
        stopping at the slab soffit must transfer its entire axial share across
        the interface; one cored through the slab and bearing on the foundation
        transfers it by direct bearing and the axial term vanishes.
    n_points, n_layers : int
        Interaction-diagram resolution.

    Returns
    -------
    dict
        ``provenance``, ``request_echo``, ``section``, ``geometry``, ``axial``,
        ``interaction``, ``monolithic``, ``shear``, ``confinement``,
        ``preload``, ``interface``, ``stiffness``, ``detailing``,
        ``advisories``, ``governing_checks``, ``unavailable``, ``adequate``.

        Forces are kN and moments kN-m throughout.  ``adequate`` covers the
        COMPUTED CHECKS ONLY — read ``advisories`` before relying on it.

    Raises
    ------
    ValueError
        For any input the model must not be pushed to (see
        :func:`validate_jacket_request`).
    """
    req = validate_jacket_request(existing, jacket, demand, construction)
    e, j, d, c = (req["existing"], req["jacket"], req["demand"],
                  req["construction"])
    col = _build_column(req)

    # Checks the engine refuses for this section are recorded, not faked.
    unavailable = []

    Pu = d["Pu"] * N_PER_KN
    Mu = d["Mu"] * NMM_PER_KNM
    Vu = d["Vu"] * N_PER_KN
    L = d["clear_height"]

    ia_e = col.interaction_code(existing_only=True, n_points=n_points,
                                n_layers=n_layers)
    ia_j = col.interaction_code(n_points=n_points, n_layers=n_layers)

    # ---------------- section ---------------------------------------------
    section = {
        "B": col.B, "H": col.H,
        "Ag_existing": col.Ag_e, "Ag_jacketed": col.Ag_tot,
        # the annulus alone, so a client never has to derive Ag_tot - Ag_e
        "Ag_jacket_annulus": col.Ag_j,
        "Ast_existing": col.Ast_e, "Ast_jacket": col.Ast_j,
        "rho_total_percent": col.Ast_tot / col.Ag_tot * 100.0,
        "n_bars_existing": len(col.bars_e), "n_bars_jacket": len(col.bars_j),
        "beta1_existing": col.conc_e.beta1, "beta1_jacket": col.conc_j.beta1,
    }

    # ---------------- axial -----------------------------------------------
    Po_e = col.Po_existing_only()
    Po_j = col.Po()
    phiPn_e = 0.65 * 0.80 * Po_e
    axial = {
        "Po_existing": Po_e / N_PER_KN,
        "Po_jacketed": Po_j / N_PER_KN,
        "phiPn_max_existing": phiPn_e / N_PER_KN,
        "phiPn_max_jacketed": col.phiPn_max() / N_PER_KN,
        "strength_gain": Po_j / Po_e,
        "Pnt_max": col.Pnt_max() / N_PER_KN,
        "clause": "ACI 318-19 22.4.2, Table 22.4.2.1, Table 21.2.2",
    }

    # ---------------- interaction ------------------------------------------
    def pts(P, M):
        order = sorted(range(len(P)), key=lambda i: -P[i])
        return [{"P": float(P[i]) / N_PER_KN, "M": float(M[i]) / NMM_PER_KNM}
                for i in order]

    Me = capacity_at_P(ia_e["phiP"], ia_e["phiM"], Pu)
    Mj = capacity_at_P(ia_j["phiP"], ia_j["phiM"], Pu)
    u_e = utilisation(ia_e["phiP"], ia_e["phiM"], Pu, Mu)
    u_j = utilisation(ia_j["phiP"], ia_j["phiM"], Pu, Mu)
    interaction = {
        "design_existing": pts(ia_e["phiP"], ia_e["phiM"]),
        "design_jacketed": pts(ia_j["phiP"], ia_j["phiM"]),
        "nominal_jacketed": pts(ia_j["P"], ia_j["M"]),
        # None, never a clamped number: the section cannot carry that axial
        # load at any eccentricity.
        "phiMn_at_Pu_existing": _finite(Me / NMM_PER_KNM if not math.isnan(Me)
                                        else float("nan")),
        "phiMn_at_Pu_jacketed": _finite(Mj / NMM_PER_KNM if not math.isnan(Mj)
                                        else float("nan")),
        "utilisation_existing": _finite(u_e),
        "utilisation_jacketed": _finite(u_j),
        "stress_block": "single block, beta1 from the extreme compression fibre",
        "clause": "ACI 318-19 22.2, Table 21.2.2",
    }

    # ---------------- shear -------------------------------------------------
    sh = _try(unavailable, "shear", lambda: col.shear_capacity(Nu=Pu))
    Mn = capacity_at_P(ia_j["P"], ia_j["M"], Pu)
    Vp = (col.shear_demand_from_hinging(Mn, Mn, L) if not math.isnan(Mn)
          else float("nan"))
    flexure_ok = (sh is not None and not math.isnan(Vp)
                  and sh["phiVn"] >= Vp)
    shear = None if sh is None else {
        "Vc": sh["Vc"] / N_PER_KN, "Vs": sh["Vs"] / N_PER_KN,
        "phiVn": sh["phiVn"] / N_PER_KN,
        "d_used": sh["d"], "fc_used": sh["fc_used"],
        "utilisation": (Vu / sh["phiVn"] if sh["phiVn"] > 0 else float("inf")),
        "Vp_at_hinging": (Vp / N_PER_KN) if not math.isnan(Vp) else 0.0,
        "phiVn_over_Vp": (sh["phiVn"] / Vp
                          if (not math.isnan(Vp) and Vp > 0) else 0.0),
        "failure_mode": "flexure-controlled" if flexure_ok else "shear-critical",
        "clause": "ACI 318-19 22.5.1.1, Table 22.5.5.1(a) SI, 22.5.8.5.3",
    }

    # ---------------- confinement -------------------------------------------
    cf = _try(unavailable, "confinement", lambda: col.confined("existing"))
    confinement = None if cf is None else {
        "rho_s_percent": col.rho_s_ties * 100.0,
        "ke_computed": col.mander_ke(),
        "fl": cf.fl, "fcc_core": cf.fcc,
        "strength_gain": cf.strength_gain,
        "eps_cu_core": cf.eps_cu,
        "eps_cu_unconfined": 0.0038,
        "clause": ("Mander et al. (1988) via fib Bulletin 14 "
                   "Eq. (6-5)/(6-35)/(6-36)"),
    }

    # ---------------- preload ------------------------------------------------
    P0 = c["P0_at_casting"] * N_PER_KN
    Pser = (c["P_service"] or d["Pu"] * 0.7) * N_PER_KN
    ss = _try(unavailable, "preload",
              lambda: col.service_stress_split(P0, Pser,
                                               c["creep_coefficient"]))
    preload = None if ss is None else {
        "eps_0": col.preload_strain(P0, c["creep_coefficient"]),
        "sigma_core": ss["sigma_core"],
        "sigma_core_monolithic": ss["sigma_core_monolithic"],
        "core_overstress": ss["core_overstress"],
        "jacket_share_percent": ss["jacket_share"] * 100.0,
        "jacket_share_monolithic_percent": ss["jacket_share_monolithic"] * 100.0,
        "note": ("Only load applied after casting is shared; P0 stays in the "
                 "core. Shoring before casting removes this effect entirely."),
    }

    # ---------------- interface ----------------------------------------------
    p_int = 2.0 * (col.b_e + col.h_e)
    row_area = math.pi * c["dowel_dia"] ** 2 / 4.0 * c["dowels_per_row"]
    rho_v = row_area / (p_int * c["dowel_row_spacing"])
    dP_j = max(Pu - phiPn_e, 0.0)
    C_j = 0.55 * (Mn if not math.isnan(Mn) else 0.0) / (0.80 * col.H) * 2
    it = _try(unavailable, "interface", lambda: col.interface_check(
        C_jacket=C_j, dP_jacket=dP_j, Vu=Vu,
        L_shear_span=L / 2.0, L_transfer=L,
        rho_v=rho_v, fy_dowel=c["dowel_fy"],
        continuity=c["continuity"], surface=c["interface_surface"],
    ))
    s_req = (row_area / (it["rho_v_req_dowel"] * p_int)
             if it is not None and it["rho_v_req_dowel"] > 0 else None)
    interface = None if it is None else {
        "continuity": c["continuity"],
        "v_u": it["v_u"], "v_flexural": it["v_flex"], "v_axial": it["v_axial"],
        "mu": it["mu"],
        "rho_v_provided_percent": rho_v * 100.0,
        # BOTH capacity routes are reported and NEITHER is picked. Crediting
        # the 1.8 MPa cohesion term is an Engineer-of-Record decision.
        "phi_vn_dowels_only": it["phi_v_n_dowel"],
        "DCR_dowels_only": it["DCR_dowel"],
        "phi_vn_with_bond": it["phi_v_n_bond"],
        "DCR_with_bond": it["DCR_bond"],
        "dowel_spacing_required_dowels_only": _finite(s_req),
        "relies_on_bond": (not it["ok_dowel_only"]) and it["ok_with_bond"],
        "aci562_bond_testing_required": it["aci562_bond_test_required"],
        "aci562_interface_reinforcement_required": it["aci562_reinf_required"],
        "clause": "ACI 318-19 16.4.4/22.9; ACI 562-16 7.4, Table 7.4.1.2",
    }

    # ---------------- geometry, for drawing ----------------------------------
    # Bar coordinates come from the engine's own bar objects, not regenerated
    # here.  A separately-derived layout could disagree with the section the
    # capacities were integrated over, and the drawing would then quietly
    # misrepresent the analysis.  The tie cage only exists where the jacket
    # wraps the column.
    _cage = _try(unavailable, "tie cage geometry", lambda: col.core_dims)
    b_c, h_c = _cage if _cage is not None else (None, None)
    y0_core, y1_core = col.y_core
    z0_core, z1_core = col.z_core
    geometry = {
        "B": col.B, "H": col.H,
        "existing_b": col.b_e, "existing_h": col.h_e,
        "jacket_thickness": col.t,
        # Where the existing core actually sits inside the composite, taken
        # from the engine rather than assumed centred. For a partial jacket it
        # is NOT centred, and a drawing that centred it would misrepresent the
        # section the capacities were integrated over.
        "existing_y0": y0_core, "existing_y1": y1_core,
        "existing_z0": z0_core, "existing_z1": z1_core,
        "tie_cage_b": b_c, "tie_cage_h": h_c,
        "interface_perimeter": p_int,
        "bars": [
            {"y": bb.y, "z": bb.z, "zone": bb.zone,
             "dia": e["bar_dia"] if bb.zone == "existing" else j["bar_dia"]}
            for bb in col.bars
        ],
    }

    # ---------------- stiffness + detailing ----------------------------------
    st = _try(unavailable, "stiffness", lambda: col.stiffness_summary())
    stiffness = None if st is None else {
        "Ag_ratio": st["Ag_ratio"], "EI_ratio": st["EI_ratio"],
        "local_period_factor": st["period_factor"],
        "note": ("The jacketed column attracts more seismic force. The demand "
                 "model must be re-run with this stiffness or the checks above "
                 "prove nothing."),
    }
    _det = _try(unavailable, "detailing",
                lambda: col.detailing_checks(db_long_min=j["bar_dia"],
                                             L_clear=L))
    detailing = [] if _det is None else [
        {"item": x["item"], "value": x["value"], "limit": x["limit"],
         "clause": x["clause"], "passes": bool(x["ok"])}
        for x in _det
    ]

    # ---------------- monolithic coefficients --------------------------------
    # REPORTED, never applied.  A jacketed column does not reach the strength
    # of a monolithic section, and K_F is that ratio -- 0.743 on the worked
    # example, which takes the utilisation from 0.97 to 1.31.  Whether K_F
    # multiplies phi*M_n or M_n is TN-RET-001 TODO C8/C9 and BLOCKING, so both
    # results are carried side by side and neither is chosen.
    monolithic = _monolithic_block(col, Pu, interaction, stiffness)

    # ---------------- governing checks + advisories --------------------------
    governing = []
    if u_j is None or math.isnan(u_j) or math.isinf(u_j) or u_j > 1.0:
        governing.append("P-M interaction")
    if shear is not None and shear["utilisation"] > 1.0:
        governing.append("shear strength")
    if shear is not None and shear["failure_mode"] == "shear-critical":
        governing.append("shear-critical failure mode")
    if it is not None and not it["ok_with_bond"]:
        governing.append("interface shear transfer")
    governing += ["detailing: {}".format(x["item"]) for x in detailing
                  if not x["passes"]]

    advisories = _advisories(req, it, st, confinement, preload, shear,
                             interface, monolithic)

    qaqc = _qaqc_block(req, section, axial, interaction, shear, confinement,
                       interface, d)

    return {
        "provenance": provenance(),
        "request_echo": req,
        "section": section,
        "geometry": geometry,
        "axial": axial,
        "interaction": interaction,
        "monolithic": monolithic,
        "shear": shear,
        "confinement": confinement,
        "preload": preload,
        "interface": interface,
        "stiffness": stiffness,
        "detailing": detailing,
        "qaqc": qaqc,
        "advisories": advisories,
        "governing_checks": governing,
        "unavailable": unavailable,
        # All COMPUTED checks satisfied. NOT the same as "the design is safe"
        # -- see advisories, and TN-RET-001 Limitations.
        "adequate": (not governing) and (not unavailable),
    }


def _qaqc_close(a, b, rel=1e-6, abs_=1e-6):
    if a is None or b is None:
        return a is None and b is None
    # bool(), because a numpy operand (e.g. d_used, traced back to the
    # engine's np.linspace bar layout) turns the comparison into numpy.bool_,
    # which the JSON encoder refuses.
    return bool(abs(a - b) <= max(abs_, rel * max(abs(a), abs(b))))


def _qaqc_block(req, section, axial, interaction, shear, confinement,
                interface, demand):
    """Independent recomputation of the report's own numbers.

    Every check re-derives a reported value from the RAW REQUEST INPUTS (or,
    where noted, from other reported values) along a separate arithmetic path,
    then compares. A pass means the numbers in this report agree with each
    other; it is NOT a check of the design. These are computed on the server
    so a client can render QAQC without holding a single equation.

    The recomputations are written out longhand on purpose -- reusing the
    engine's own helpers here would verify nothing.
    """
    e, j = req["existing"], req["jacket"]
    checks = []

    def add(name, method, expected, computed, rel=1e-6):
        checks.append({
            "name": name, "method": method,
            "expected": _finite(expected), "computed": _finite(computed),
            "pass": _qaqc_close(expected, computed, rel=rel),
        })

    # -- geometry, from the raw inputs and the sides mapping ---------------
    t = j["thickness"]
    t_top, t_bot, t_left, t_right = {
        "four":  (t, t, t, t),
        "three": (t, t, t, 0.0),
        "two":   (t, t, 0.0, 0.0),
        "one":   (t, 0.0, 0.0, 0.0),
    }[j["sides"]]
    B = e["width"] + t_left + t_right
    H = e["depth"] + t_top + t_bot
    add("Overall width B",
        "b_e + t_left + t_right from the inputs and the faces mapping",
        B, section["B"])
    add("Overall depth H",
        "h_e + t_top + t_bot from the inputs and the faces mapping",
        H, section["H"])
    add("Gross area of composite",
        "B x H from the recomputed overall dimensions",
        B * H, section["Ag_jacketed"])

    # -- reinforcement, from bar counts and diameters ----------------------
    n_e = 2 * e["bars_per_face_width"] + 2 * e["bars_per_face_depth"] - 4
    Ast_e = n_e * math.pi / 4.0 * e["bar_dia"] ** 2
    add("Existing bar count",
        "perimeter cage: 2*n_width + 2*n_depth - 4 (corners once)",
        n_e, section["n_bars_existing"])
    add("Existing steel area",
        "count x (pi/4) d_b^2 from the recomputed count",
        Ast_e, section["Ast_existing"])

    four_sided = j["sides"] == "four"
    if four_sided:
        n_j = 2 * j["bars_per_face_width"] + 2 * j["bars_per_face_depth"] - 4
        Ast_j = n_j * math.pi / 4.0 * j["bar_dia"] ** 2
        add("Jacket bar count",
            "perimeter cage: 2*n_width + 2*n_depth - 4 (corners once)",
            n_j, section["n_bars_jacket"])
        add("Jacket steel area",
            "count x (pi/4) d_b^2 from the recomputed count",
            Ast_j, section["Ast_jacket"])

        # -- axial, longhand from ACI 318-19 Eq. (22.4.2.2) ----------------
        Ag_e = e["width"] * e["depth"]
        Ag_j = B * H - Ag_e
        fy_e = min(e["fy"], 550.0)
        fy_j = min(j["fy"], 550.0)
        Po = (0.85 * e["fc"] * (Ag_e - Ast_e) + fy_e * Ast_e
              + 0.85 * j["fc"] * (Ag_j - Ast_j) + fy_j * Ast_j)
        add("Nominal axial strength Po",
            "0.85 f'c (Ag - Ast) + fy Ast per concrete region, fy capped at "
            "550 MPa, every area re-derived from the raw inputs",
            Po / N_PER_KN, axial["Po_jacketed"])
        phi = 0.75 if j["spiral"] else 0.65
        alpha = 0.85 if j["spiral"] else 0.80
        add("Design axial ceiling phi*Pn,max",
            "phi x alpha x Po with the recomputed Po "
            "(ACI 318-19 22.4.2.1, Table 21.2.2)",
            phi * alpha * Po / N_PER_KN, axial["phiPn_max_jacketed"])

    add("Total steel ratio",
        "(Ast_e + Ast_j) / (B x H), areas as reported",
        (section["Ast_existing"] + section["Ast_jacket"])
        / section["Ag_jacketed"] * 100.0,
        section["rho_total_percent"])
    add("Axial strength gain",
        "Po_jacketed / Po_existing, both as reported",
        axial["Po_jacketed"] / axial["Po_existing"], axial["strength_gain"])

    # -- interaction: interpolate the returned curve ourselves -------------
    pts = sorted(interaction["design_jacketed"], key=lambda p: p["P"])
    Pu_kN = demand["Pu"]
    phiMn = None
    if pts and pts[0]["P"] <= Pu_kN <= pts[-1]["P"]:
        for lo, hi in zip(pts, pts[1:]):
            if lo["P"] <= Pu_kN <= hi["P"]:
                f = (0.0 if hi["P"] == lo["P"]
                     else (Pu_kN - lo["P"]) / (hi["P"] - lo["P"]))
                phiMn = lo["M"] + f * (hi["M"] - lo["M"])
                break
    add("phi*Mn at Pu",
        "piecewise-linear interpolation of the returned design_jacketed "
        "points at Pu -- the same points the chart draws",
        phiMn, interaction["phiMn_at_Pu_jacketed"], rel=1e-5)
    util = (abs(demand["Mu"]) / phiMn if phiMn else None)
    add("Utilisation at (Pu, Mu)",
        "Mu / (phi*Mn at Pu) with the independently interpolated capacity",
        util, interaction["utilisation_jacketed"], rel=1e-5)

    # -- shear: full longhand recomputation --------------------------------
    if shear is not None:
        fc = shear["fc_used"]
        d_ = shear["d_used"]
        bw = section["B"]
        Ag = section["Ag_jacketed"]
        Nu = max(demand["Pu"] * N_PER_KN, 0.0)
        ax_term = min(Nu / (6.0 * Ag), 0.05 * fc)
        Vc = min((0.17 * math.sqrt(fc) + ax_term) * bw * d_,
                 0.42 * math.sqrt(fc) * bw * d_)
        Av = math.pi / 4.0 * j["tie_dia"] ** 2 * j["tie_legs_each_way"]
        Vs = Av * j["tie_fy"] * d_ / j["tie_spacing"]
        Vn = min(Vc + Vs, Vc + 0.66 * math.sqrt(fc) * bw * d_)
        add("Design shear strength phi*Vn",
            "Table 22.5.5.1(a) SI longhand: Vc with the axial term capped at "
            "0.05 f'c and the 0.42 sqrt(f'c) ceiling, Vs = Av fyt d / s, "
            "Vn capped per 22.5.1.2, phi = 0.75",
            0.75 * Vn / N_PER_KN, shear["phiVn"])

    # -- interface: the demand must be the sum of its parts ----------------
    if interface is not None:
        add("Interface demand v_u",
            "v_flexural + v_axial, both as reported",
            interface["v_flexural"] + interface["v_axial"], interface["v_u"])

    # -- confinement: Mander closed form from the reported f_l -------------
    if confinement is not None:
        fco = e["fc"]
        r = confinement["fl"] / fco
        fcc = max(fco * (2.254 * math.sqrt(1.0 + 7.94 * r) - 2.0 * r - 1.254),
                  fco)
        add("Confined strength f'cc",
            "fib B14 Eq. (6-5) evaluated from the reported f_l and the "
            "existing f'co, with the f'cc >= f'co clamp",
            fcc, confinement["fcc_core"])

    return {
        "checks": checks,
        "n_pass": sum(1 for c in checks if c["pass"]),
        "all_pass": all(c["pass"] for c in checks),
        "note": ("Each row re-derives a reported value along an independent "
                 "arithmetic path and compares. A pass means this report is "
                 "internally consistent; it is not a statement about the "
                 "design."),
    }


def _monolithic_block(col, Pu, interaction, stiffness):
    """Lampropoulos et al. (2012) coefficients, reported and never applied."""
    raw = col.monolithic_coefficients(Pu)
    do_not_apply_to = [
        "capacity design: V_e = (M_pr,top + M_pr,bot) / l_u needs an UPPER "
        "bound on M_pr. K_F <= 1 always, so applying it LOWERS V_e -- "
        "conservative-looking and unconservative in fact.",
        "any member whose pushover peak is not flexural. The source never "
        "separates a shear- or slip-governed peak, so there is no basis for "
        "K_F at all; fix the deficiency first.",
    ]
    clause = ("Lampropoulos, Tsioulou & Dritsos (2012) Tables 2-5; "
              "TN-RET-001 Core Technical Requirements 4")

    if not raw["applicable"]:
        return {"applicable": False, "n_faces": raw["n_faces"],
                "reason": raw["reason"], "do_not_apply_to": do_not_apply_to,
                "clause": clause}

    kf = raw["coefficients"]["K_F"]
    kk = raw["coefficients"]["K_k"]
    phiMn_mono = interaction["phiMn_at_Pu_jacketed"]
    u_mono = interaction["utilisation_jacketed"]
    out = {
        "applicable": True, "n_faces": raw["n_faces"], "reason": raw["reason"],
        "nu": raw["nu"], "r": raw["r"],
        "out_of_range": raw["out_of_range"], "orders_agree": raw["orders_agree"],
        "phiMn_monolithic": phiMn_mono,
        "phiMn_x_KF": (phiMn_mono * kf) if phiMn_mono is not None else None,
        "utilisation_monolithic": u_mono,
        "utilisation_x_KF": (u_mono / kf) if u_mono is not None else None,
        "EI_ratio_x_Kk": ((stiffness["EI_ratio"] * kk)
                          if stiffness is not None else None),
        "ductility_factor": raw["coefficients"]["K_du"] / raw["coefficients"]["K_dy"],
        "large_jacket_caution": raw["r"] > 2.0,
        "do_not_apply_to": do_not_apply_to,
        "clause": clause,
    }
    out.update(raw["coefficients"])
    return out


def _advisories(req, it, st, confinement, preload, shear, interface,
                monolithic):
    """The findings a number alone does not carry.

    Each one here is a real trap that the arithmetic cannot rule out.  They are
    part of the answer — never drop, downgrade or make them opt-in.
    """
    out = []

    if monolithic["applicable"]:
        # The effect is quantified only where there IS a capacity at P_u.
        # Beyond the section's axial range phi*M_n is None, and formatting it
        # would crash on exactly the overloaded column this advisory matters
        # most for.
        if monolithic["phiMn_monolithic"] is None:
            effect = (
                "The effect cannot be quantified here: P_u lies outside the "
                "section's axial range, so there is no phi*M_n at this axial "
                "load to correct."
            )
        else:
            effect = (
                "Applying K_F takes phi*M_n from {:.0f} to {:.0f} kN.m".format(
                    monolithic["phiMn_monolithic"], monolithic["phiMn_x_KF"])
                + (" and the utilisation from {:.2f} to {:.2f}.".format(
                    monolithic["utilisation_monolithic"],
                    monolithic["utilisation_x_KF"])
                   if monolithic["utilisation_monolithic"] is not None else ".")
            )
        out.append(_advisory(
            "critical", "MONOLITHIC_COEFFICIENT_NOT_APPLIED",
            "A jacketed column does not reach the strength of a monolithic "
            "section cast at once. K_F = {:.3f} and K_k = {:.3f} for this "
            "section. {} THIS MODULE HAS NOT APPLIED IT: whether K_F "
            "multiplies phi*M_n or M_n, and whether phi and K_F may both be "
            "applied, is TN-RET-001 TODO C8/C9 and BLOCKING. Both results are "
            "in `monolithic`; choosing between them is an Engineer-of-Record "
            "decision.".format(monolithic["K_F"], monolithic["K_k"], effect),
            "Lampropoulos et al. (2012) Tables 2-5; TN-RET-001 CTR 4"))
        out.append(_advisory(
            "critical", "K_NOT_FOR_CAPACITY_DESIGN",
            "K_F is a strength REDUCTION, so it is safe only where a smaller "
            "strength is the safe answer. Do NOT apply it to a capacity-design "
            "quantity: V_e = (M_pr,top + M_pr,bot)/l_u needs an UPPER bound on "
            "the probable moment, and reducing M_pr lowers V_e -- it looks "
            "conservative and is not. The same applies to any member whose "
            "pushover peak is not flexural, where the source gives no basis "
            "for K_F at all. This module reports Vp_at_hinging WITHOUT K_F "
            "for that reason.",
            "Lampropoulos et al. (2012) sect. 4; ACI 318-19 18.7.6"))
        out.append(_advisory(
            "warning", "DEFORMATION_PENALTY_IS_PROCEDURE_DEPENDENT",
            "K_du/K_dy = {:.3f}: the deformation pair costs {:.0f}% of the "
            "ductility while plastic rotation barely moves. An ASCE 41 LINEAR "
            "procedure accepts on m-factors, which are ductility-like, and "
            "absorbs the whole penalty; a NONLINEAR procedure accepts on "
            "plastic rotation and absorbs almost none. The two will reach "
            "different verdicts on the same column and the linear one is "
            "pessimistic. Report both; if the scheme only closes under NSP, "
            "say that the procedure choice is doing the work.".format(
                monolithic["ductility_factor"],
                (1 - monolithic["ductility_factor"]) * 100),
            "ASCE 41-17 sect. 7.5, Tables 10-7/10-8"))
        if monolithic["large_jacket_caution"]:
            out.append(_advisory(
                "warning", "LARGE_JACKET_COEFFICIENTS_NON_MONOTONIC",
                "A_cj/A_co = {:.2f} is past r = 2, where the tabulated "
                "behaviour reverses: the ductility factor crosses 1.0 near "
                "r = 2.04 and K_k overtakes K_F by r = 3.5. Both reversals are "
                "driven by Table 3's single r = 4.0 anchor "
                "(0.75, 0.75, 1.05, 2.85), whose K_du = 2.85 is an outlier "
                "against about 1.0 everywhere else. State which anchor the "
                "result leans on.".format(monolithic["r"]),
                "Lampropoulos et al. (2012) Table 3"))
        if monolithic["out_of_range"]:
            out.append(_advisory(
                "warning", "K_OUT_OF_RANGE",
                "nu = {:.3f}, r = {:.3f} lie outside the published range "
                "(nu 0.1-0.4, r 0.5-4.0). The coefficients are clamped at the "
                "table anchors, not extrapolated.".format(
                    monolithic["nu"], monolithic["r"]),
                "Lampropoulos et al. (2012)"))
    else:
        out.append(_advisory(
            "critical", "NO_MONOLITHIC_COEFFICIENTS",
            "No published monolithic coefficients cover this jacket. "
            + monolithic["reason"],
            "TN-RET-002 TODO C10 -- BLOCKING"))

    # A partial jacket has no confinement, shear, interface or stiffness model
    # yet.  The advisories that read those are skipped rather than invented;
    # `unavailable` is what records the absence.
    if any(x is None for x in (it, st, confinement, preload, shear, interface)):
        out.append(_advisory(
            "critical", "PARTIAL_JACKET_INCOMPLETE",
            "This section is jacketed on fewer than four faces. Its geometry, "
            "axial strength and P-M interaction are computed. Its confinement, "
            "shear, interface, stiffness and detailing are NOT -- refused "
            "rather than answered with four-sided behaviour; see "
            "`unavailable`. TWO FURTHER OMISSIONS THAT ARE NOT IN THAT LIST. "
            "(1) The axial load arrives through the EXISTING column's centroid "
            "and is resisted about the composite plastic centroid, so it "
            "induces P*e. The engine can compute that eccentricity; this "
            "module does NOT add it to M_u, so the reported DCR is low -- on "
            "TN-RET-002's one-sided case it is 0.82 where the note's own table "
            "gives 1.09, a pass reported where a fail governs. (2) An "
            "unsymmetric section has TWO interaction diagrams and the weaker "
            "governs; only the as-built sense is computed. Both are "
            "TN-RET-002 Common Error 3 and are open findings, not modelled "
            "behaviour. This section has NOT been shown adequate.",
            "TN-RET-002; CODE_IMPLEMENTATION_PARTIAL_JACKET.md s10"))
        return out

    out.append(_advisory(
        "warning", "STIFFNESS_FEEDBACK",
        "Gross area rose x{:.2f} but flexural stiffness rose x{:.2f}. The "
        "jacketed column will attract proportionally more seismic force. "
        "Re-run the demand model with this stiffness; a check against the "
        "pre-retrofit demand proves nothing. Jacketing some columns in a "
        "storey and not others redistributes demand and can create an "
        "irregularity.".format(st["Ag_ratio"], st["EI_ratio"])))

    if req["construction"]["continuity"] == "discontinuous":
        out.append(_advisory(
            "critical", "JACKET_DISCONTINUOUS",
            "The jacket is not continuous, so its entire axial share must "
            "cross the interface: {:.2f} MPa of the {:.2f} MPa total demand. "
            "Coring the slab and bearing the jacket on the foundation removes "
            "this term entirely and is almost always cheaper than the dowels "
            "it replaces.".format(interface["v_axial"], interface["v_u"])))

    if interface["relies_on_bond"]:
        s_req = interface["dowel_spacing_required_dowels_only"]
        out.append(_advisory(
            "critical", "RELIES_ON_INTERFACE_BOND",
            "This design does not close on dowels alone"
            + (" (they would need {:.0f} mm spacing)".format(s_req)
               if s_req else "")
            + ". It relies on the 1.8 MPa cohesion term, which is a BOND "
            "contribution available only from a clean, roughened, well-cured "
            "interface and lost if the interface cracks. Crediting it is an "
            "Engineer-of-Record decision: record it in the calculation, "
            "specify the roughening on the drawings (ACI 318-19 26.5.6.1(c)), "
            "and verify it by pull-off test.",
            "ACI 318-19 Table 16.4.4.2; ACI 562-16 7.4.3"))

    if interface["aci562_interface_reinforcement_required"]:
        out.append(_advisory(
            "warning", "ACI562_INTERFACE_TESTING",
            "v_u = {:.2f} MPa exceeds 0.41 MPa, so interface reinforcement is "
            "required and quantitative bond testing with it: ASTM C1583 "
            "pull-off, minimum 3 tests.".format(interface["v_u"]),
            "ACI 562-16 Table 7.4.1.2, 7.4.3, 7.4.4"))

    out.append(_advisory(
        "warning", "SHRINKAGE_NOT_INCLUDED",
        "Restrained volume change is NOT included in the interface demand. "
        "ACI 562-16 7.4.1 requires it. Jacket shrinkage against the "
        "restraining core puts the interface into tension exactly where the "
        "bond term is being relied on. Assess it from the mix, curing and age "
        "and add it to v_u.",
        "ACI 562-16 7.4.1"))

    if req["construction"]["P0_at_casting"] > 0:
        out.append(_advisory(
            "warning", "UNSHORED_JACKET",
            "The jacket is cast on a loaded column, so it picks up only "
            "{:.0f}% of the axial load rather than the {:.0f}% a monolithic "
            "split implies, and the existing core runs at x{:.2f} the "
            "monolithic stress. Shoring the column before casting removes "
            "this.".format(preload["jacket_share_percent"],
                           preload["jacket_share_monolithic_percent"],
                           preload["core_overstress"])))

    if req["jacket"]["bars_restrained_per_face"] <= 2:
        out.append(_advisory(
            "warning", "NO_CROSSTIES",
            "Only the corner bars are laterally restrained, giving "
            "k_e = {:.2f}. Adding crossties would raise the confined strength "
            "of the existing core materially. If the detail actually has "
            "crossties, set bars_restrained_per_face to match, because "
            "assuming k_e = 0.80 without them overstates f'cc by roughly "
            "12%.".format(confinement["ke_computed"]),
            "ACI 318-19 25.7.2.3"))

    if shear["failure_mode"] == "shear-critical":
        out.append(_advisory(
            "critical", "SHEAR_CRITICAL",
            "phi*Vn = {:.0f} kN is below the shear that comes with flexural "
            "hinging, Vp = {:.0f} kN. Raising Mn without raising Vn converts a "
            "ductile failure into a brittle one. Reduce jacket tie spacing or "
            "add legs.".format(shear["phiVn"], shear["Vp_at_hinging"])))
    elif 1.0 <= shear["phiVn_over_Vp"] < 1.2:
        out.append(_advisory(
            "info", "SHEAR_MARGIN_THIN",
            "phi*Vn/Vp = {:.2f}. Flexure-controlled, but the margin is thin; "
            "it will not survive an increase in Mn.".format(
                shear["phiVn_over_Vp"])))

    out.append(_advisory(
        "info", "FOUNDATION_NOT_CHECKED",
        "The jacketed column delivers more load, and often more moment, to a "
        "footing sized for less. The foundation is outside this module's scope "
        "and must be checked separately."))

    out.append(_advisory(
        "info", "NSCP_CLAUSE_NUMBERING",
        "Clause numbers cited are ACI 318-19. NSCP 2015 Chapter 4 is the "
        "governing Philippine code and adopts ACI 318-14: the values agree, "
        "the numbering does not. Confirm NSCP numbers before issue."))

    return out
