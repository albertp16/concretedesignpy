# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam and joint calculation sheets
=================================

Turns :mod:`beam_moment`, :mod:`beam_shear` and :mod:`joint_shear` into
printable calculation sheets, and is the ONLY place units are converted
for these three.

    calculators : mm, MPa, and whatever each function already takes
    this API    : mm, MPa, kN, kN*m in every reported value

Four properties are load-bearing, not stylistic. They are the same four
that govern :mod:`column_jacket_design`, restated because a reader who
lands here should not have to go and find them.

1.  **QAQC is computed on the server.** Every report carries a ``qaqc``
    block that re-derives its own reported values along a separate
    arithmetic path, written longhand from the printed clause. The client
    renders a table; it holds no equation and compares nothing. Reusing a
    calculator's own helper inside a check would verify nothing, so the
    checks do not.
2.  **Advisories are part of the answer.** ``adequate`` is a statement
    about the COMPUTED CHECKS ONLY. It is deliberately not called
    ``safe`` or ``passes``. A caller rendering only the D/C ratios is
    using this module wrong.
3.  **A check with no model is recorded, not faked.** What a calculator
    does not check lands in ``unavailable`` and holds ``adequate`` false
    when it bears on the verdict. A section whose ``As,min`` was never
    checked has not been shown adequate.
4.  **A demand is optional; a verdict is not invented without one.** Pass
    no demand and the sheet reports capacity, ``adequate`` is ``None``,
    and it says so. It does not quietly score 0/0 as passing.

Governing edition: **NSCP 2015 (= ACI 318M-14)**. See ``CLAUSES.md`` for
the clause register and the page each provision was verified on.
"""

from __future__ import annotations

import math

from concretedesignpy.calculators.beam_moment import calculate_beam_moment
from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength,
    compute_steel_shear_strength,
    shear_design,
)
from concretedesignpy.calculators.joint_shear import (
    joint_shear_check,
    NSCP_2015_TABLE_418_8_4_1,
    PHI_SMF_JOINT,
)

__all__ = [
    "beam_flexure_report",
    "beam_shear_report",
    "joint_shear_report",
    "provenance",
    "N_PER_KN",
    "NMM_PER_KNM",
]

# ──────────────────────────────────────────────
# Unit boundary — the only two conversions in the module
# ──────────────────────────────────────────────
N_PER_KN = 1e3
NMM_PER_KNM = 1e6

REPORT_VERSION = "beam-report-1.0 (note 14 rectification, 2026-08-15)"

CODE_BASIS = (
    "NSCP 2015 (= ACI 318M-14). Clause numbers and printed pages quoted "
    "below are ACI 318-25M, which is the readable, page-verifiable "
    "statement of the same provisions; where the two round differently the "
    "SI value printed in ACI 318M is used."
)


def provenance():
    """Where every number in a report came from."""
    return {
        "report_version": REPORT_VERSION,
        "units": "mm, MPa internally; kN and kN.m in every reported value",
        "code_basis": CODE_BASIS,
        "clause_register": "CLAUSES.md",
        "disclaimer": (
            "Implements code equations only. Does not exercise engineering "
            "judgement and does not replace the Engineer of Record. "
            "'adequate' covers the computed checks only -- read the "
            "advisories and the unavailable list."
        ),
    }


# ──────────────────────────────────────────────
# Small shared helpers
# ──────────────────────────────────────────────

def _f(x):
    """None for a non-finite value rather than a NaN the JSON encoder emits
    as the literal NaN, which is not valid JSON and which every downstream
    parser then disagrees about."""
    if x is None:
        return None
    x = float(x)
    return x if math.isfinite(x) else None


def _r(x, n=2):
    v = _f(x)
    return None if v is None else round(v, n)


def _step(ref, desc, eq, result, status=None):
    """One REFERENCES | CALCULATIONS | RESULT row of the sheet."""
    return {"ref": ref, "desc": desc, "eq": eq, "result": result,
            "status": status}


def _dc(demand, capacity):
    """Demand/capacity, or None when either side is missing.

    Never returns 0.0 for a zero capacity. That substitution is the exact
    defect this audit family keeps finding: a zero capacity scored as
    D/C = 0 reads as passing.
    """
    if demand is None or capacity is None:
        return None
    if capacity <= 0:
        return None
    return demand / capacity


def _verdict(dc):
    if dc is None:
        return None
    return dc <= 1.0


class _Checks:
    """Accumulator for the server-side QAQC block."""

    def __init__(self):
        self.items = []

    def add(self, name, method, expected, computed, rel=1e-6, abs_=1e-9):
        e, c = _f(expected), _f(computed)
        if e is None or c is None:
            ok = e is None and c is None
        else:
            ok = abs(e - c) <= max(abs_, rel * max(abs(e), abs(c)))
        self.items.append({"name": name, "method": method,
                           "expected": e, "computed": c, "pass": bool(ok)})

    def block(self):
        n = sum(1 for c in self.items if c["pass"])
        return {
            "checks": self.items,
            "n_pass": n,
            "n_total": len(self.items),
            "all_pass": n == len(self.items),
            "note": (
                "Each row re-derives a reported value along an independent "
                "arithmetic path written from the printed clause, and "
                "compares. A pass means this report is internally consistent "
                "with the provision as printed. It is NOT a statement about "
                "the design, and it is not a substitute for the textbook "
                "benchmarks in the test suite -- an implementation and a "
                "re-derivation can agree and both be wrong about the Code."
            ),
        }


def _adv(code, severity, text, clause=None):
    return {"code": code, "severity": severity, "text": text, "clause": clause}


def _assemble(title, basis, req, sections, summary, checks, advisories,
              governing, unavailable, has_demand):
    """Common envelope.

    Two verdict fields, and they mean different things on purpose.

    ``adequate``
        Every check this sheet ACTUALLY COMPUTED is satisfied. None when no
        demand was supplied, because a capacity with nothing to compare it
        against has no verdict and this module will not invent one.

    ``complete``
        Every check the provision requires was actually performed. It is
        ``False`` whenever ``unavailable`` is non-empty, which for these
        three sheets is always -- ``As,min``, joint transverse
        reinforcement and the rest have no model here.

    Collapsing the two into one boolean is what makes a report lie. An
    ``adequate: true`` on a section whose ``As,min`` was never checked is
    true about the arithmetic and false about the section, so the caller
    gets both fields and the renderer shows both.
    """
    return {
        "provenance": provenance(),
        "title": title,
        "basis": basis,
        "request_echo": req,
        "sections": sections,
        "summary": summary,
        "qaqc": checks.block(),
        "advisories": advisories,
        "governing_checks": governing,
        "unavailable": unavailable,
        "adequate": None if not has_demand else not governing,
        "complete": not unavailable,
        "has_demand": has_demand,
    }


# ══════════════════════════════════════════════════════════════════════
# 1 — Beam flexure
# ══════════════════════════════════════════════════════════════════════

def beam_flexure_report(rebar_list, fc, fy, b, h, es=200000.0, mu_demand=None):
    """Calculation sheet for the flexural capacity of an RC beam section.

    Parameters
    ----------
    rebar_list : list of dict
        ``d`` (mm from the compression face), ``diam`` (mm), ``num``.
    fc, fy : float
        Cylinder strength and bar yield, MPa.
    b, h : float
        Width and overall depth, mm.
    es : float
        Steel modulus, MPa.
    mu_demand : float or None
        Factored moment, kN-m. **Optional.** With no demand the sheet
        reports capacity and ``adequate`` is None -- it does not invent a
        verdict.

    Returns
    -------
    dict
        ``provenance``, ``request_echo``, ``sections``, ``summary``,
        ``qaqc``, ``advisories``, ``governing_checks``, ``unavailable``,
        ``adequate``.
    """
    req = {"rebar_list": rebar_list, "fc": fc, "fy": fy, "b": b, "h": h,
           "es": es, "mu_demand": mu_demand}
    r = calculate_beam_moment(rebar_list=rebar_list, fc=fc, fy=fy, b=b, h=h,
                              es=es)

    as_total = sum(math.pi / 4 * bar["diam"] ** 2 * bar["num"]
                   for bar in rebar_list)
    d_max = max(bar["d"] for bar in rebar_list)
    beta1 = r["beta1"]
    c = r["neutral_axis"]
    a = r["a"]
    eps_ty = fy / es
    dc = _dc(mu_demand, r["mu"])

    # ── the sheet ──────────────────────────────────────────────────────
    inputs = [
        _step("--", "Concrete cylinder strength",
              r"\( f'_c = %s \text{ MPa} \)" % _r(fc, 1), "%s MPa" % _r(fc, 1)),
        _step("--", "Reinforcement yield strength",
              r"\( f_y = %s \text{ MPa} \)" % _r(fy, 1), "%s MPa" % _r(fy, 1)),
        _step("--", "Section width and overall depth",
              r"\( b = %s,\ h = %s \text{ mm} \)" % (_r(b, 1), _r(h, 1)),
              "%s &times; %s mm" % (_r(b, 1), _r(h, 1))),
        _step("--", "Total steel area and extreme tension bar depth",
              r"\( A_s = %s \text{ mm}^2,\ d_t = %s \text{ mm} \)"
              % (_r(as_total, 1), _r(d_max, 1)),
              "%s mm&sup2;" % _r(as_total, 1)),
    ]

    materials = [
        _step("ACI 318-25M Table 22.2.2.4.3, printed 438-439",
              "Whitney stress-block factor. 0.85 up to f'c = 28 MPa, then "
              "reducing 0.05 per 7 MPa, with a floor of 0.65 at f'c &ge; 55.",
              r"\( \beta_1 = %s \)" % _r(beta1, 4), "%s" % _r(beta1, 4)),
        _step("ACI 318-25M Section 22.2.2.4.1, printed 438",
              "Equivalent rectangular concrete stress of 0.85 f'c over a "
              "depth a = &beta;&#8321; c.",
              r"\( 0.85 f'_c = %s \text{ MPa} \)" % _r(0.85 * fc, 2),
              "%s MPa" % _r(0.85 * fc, 2)),
        _step("--", "Yield strain of the reinforcement",
              r"\( \varepsilon_{ty} = f_y / E_s = %s / %s = %s \)"
              % (_r(fy, 1), _r(es, 0), _r(eps_ty, 5)), "%s" % _r(eps_ty, 5)),
    ]

    solve = [
        _step("Strain compatibility, force equilibrium",
              "Neutral axis found by iteration so that the concrete "
              "compression plus the compression steel equals the tension "
              "steel. Non-convergence raises rather than reporting a zero.",
              r"\( c = %s \text{ mm} \)" % _r(c, 2), "c = %s mm" % _r(c, 2)),
        _step("ACI 318-25M Section 22.2.2.4.1, printed 438",
              "Depth of the equivalent rectangular stress block.",
              r"\( a = \beta_1 c = %s \times %s = %s \text{ mm} \)"
              % (_r(beta1, 4), _r(c, 2), _r(a, 2)), "a = %s mm" % _r(a, 2)),
        _step("W&amp;M Ex 4-4 step 5, printed 170",
              "Bars inside depth a displace concrete the stress block has "
              "already counted at 0.85 f'c, so their force is "
              "A's(f's &minus; 0.85 f'c). A bar in compression but below a "
              "sits in cracked concrete and displaces nothing.",
              r"\( C_c = 0.85 f'_c a b = %s \text{ kN} \)"
              % _r(r["fc_concrete"], 1), "Cc = %s kN" % _r(r["fc_concrete"], 1)),
        _step("--", "Net tensile strain at the extreme tension bar",
              r"\( \varepsilon_t = 0.003\,\dfrac{d_t - c}{c} = %s \)"
              % _r(r["epsilon_t"], 5), "%s" % _r(r["epsilon_t"], 5)),
        _step("NSCP 2015 / ACI 318M-14",
              "Strength reduction factor. Tension-controlled at "
              "&epsilon;t &ge; 0.005. THIS IS THE PRE-2019 RULE AND IT IS "
              "THE CORRECT ONE UNDER NSCP 2015. ACI 318-19/-25 Table 21.2.2 "
              "(printed 430) replaced it with &epsilon;ty + 0.003.",
              r"\( \phi = %s \)" % _r(r["phi"], 4),
              "%s <span class='dim'>%s</span>" % (_r(r["phi"], 4),
                                                  r["classification"])),
    ]

    capacity = [
        _step("Moments about the compression-block centroid",
              "Nominal flexural strength.",
              r"\( M_n = %s \text{ kN} \cdot \text{m} \)" % _r(r["mn"], 2),
              "Mn = %s kN&middot;m" % _r(r["mn"], 2)),
        _step("ACI 318-25M Section 21.2.1, printed 429",
              "Design flexural strength.",
              r"\( \phi M_n = %s \times %s = %s \text{ kN} \cdot \text{m} \)"
              % (_r(r["phi"], 4), _r(r["mn"], 2), _r(r["mu"], 2)),
              "&phi;Mn = %s kN&middot;m" % _r(r["mu"], 2)),
    ]

    if mu_demand is not None:
        capacity.append(_step(
            "ACI 318-25M Section 9.5.1.1, printed 148",
            "Flexural demand/capacity ratio.",
            r"\( M_u / \phi M_n = %s / %s = %s \)"
            % (_r(mu_demand, 2), _r(r["mu"], 2), _r(dc, 3)),
            "D/C = %s" % _r(dc, 3),
            "ok" if _verdict(dc) else "fail"))

    sections = [
        {"heading": "1. Design inputs", "steps": inputs},
        {"heading": "2. Material constants", "steps": materials},
        {"heading": "3. Neutral axis and strain state", "steps": solve},
        {"heading": "4. Flexural capacity", "steps": capacity},
    ]

    summary = [
        {"label": "Neutral axis, c", "value": "%s mm" % _r(c, 2), "note": ""},
        {"label": "Stress-block depth, a", "value": "%s mm" % _r(a, 2),
         "note": "&beta;&#8321; = %s" % _r(beta1, 4)},
        {"label": "Net tensile strain, &epsilon;t",
         "value": "%s" % _r(r["epsilon_t"], 5), "note": r["classification"]},
        {"label": "Nominal moment, Mn",
         "value": "%s kN&middot;m" % _r(r["mn"], 2), "note": ""},
        {"label": "Design moment, &phi;Mn",
         "value": "%s kN&middot;m" % _r(r["mu"], 2),
         "note": "&phi; = %s" % _r(r["phi"], 4)},
        {"label": "Flexural D/C",
         "value": "not evaluated" if dc is None else "%s" % _r(dc, 3),
         "note": "no demand supplied" if dc is None else "Mu = %s kN&middot;m"
                 % _r(mu_demand, 2)},
    ]

    # ── QAQC ───────────────────────────────────────────────────────────
    checks = _Checks()
    beta1_hand = (0.85 if fc <= 28
                  else (0.65 if fc >= 55
                        else max(0.85 - 0.05 * (fc - 28) / 7, 0.65)))
    checks.add("beta1", "Table 22.2.2.4.3 longhand from f'c",
               beta1_hand, beta1)
    checks.add("Stress-block depth a", "beta1 * c from the reported c",
               beta1 * c, a, rel=1e-3)
    checks.add("Concrete compression Cc",
               "0.85 f'c a b, longhand, in kN",
               0.85 * fc * a * b / N_PER_KN, r["fc_concrete"], rel=1e-3)
    checks.add("Net tensile strain eps_t",
               "0.003 (d_t - c)/c from the reported c",
               0.003 * (d_max - c) / c, r["epsilon_t"], rel=1e-3)
    checks.add("Yield strain eps_ty",
               "fy / Es, longhand from the raw material inputs",
               eps_ty, r["epsilon_ty"], rel=1e-6)
    # phi from the pre-2019 law, longhand
    if r["epsilon_t"] >= 0.005:
        phi_hand = 0.90
    elif r["epsilon_t"] <= eps_ty:
        phi_hand = 0.65
    else:
        phi_hand = 0.65 + 0.25 * (r["epsilon_t"] - eps_ty) / (0.005 - eps_ty)
    checks.add("Strength reduction factor phi",
               "NSCP 2015 / ACI 318M-14 pre-2019 law, longhand from eps_t",
               phi_hand, r["phi"], rel=1e-3)
    checks.add("Design moment phi*Mn", "phi * Mn from the two reported values",
               r["phi"] * r["mn"], r["mu"], rel=1e-3)
    # equilibrium: sum of reported bar forces must equal the concrete force
    f_net = sum(rf["force"] for rf in r["rebar_forces"]) / N_PER_KN
    checks.add("Section equilibrium", "sum of the reported bar forces vs Cc",
               r["fc_concrete"], f_net, rel=2e-3, abs_=1.0)
    checks.add("Total steel area As", "sum of pi/4 d^2 n over the input bars",
               as_total, sum(rf["area"] for rf in r["rebar_forces"]), rel=1e-6)
    if mu_demand is not None:
        checks.add("Flexural D/C", "Mu / phi*Mn from the two reported values",
                   mu_demand / r["mu"], dc, rel=1e-9)

    # ── advisories, governing, unavailable ─────────────────────────────
    advisories = [
        _adv("F-EDITION", "info",
             "This sheet is NSCP 2015 (= ACI 318M-14), not ACI 318-19. The "
             "phi law used is tension-controlled at eps_t >= 0.005. Under "
             "ACI 318-19/-25 the limit is eps_ty + 0.003, which for Grade "
             "550 and 690 bars differs by +7.5% and +15.5% on phi. Confirm "
             "the governing edition for the project before filing this.",
             "ACI 318-25M Table 21.2.2, printed 430"),
        _adv("F-ASMIN", "critical",
             "As,min is NOT checked by this sheet. A section can report a "
             "flexural capacity here and still be inadmissible.",
             "ACI 318-25M Section 9.6.1.2, printed 149"),
        _adv("F-RHO", "warning",
             "The rho <= 0.025 cap and the SMF relations Mn+ >= 0.5 Mn- at a "
             "joint face are NOT checked. If this beam is part of a special "
             "moment frame, those are separate verifications.",
             "ACI 318-25M Sections 18.6.3.1 and 18.6.3.2"),
        _adv("F-DETAIL", "warning",
             "Bar spacing, cover, and development are NOT checked. The bar "
             "depths supplied are taken as given; nothing verifies that the "
             "bars fit in the width, or that a layer is developed.",
             "ACI 318-25M Chapter 25"),
    ]
    if r["classification"] != "tension-controlled":
        advisories.insert(0, _adv(
            "F-BRITTLE", "critical",
            "The section is %s (eps_t = %s). ACI 318-25M Section 9.3.3.1 "
            "requires eps_t >= eps_ty + 0.003 for beams -- a compression-"
            "controlled beam is not merely penalised by phi, it is "
            "prohibited. This sheet applies the phi penalty and does NOT "
            "apply the prohibition."
            % (r["classification"], _r(r["epsilon_t"], 5)),
            "ACI 318-25M Section 9.3.3.1"))

    unavailable = [
        {"check": "As,min", "why": "not implemented",
         "clause": "ACI 318-25M Section 9.6.1.2, printed 149"},
        {"check": "rho <= 0.025 and Mn+ >= 0.5 Mn-", "why": "not implemented",
         "clause": "ACI 318-25M Sections 18.6.3.1 / 18.6.3.2"},
        {"check": "bar spacing, cover, development", "why": "not implemented",
         "clause": "ACI 318-25M Chapter 25"},
    ]

    governing = []
    if dc is not None and dc > 1.0:
        governing.append({
            "check": "Flexural strength", "dc": _r(dc, 3),
            "detail": "Mu = %s kN.m exceeds phi*Mn = %s kN.m"
                      % (_r(mu_demand, 2), _r(r["mu"], 2)),
            "clause": "ACI 318-25M Section 9.5.1.1"})

    out = _assemble(
        "Beam Flexural Capacity",
        "NSCP 2015 Section 422.2 (= ACI 318M-14)",
        req, sections, summary, checks, advisories, governing, unavailable,
        has_demand=mu_demand is not None)
    # The solver's iteration log is several hundred rows and the sheet uses
    # none of it. Shipping it put ~40 kB of neutral-axis trial values into
    # every response, and carried the string "Infinity" for the ratio column
    # of every pre-convergence row -- valid JSON, but a token no consumer
    # should have to reason about. Kept out of the report; the calculator
    # endpoint that wants it still returns it.
    out["raw"] = {k: v for k, v in r.items() if k != "iterations"}
    return out


# ══════════════════════════════════════════════════════════════════════
# 2 — Beam shear
# ══════════════════════════════════════════════════════════════════════

def beam_shear_report(fc, fyv, bw, h, d, vu, s_chosen, n_legs, db_stirrup,
                      cc=40.0, c=65.0, nu=0.0, phi=0.75):
    """Calculation sheet for beam shear design.

    Parameters
    ----------
    fc, fyv : float   Cylinder strength and stirrup yield, MPa.
    bw, h, d : float  Web width, overall depth, effective depth, mm.
    vu : float        Factored shear, kN. Required -- shear design without a
                      demand is meaningless, unlike a flexural capacity.
    s_chosen : float  Trial stirrup spacing, mm.
    n_legs : int      Stirrup legs crossing the crack.
    db_stirrup : float  Stirrup bar diameter, mm.
    cc, c : float     Clear cover to links, and cover to the bar centroid, mm.
    nu : float        Factored axial force, kN, positive = compression.
    phi : float       Strength reduction factor for shear.
    """
    req = {"fc": fc, "fyv": fyv, "bw": bw, "h": h, "d": d, "vu": vu,
           "s_chosen": s_chosen, "n_legs": n_legs, "db_stirrup": db_stirrup,
           "cc": cc, "c": c, "nu": nu, "phi": phi}
    sd = shear_design(fc=fc, fyv=fyv, phi=phi, bw=bw, h=h, cc=cc, c=c, d=d,
                      vu=vu, nu=nu, s_chosen=s_chosen, n_legs=n_legs,
                      db_stirrup=db_stirrup)

    av_leg = math.pi / 4 * db_stirrup ** 2
    av = n_legs * av_leg
    vs_provided = av * fyv * d / s_chosen / N_PER_KN
    phi_vn = phi * (sd["vc"] + min(vs_provided, sd["vs_max"]))
    dc = _dc(vu, phi_vn)
    ag = bw * h

    inputs = [
        _step("--", "Concrete and stirrup material strengths",
              r"\( f'_c = %s,\ f_{yt} = %s \text{ MPa} \)"
              % (_r(fc, 1), _r(fyv, 1)),
              "%s / %s MPa" % (_r(fc, 1), _r(fyv, 1))),
        _step("--", "Web width, overall depth, effective depth",
              r"\( b_w = %s,\ h = %s,\ d = %s \text{ mm} \)"
              % (_r(bw, 1), _r(h, 1), _r(d, 1)),
              "%s &times; %s, d = %s mm" % (_r(bw, 1), _r(h, 1), _r(d, 1))),
        _step("--", "Factored shear and axial force at the section",
              r"\( V_u = %s \text{ kN},\ N_u = %s \text{ kN} \)"
              % (_r(vu, 2), _r(nu, 2)),
              "Vu = %s kN" % _r(vu, 2)),
        _step("--", "Trial stirrups",
              r"\( %d \text{ legs} \times \varnothing %s \text{ at } s = %s "
              r"\text{ mm},\ A_v = %s \text{ mm}^2 \)"
              % (n_legs, _r(db_stirrup, 1), _r(s_chosen, 1), _r(av, 1)),
              "Av = %s mm&sup2;" % _r(av, 1)),
    ]

    if nu > 0:
        vc_ref = "W&amp;M Eq. (6-13aM), printed 282 (ACI Section 22.5.6.1)"
        vc_desc = ("Concrete shear strength with axial compression. There is "
                   "NO upper cap on this expression -- an unsupported "
                   "min() branch was removed, and it never governed for any "
                   "Nu &ge; 0 in any case.")
        vc_eq = (r"\( V_c = \left(1 + \dfrac{N_u}{14 A_g}\right)"
                 r"\dfrac{\sqrt{f'_c}}{6} b_w d = %s \text{ kN} \)"
                 % _r(sd["vc"], 2))
    elif nu < 0:
        vc_ref = "ACI 318M-14 Section 22.5.7.1"
        vc_desc = ("Concrete shear strength with net axial tension, floored "
                   "at zero.")
        vc_eq = (r"\( V_c = \left(1 + \dfrac{0.29 N_u}{A_g}\right)"
                 r"\dfrac{\sqrt{f'_c}}{6} b_w d = %s \text{ kN} \)"
                 % _r(sd["vc"], 2))
    else:
        vc_ref = "NSCP 2015 Section 422.5 (= ACI 318M-14)"
        vc_desc = ("Concrete shear strength, no axial load. This is the "
                   "318-14 law. ACI 318-19/-25 Table 22.5.5.1 restructured "
                   "it and adds a size-effect factor lambda_s for sections "
                   "with less than Av,min -- deliberately absent here.")
        vc_eq = (r"\( V_c = \dfrac{\sqrt{f'_c}}{6} b_w d = %s \text{ kN} \)"
                 % _r(sd["vc"], 2))

    concrete = [_step(vc_ref, vc_desc, vc_eq, "Vc = %s kN" % _r(sd["vc"], 2))]

    steel = [
        _step("ACI 318-25M Section 22.5.10.5.3",
              "Shear the stirrups must carry.",
              r"\( V_{s,req} = \dfrac{V_u - \phi V_c}{\phi} = %s \text{ kN} \)"
              % _r(sd["vs_req"], 2), "Vs,req = %s kN" % _r(sd["vs_req"], 2)),
        _step("ACI 318-25M Section 22.5.1.2, printed 442",
              "Cross-sectional limit. Above this the section itself is "
              "inadequate and no stirrup spacing rescues it.",
              r"\( V_{s,max} = \tfrac{2}{3}\sqrt{f'_c}\,b_w d = %s "
              r"\text{ kN} \)" % _r(sd["vs_max"], 2),
              "Vs,max = %s kN" % _r(sd["vs_max"], 2)),
        _step("ACI 318-25M Section 22.5.1.2, printed 442",
              "Governing upper bound on the factored shear.",
              r"\( \phi(V_c + V_{s,max}) = %s \text{ kN} \)"
              % _r(sd["vu_max"], 2),
              "%s <span class='dim'>Vu,max = %s kN</span>"
              % (sd["shear_status"], _r(sd["vu_max"], 2)),
              "ok" if sd["shear_status"] == "SAFE" else "fail"),
        _step("ACI 318-25M Section 22.5.10.5.3",
              "Transverse steel required per unit length for strength.",
              r"\( (A_v/s)_{req} = \dfrac{V_u - \phi V_c}{\phi f_{yt} d} = %s "
              r"\text{ mm}^2/\text{mm} \)" % _r(sd["av_s_req"], 4),
              "%s mm&sup2;/mm" % _r(sd["av_s_req"], 4)),
        _step("ACI 318-25M Table 9.6.3.4, printed 151",
              "Minimum transverse steel, the greater of expressions (a) and "
              "(b). One constant, 0.062, is used everywhere in this package.",
              r"\( (A_v/s)_{min} = \max\left(0.062\sqrt{f'_c}\dfrac{b_w}"
              r"{f_{yt}},\ 0.35\dfrac{b_w}{f_{yt}}\right) = %s \)"
              % _r(sd["av_min"], 4), "%s mm&sup2;/mm" % _r(sd["av_min"], 4)),
        _step("--", "Governing transverse steel requirement.",
              r"\( (A_v/s)_{gov} = \max = %s \text{ mm}^2/\text{mm} \)"
              % _r(sd["av_s_govern"], 4),
              "%s mm&sup2;/mm" % _r(sd["av_s_govern"], 4)),
    ]

    detailing = [
        _step("ACI 318-25M Table 9.7.6.2.2, printed 160",
              "Maximum stirrup spacing, halved once Vs exceeds "
              "(1/3)&radic;f'c bw d.",
              r"\( s_{max} = %s \text{ mm} \)" % _r(sd["smax"], 0),
              "s,max = %s mm" % _r(sd["smax"], 0)),
        _step("--", "Provided transverse steel at the trial spacing.",
              r"\( A_v/s = \dfrac{%d \times %s}{%s} = %s \text{ mm}^2"
              r"/\text{mm} \)"
              % (n_legs, _r(av_leg, 1), _r(s_chosen, 1),
                 _r(sd["av_provided"], 4)),
              "%s mm&sup2;/mm &mdash; %s"
              % (_r(sd["av_provided"], 4), sd["rft_ok"]),
              "ok" if sd["rft_ok"] == "OK" else "fail"),
        _step("--", "Trial spacing against the maximum.",
              r"\( s = %s \le s_{max} = %s \text{ mm} \)"
              % (_r(s_chosen, 1), _r(sd["smax"], 0)),
              sd["spacing_ok"],
              "ok" if sd["spacing_ok"] == "OK" else "fail"),
    ]

    cap = [
        _step("ACI 318-25M Section 22.5.1.1",
              "Design shear strength at the trial spacing, with Vs capped at "
              "the Section 22.5.1.2 limit.",
              r"\( \phi V_n = \phi (V_c + V_s) = %s \text{ kN} \)"
              % _r(phi_vn, 2), "&phi;Vn = %s kN" % _r(phi_vn, 2)),
        _step("ACI 318-25M Section 22.5.1.1",
              "Shear demand/capacity ratio.",
              r"\( V_u / \phi V_n = %s / %s = %s \)"
              % (_r(vu, 2), _r(phi_vn, 2), _r(dc, 3)),
              "D/C = %s" % _r(dc, 3), "ok" if _verdict(dc) else "fail"),
    ]

    sections = [
        {"heading": "1. Design inputs", "steps": inputs},
        {"heading": "2. Concrete shear strength", "steps": concrete},
        {"heading": "3. Transverse reinforcement", "steps": steel},
        {"heading": "4. Detailing limits", "steps": detailing},
        {"heading": "5. Design shear strength", "steps": cap},
    ]

    summary = [
        {"label": "Concrete, Vc", "value": "%s kN" % _r(sd["vc"], 2),
         "note": sd["vc_note"]},
        {"label": "Required, Vs,req", "value": "%s kN" % _r(sd["vs_req"], 2),
         "note": "limit %s kN" % _r(sd["vs_max"], 2)},
        {"label": "Governing Av/s",
         "value": "%s mm&sup2;/mm" % _r(sd["av_s_govern"], 4),
         "note": "Av,min = %s" % _r(sd["av_min"], 4)},
        {"label": "Provided Av/s",
         "value": "%s mm&sup2;/mm" % _r(sd["av_provided"], 4),
         "note": "%d legs &varnothing;%s at %s mm"
                 % (n_legs, _r(db_stirrup, 1), _r(s_chosen, 1))},
        {"label": "Maximum spacing", "value": "%s mm" % _r(sd["smax"], 0),
         "note": sd["spacing_ok"]},
        {"label": "Design shear, &phi;Vn", "value": "%s kN" % _r(phi_vn, 2),
         "note": "&phi; = %s" % _r(phi, 2)},
        {"label": "Shear D/C", "value": "%s" % _r(dc, 3),
         "note": "Vu = %s kN" % _r(vu, 2)},
    ]

    # ── QAQC ───────────────────────────────────────────────────────────
    checks = _Checks()
    if nu > 0:
        vc_hand = ((1 + nu * N_PER_KN / (14 * ag)) * math.sqrt(fc) / 6
                   * bw * d / N_PER_KN)
    elif nu < 0:
        vc_hand = max(0.0, (1 + 0.29 * nu * N_PER_KN / ag)
                      * math.sqrt(fc) / 6 * bw * d / N_PER_KN)
    else:
        vc_hand = math.sqrt(fc) / 6 * bw * d / N_PER_KN
    checks.add("Concrete shear Vc", "Section 22.5.5/22.5.6/22.5.7 longhand",
               vc_hand, sd["vc"], rel=1e-3)
    checks.add("Vs required", "(Vu - phi Vc)/phi from the reported Vc",
               max((vu - phi * sd["vc"]) / phi, 0.0), sd["vs_req"], rel=1e-3)
    checks.add("Vs,max", "(2/3) sqrt(f'c) bw d, Section 22.5.1.2 longhand",
               2 / 3 * math.sqrt(fc) * bw * d / N_PER_KN, sd["vs_max"],
               rel=1e-3)
    checks.add("Vu,max", "phi (Vc + Vs,max) from the two reported values",
               phi * (sd["vc"] + sd["vs_max"]), sd["vu_max"], rel=1e-3)
    checks.add("Av,min per unit length",
               "max(0.062 sqrt(f'c) bw/fyt, 0.35 bw/fyt), Table 9.6.3.4",
               max(0.062 * math.sqrt(fc) * bw / fyv, 0.35 * bw / fyv),
               sd["av_min"], rel=1e-3)
    av_s_req_hand = ((vu * N_PER_KN - phi * sd["vc"] * N_PER_KN)
                     / (phi * fyv * d)) if sd["vs_req"] > 0 else 0.0
    checks.add("Av/s required", "(Vu - phi Vc)/(phi fyt d) longhand",
               max(av_s_req_hand, 0.0), sd["av_s_req"], rel=1e-3)
    checks.add("Av/s governing",
               "greater of the strength requirement and Table 9.6.3.4 minimum",
               max(sd["av_s_req"], sd["av_min"]), sd["av_s_govern"], rel=1e-6)
    checks.add("Av provided per unit length",
               "n_legs * pi/4 * db^2 / s, longhand from the raw inputs",
               av / s_chosen, sd["av_provided"], rel=1e-3)
    smax_hand = (min(300.0, d / 4)
                 if sd["vs"] > (1 / 3) * math.sqrt(fc) * bw * d / N_PER_KN
                 else min(600.0, d / 2))
    # abs_ = 0.51 mm, not a relative tolerance: shear_design rounds s,max to
    # whole millimetres, and round-half-even turns an exact 152.5 into 152.
    # A relative tolerance flags that as a defect. Half a millimetre is the
    # rounding granularity itself, so it is the correct band -- anything
    # wider would stop the check from catching a real error.
    checks.add("Maximum spacing s,max",
               "Table 9.7.6.2.2 longhand, halved above (1/3)sqrt(f'c) bw d; "
               "compared to the nearest whole mm, the reported precision",
               smax_hand, sd["smax"], rel=0.0, abs_=0.51)
    checks.add("Design shear phi*Vn",
               "phi (Vc + min(Vs provided, Vs,max)) from reported values",
               phi * (sd["vc"] + min(vs_provided, sd["vs_max"])), phi_vn,
               rel=1e-9)
    checks.add("Shear D/C", "Vu / phi*Vn from the two reported values",
               vu / phi_vn if phi_vn > 0 else None, dc, rel=1e-9)

    # ── advisories, governing, unavailable ─────────────────────────────
    advisories = [
        _adv("V-EDITION", "info",
             "This sheet is NSCP 2015 (= ACI 318M-14). Vc has no size-effect "
             "factor lambda_s. Under ACI 318-19/-25 a member with less than "
             "Av,min and d = 900 mm carries lambda_s = 0.66 -- a third off "
             "Vc. Confirm the governing edition before filing this.",
             "ACI 318-25M Table 22.5.5.1 and Section 22.5.5.1.3, printed 444-445"),
        _adv("V-SHEARONLY", "warning",
             "TORSION IS NOT CONSIDERED HERE. If the member carries a "
             "torsional moment at or above the threshold, the stirrups must "
             "carry Av + 2At and the longitudinal torsion steel is a "
             "separate requirement. Use the combined shear+torsion sheet.",
             "ACI 318-25M Section 22.7"),
        _adv("V-DEEP", "warning",
             "Deep-beam provisions are NOT applied. If the clear span is "
             "less than four times the effective depth, or a concentrated "
             "load sits within twice the depth of a support, this sheet is "
             "the wrong model.",
             "ACI 318-25M Section 9.9"),
        _adv("V-SMF", "warning",
             "The special-moment-frame capacity-design shear is NOT computed "
             "here. Vu must be supplied by the caller; if this is an SMF "
             "beam, Vu must come from the probable moment strengths at both "
             "ends, and Vc must be taken as zero in the hinge region.",
             "ACI 318-25M Sections 18.6.5.1 and 18.6.5.2"),
    ]
    if nu < 0:
        advisories.insert(0, _adv(
            "V-TENSION", "critical",
            "The section carries net axial tension (Nu = %s kN). Vc is "
            "reduced accordingly, but confirm the tension is real and "
            "factored -- this is the branch where an optimistic sign "
            "convention becomes unconservative." % _r(nu, 2),
            "ACI 318M-14 Section 22.5.7.1"))

    unavailable = [
        {"check": "torsion", "why": "not part of this sheet",
         "clause": "ACI 318-25M Section 22.7"},
        {"check": "deep-beam provisions", "why": "not implemented",
         "clause": "ACI 318-25M Section 9.9"},
        {"check": "shear friction and shear at discontinuities",
         "why": "not implemented", "clause": "ACI 318-25M Sections 22.9, 16.4"},
    ]

    governing = []
    if sd["shear_status"] != "SAFE":
        governing.append({
            "check": "Section 22.5.1.2 cross-sectional limit",
            "dc": _r(vu / sd["vu_max"], 3) if sd["vu_max"] else None,
            "detail": "Vu = %s kN exceeds phi(Vc + Vs,max) = %s kN. Enlarge "
                      "the section." % (_r(vu, 2), _r(sd["vu_max"], 2)),
            "clause": "ACI 318-25M Section 22.5.1.2, printed 442"})
    if sd["rft_ok"] != "OK":
        governing.append({
            "check": "Transverse reinforcement provided",
            "dc": _r(sd["av_s_govern"] / sd["av_provided"], 3)
                  if sd["av_provided"] else None,
            "detail": "Provided Av/s = %s is below the governing %s mm2/mm"
                      % (_r(sd["av_provided"], 4), _r(sd["av_s_govern"], 4)),
            "clause": "ACI 318-25M Section 22.5.10.5.3 / Table 9.6.3.4"})
    if sd["spacing_ok"] != "OK":
        governing.append({
            "check": "Stirrup spacing",
            "dc": _r(s_chosen / sd["smax"], 3) if sd["smax"] else None,
            "detail": "s = %s mm exceeds s,max = %s mm"
                      % (_r(s_chosen, 1), _r(sd["smax"], 0)),
            "clause": "ACI 318-25M Table 9.7.6.2.2, printed 160"})
    if dc is not None and dc > 1.0:
        governing.append({
            "check": "Shear strength", "dc": _r(dc, 3),
            "detail": "Vu = %s kN exceeds phi*Vn = %s kN"
                      % (_r(vu, 2), _r(phi_vn, 2)),
            "clause": "ACI 318-25M Section 22.5.1.1"})

    out = _assemble(
        "Beam Shear Design",
        "NSCP 2015 Section 422.5 (= ACI 318M-14)",
        req, sections, summary, checks, advisories, governing, unavailable,
        has_demand=True)
    out["raw"] = sd
    return out


# ══════════════════════════════════════════════════════════════════════
# 3 — Beam-column joint shear
# ══════════════════════════════════════════════════════════════════════

def joint_shear_report(ve, as1, n_bars1, as2, n_bars2, fy, fc,
                       beam_width, joint_depth, column_width,
                       perpendicular_dist=0.0, joint_config=1, lamda=1.0):
    """Calculation sheet for special-moment-frame joint shear.

    Parameters
    ----------
    ve : float          Column shear consistent with the beam probable
                        moments, kN.
    as1, n_bars1 : float, int   Top bar group: area per bar and count.
    as2, n_bars2 : float, int   Bottom bar group.
    fy, fc : float      Bar yield and cylinder strength, MPa.
    beam_width : float  Beam width b, mm.
    joint_depth : float Column depth h in the direction of joint shear, mm.
    column_width : float  Column width perpendicular to the joint shear, mm.
    perpendicular_dist : float  x, column face to beam face, mm.
    joint_config : int  Row of NSCP 2015 Table 418.8.4.1.
    lamda : float       0.75 lightweight, 1.0 normalweight.

    phi is NOT a parameter. Section 21.2.4.4 fixes it at 0.85.
    """
    req = {"ve": ve, "as1": as1, "n_bars1": n_bars1, "as2": as2,
           "n_bars2": n_bars2, "fy": fy, "fc": fc, "beam_width": beam_width,
           "joint_depth": joint_depth, "column_width": column_width,
           "perpendicular_dist": perpendicular_dist,
           "joint_config": joint_config, "lamda": lamda}
    j = joint_shear_check(
        ve=ve, as1=as1, n_bars1=n_bars1, as2=as2, n_bars2=n_bars2,
        fy=fy, fc=fc, beam_width=beam_width, joint_depth=joint_depth,
        column_width=column_width, perpendicular_dist=perpendicular_dist,
        joint_config=joint_config, lamda=lamda)

    opt_a = beam_width + joint_depth
    opt_b = beam_width + 2 * perpendicular_dist
    gamma = j["gamma"]
    dc = _dc(j["v_joint"], j["phi_vn"])
    config_name = {1: "confined on all four faces",
                   2: "confined on three faces or two opposite faces",
                   3: "other joints"}[joint_config]

    inputs = [
        _step("--", "Concrete and reinforcement strengths",
              r"\( f'_c = %s \text{ MPa},\ f_y = %s \text{ MPa} \)"
              % (_r(fc, 1), _r(fy, 1)),
              "%s / %s MPa" % (_r(fc, 1), _r(fy, 1))),
        _step("--", "Beam width, column depth in the direction of shear, and "
                    "column width perpendicular to it",
              r"\( b = %s,\ h = %s,\ b_{col} = %s \text{ mm} \)"
              % (_r(beam_width, 1), _r(joint_depth, 1), _r(column_width, 1)),
              "b %s / h %s / bcol %s mm"
              % (_r(beam_width, 1), _r(joint_depth, 1), _r(column_width, 1))),
        _step("--", "Beam reinforcement crossing the joint",
              r"\( n_1 = %d \times %s,\ n_2 = %d \times %s \text{ mm}^2 \)"
              % (n_bars1, _r(as1, 1), n_bars2, _r(as2, 1)),
              "%s / %s mm&sup2;"
              % (_r(as1 * n_bars1, 1), _r(as2 * n_bars2, 1))),
        _step("--", "Column shear consistent with the beam probable moments",
              r"\( V_{col} = %s \text{ kN} \)" % _r(ve, 2),
              "Vcol = %s kN" % _r(ve, 2)),
    ]

    demand = [
        _step("ACI 318-25M Section 18.8.2.1, printed 342",
              "Beam bars at a joint are assumed stressed to 1.25 fy, not fy "
              "-- the joint must survive the beam reaching its probable "
              "strength, so the overstrength is deliberate.",
              r"\( T_1 = 1.25 A_{s1} f_y = %s \text{ kN} \)" % _r(j["t1"], 2),
              "T&#8321; = %s kN" % _r(j["t1"], 2)),
        _step("ACI 318-25M Section 18.8.2.1, printed 342",
              "Second bar group, same overstrength. The compression force "
              "on one side equals the tension force on the other.",
              r"\( C_2 = T_2 = 1.25 A_{s2} f_y = %s \text{ kN} \)"
              % _r(j["t2"], 2), "T&#8322; = %s kN" % _r(j["t2"], 2)),
        _step("ACI 318-25M Section 18.8.4.1, printed 342",
              "Joint shear on a horizontal plane at mid-height.",
              r"\( V_u = T_1 + C_2 - V_{col} = %s + %s - %s = %s "
              r"\text{ kN} \)"
              % (_r(j["t1"], 2), _r(j["t2"], 2), _r(ve, 2),
                 _r(j["v_joint"], 2)),
              "Vu,j = %s kN" % _r(j["v_joint"], 2)),
    ]

    area = [
        _step("ACI 318-25M Section 15.5.2.2, printed 232",
              "Joint depth is the overall depth of the column in the "
              "direction of the joint shear considered.",
              r"\( h = %s \text{ mm} \)" % _r(joint_depth, 1),
              "h = %s mm" % _r(joint_depth, 1)),
        _step("ACI 318-25M Section 15.5.2.2(a), printed 232",
              "First limb: beam width plus joint depth.",
              r"\( b + h = %s + %s = %s \text{ mm} \)"
              % (_r(beam_width, 1), _r(joint_depth, 1), _r(opt_a, 1)),
              "%s mm" % _r(opt_a, 1)),
        _step("ACI 318-25M Section 15.5.2.2(b), printed 232",
              "Second limb: twice the perpendicular distance from the beam "
              "axis to the nearest column side face. That is 2(b/2 + x), "
              "which is identically b + 2x with x measured beam face to "
              "column face -- the form used here, and the form Fig. "
              "R15.5.2.2 draws.",
              r"\( b + 2x = %s + 2(%s) = %s \text{ mm} \)"
              % (_r(beam_width, 1), _r(perpendicular_dist, 1), _r(opt_b, 1)),
              "%s mm" % _r(opt_b, 1)),
        _step("ACI 318-25M Section 15.5.2.2 + R15.5.2.2, printed 231-232",
              "Effective joint width is the overall column width where the "
              "beam is wider than the column; otherwise the lesser of (a) "
              "and (b). R15.5.2.2: in no case is Aj greater than the column "
              "cross-sectional area. All three are taken together.",
              r"\( b_j = \min(b_{col},\ b+h,\ b+2x) = \min(%s, %s, %s) = %s "
              r"\text{ mm} \)"
              % (_r(column_width, 1), _r(opt_a, 1), _r(opt_b, 1),
                 _r(j["joint_width"], 1)),
              "bj = %s mm" % _r(j["joint_width"], 1)),
        _step("ACI 318-25M Section 15.5.2.2, printed 232",
              "Effective joint area.",
              r"\( A_j = b_j h = %s \times %s = %s \text{ mm}^2 \)"
              % (_r(j["joint_width"], 1), _r(joint_depth, 1), _r(j["aj"], 0)),
              "Aj = %s mm&sup2;" % _r(j["aj"], 0)),
    ]

    strength = [
        _step("NSCP 2015 Table 418.8.4.1 (= ACI 318M-14)",
              "Joint shear strength coefficient for a %s. THIS IS A "
              "THREE-ROW TABLE. ACI 318-19/-25 Table 18.8.4.3 (printed 342) "
              "replaced it with eight rows, 1.7/1.3/1.3/1.0/1.3/1.0/1.0/0.7 "
              "-- the middle value moved 1.2 to 1.3 and a 0.7 row was added "
              "that this enum cannot express." % config_name,
              r"\( \gamma = %s \)" % _r(gamma, 2), "&gamma; = %s" % _r(gamma, 2)),
        _step("ACI 318-25M Table 15.5.2.1 footnote [1], printed 231",
              "Lightweight modification factor. Two legal values only: 0.75 "
              "for lightweight and 1.0 for normalweight.",
              r"\( \lambda = %s \)" % _r(lamda, 2), "&lambda; = %s" % _r(lamda, 2)),
        _step("NSCP 2015 Table 418.8.4.1 (= ACI 318M-14)",
              "Nominal joint shear strength.",
              r"\( V_n = \gamma \lambda \sqrt{f'_c}\, A_j = %s \text{ kN} \)"
              % _r(j["vn"], 2), "Vn = %s kN" % _r(j["vn"], 2)),
        _step("ACI 318-25M Section 21.2.4.4, printed 435",
              "Strength reduction factor for the joint shear of a special "
              "moment frame. Fixed at 0.85 -- the Code gives no caller "
              "discretion, and this sheet accepts no override.",
              r"\( \phi = %s \)" % _r(PHI_SMF_JOINT, 2),
              "&phi; = %s" % _r(PHI_SMF_JOINT, 2)),
        _step("ACI 318-25M Section 18.8.4.3",
              "Design joint shear strength.",
              r"\( \phi V_n = %s \times %s = %s \text{ kN} \)"
              % (_r(PHI_SMF_JOINT, 2), _r(j["vn"], 2), _r(j["phi_vn"], 2)),
              "&phi;Vn = %s kN" % _r(j["phi_vn"], 2)),
        _step("ACI 318-25M Section 18.8.4",
              "Joint shear demand/capacity ratio.",
              r"\( V_u / \phi V_n = %s / %s = %s \)"
              % (_r(j["v_joint"], 2), _r(j["phi_vn"], 2), _r(dc, 3)),
              "D/C = %s &mdash; %s" % (_r(dc, 3), j["status"]),
              "ok" if _verdict(dc) else "fail"),
    ]

    sections = [
        {"heading": "1. Design inputs", "steps": inputs},
        {"heading": "2. Joint shear demand", "steps": demand},
        {"heading": "3. Effective joint area", "steps": area},
        {"heading": "4. Joint shear strength", "steps": strength},
    ]

    summary = [
        {"label": "T&#8321;", "value": "%s kN" % _r(j["t1"], 2),
         "note": "1.25 fy overstrength"},
        {"label": "T&#8322; = C&#8322;", "value": "%s kN" % _r(j["t2"], 2),
         "note": ""},
        {"label": "Joint shear, Vu",
         "value": "%s kN" % _r(j["v_joint"], 2),
         "note": "less Vcol = %s kN" % _r(ve, 2)},
        {"label": "Effective width, bj",
         "value": "%s mm" % _r(j["joint_width"], 1),
         "note": "min(bcol %s, b+h %s, b+2x %s)"
                 % (_r(column_width, 0), _r(opt_a, 0), _r(opt_b, 0))},
        {"label": "Effective area, Aj",
         "value": "%s mm&sup2;" % _r(j["aj"], 0), "note": ""},
        {"label": "Design strength, &phi;Vn",
         "value": "%s kN" % _r(j["phi_vn"], 2),
         "note": "&gamma; = %s, &phi; = 0.85" % _r(gamma, 2)},
        {"label": "Joint D/C", "value": "%s" % _r(dc, 3),
         "note": j["status"]},
    ]

    # ── QAQC ───────────────────────────────────────────────────────────
    checks = _Checks()
    checks.add("Tension force T1", "1.25 n As fy / 1000, longhand",
               1.25 * n_bars1 * as1 * fy / N_PER_KN, j["t1"], rel=1e-3)
    checks.add("Tension force T2", "1.25 n As fy / 1000, longhand",
               1.25 * n_bars2 * as2 * fy / N_PER_KN, j["t2"], rel=1e-3)
    checks.add("Joint shear Vu", "T1 + T2 - Vcol from the reported forces",
               j["t1"] + j["t2"] - ve, j["v_joint"], rel=1e-3)
    # R15.5.2.2, printed 231-232: "In no case is Aj greater than the column
    # cross-sectional area." Asserted as an invariant on the reported Aj,
    # not recomputed from the same min() that produced it -- this is the
    # check that would have caught the missing column cap.
    checks.add("Aj within the column cross-section",
               "R15.5.2.2 invariant: reported Aj vs column_width * h",
               min(j["aj"], column_width * joint_depth), j["aj"], rel=1e-9)
    checks.add("Effective width bj",
               "min(column width, b+h, b+2x) longhand from the raw inputs",
               min(column_width, opt_a, opt_b), j["joint_width"], rel=1e-9)
    checks.add("Effective area Aj", "bj * h from the reported bj",
               j["joint_width"] * joint_depth, j["aj"], rel=1e-6)
    checks.add("Gamma from the NSCP table",
               "NSCP 2015 Table 418.8.4.1 lookup on joint_config",
               NSCP_2015_TABLE_418_8_4_1[joint_config], gamma, rel=1e-9)
    checks.add("Nominal strength Vn",
               "gamma lambda sqrt(f'c) Aj / 1000, longhand",
               gamma * lamda * math.sqrt(fc) * j["aj"] / N_PER_KN, j["vn"],
               rel=1e-3)
    checks.add("Design strength phi*Vn",
               "0.85 * Vn, Section 21.2.4.4 fixed value",
               PHI_SMF_JOINT * j["vn"], j["phi_vn"], rel=1e-3)
    checks.add("Phi is the fixed code value",
               "Section 21.2.4.4 permits no other value",
               PHI_SMF_JOINT, j["phi"], rel=1e-9)
    checks.add("Joint D/C", "Vu / phi*Vn from the two reported values",
               j["v_joint"] / j["phi_vn"] if j["phi_vn"] > 0 else None, dc,
               rel=1e-9)

    # ── advisories, governing, unavailable ─────────────────────────────
    advisories = [
        _adv("J-EXTERIOR", "critical",
             "THIS SHEET ALWAYS SUMS BOTH BAR GROUPS. At an exterior or "
             "corner joint only one beam frames in, and the caller must pass "
             "zeros for the absent group. Nothing in the inputs detects "
             "this, and a two-sided demand at a one-sided joint is "
             "unconservative by roughly the second group's force.",
             "ACI 318-25M Section 18.8.4.1"),
        _adv("J-NSCPPAGE", "warning",
             "The gamma values 1.7 / 1.2 / 1.0 are NSCP 2015 Table "
             "418.8.4.1, printed 4-118. That page is CARRIED FROM A PRIOR "
             "REVIEW AND HAS NOT BEEN RE-READ -- the NSCP 2015 reference "
             "copy is image-only. The values are corroborated by ACI "
             "318M-14; the page itself is unverified.",
             "NSCP 2015 Table 418.8.4.1, printed 4-118"),
        _adv("J-EDITION", "info",
             "The three-row table above is correct under NSCP 2015. ACI "
             "318-19/-25 Table 18.8.4.3 has eight rows and is keyed on "
             "column continuity x beam continuity x transverse-beam "
             "confinement. A discontinuous unconfined joint scored here as "
             "'other' takes 1.0 where that table gives 0.7 -- 43% "
             "unconservative against the newer edition.",
             "ACI 318-25M Table 18.8.4.3, printed 342"),
        _adv("J-TRANSVERSE", "critical",
             "JOINT TRANSVERSE REINFORCEMENT IS NOT CHECKED. The gamma "
             "values assume at least the minimum confinement is present in "
             "the joint; this sheet verifies shear only.",
             "ACI 318-25M Section 18.8.3"),
        _adv("J-DEVELOPMENT", "warning",
             "Development of the beam bars into the joint is NOT checked, "
             "nor is the column-to-beam flexural strength ratio.",
             "ACI 318-25M Sections 18.8.5 and 18.7.3"),
    ]
    if joint_config == 1 and perpendicular_dist == 0 and \
            beam_width < column_width:
        advisories.append(_adv(
            "J-CONFINE", "warning",
            "Config 1 claims confinement on all four faces, which requires "
            "transverse beams satisfying Section 15.5.2.5 -- covering at "
            "least three quarters of the column face, on both sides. "
            "Confirm the transverse beams actually exist and are deep "
            "enough; the gamma of 1.7 is the largest in the table.",
            "ACI 318-25M Section 15.5.2.5"))

    unavailable = [
        {"check": "joint transverse reinforcement", "why": "not implemented",
         "clause": "ACI 318-25M Section 18.8.3"},
        {"check": "development of beam bars into the joint",
         "why": "not implemented", "clause": "ACI 318-25M Section 18.8.5"},
        {"check": "column-to-beam flexural strength ratio",
         "why": "not implemented", "clause": "ACI 318-25M Section 18.7.3.2"},
    ]

    governing = []
    if dc is not None and dc > 1.0:
        governing.append({
            "check": "Joint shear strength", "dc": _r(dc, 3),
            "detail": "Vu = %s kN exceeds phi*Vn = %s kN. Enlarge the "
                      "joint, raise f'c, or revise the beam reinforcement."
                      % (_r(j["v_joint"], 2), _r(j["phi_vn"], 2)),
            "clause": "ACI 318-25M Section 18.8.4.3"})

    out = _assemble(
        "Beam-Column Joint Shear",
        "NSCP 2015 Section 418.8.4 (= ACI 318M-14) &mdash; Special Moment Frames",
        req, sections, summary, checks, advisories, governing, unavailable,
        has_demand=True)
    out["raw"] = j
    return out
