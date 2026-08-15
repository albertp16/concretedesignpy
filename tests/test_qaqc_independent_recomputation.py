# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
QAQC -- independent re-derivation of every reported value
=========================================================

Same idiom as ``column_jacket_design._selfcheck``: each row re-derives a
value the module reports along an INDEPENDENT arithmetic path, written
straight from the printed clause, and compares.

A pass means the module is internally consistent with the provision AS
PRINTED. It is NOT a statement about any design, and it does not replace
the textbook pins in the other four test files -- an implementation and a
re-derivation can agree and both be wrong about the Code. The textbook
pins test that; this file tests that the number a caller receives is the
number the formula produces.

Run standalone for the tabular sheet:

    python tests/test_qaqc_independent_recomputation.py

Basis pages are ACI 318-25M and Wight & MacGregor 7th, printed folios.
"""

import math

import pytest

from concretedesignpy.calculators.beam_moment import calculate_beam_moment
from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength, compute_steel_shear_strength,
    compute_shear_spacing, shear_design, shear_torsion_design)
from concretedesignpy.calculators.joint_shear import joint_shear_check

IN, KSI, KIP, KIPFT = 25.4, 6.894757, 4.44822, 1.35582
CHECKS = []


def add(name, basis, expected, reported, rel=1e-3):
    """Record one independent re-derivation against one reported value."""
    dev = abs(reported / expected - 1.0) if expected else abs(reported)
    CHECKS.append(dict(name=name, basis=basis, expected=expected,
                       reported=reported, dev=dev, tol=rel, ok=dev <= rel))


# ── flexure ────────────────────────────────────────────────────────────
r = calculate_beam_moment(rebar_list=[{"d": 500, "diam": 25.4692, "num": 3}],
                          fc=20, fy=420, b=250, h=565)
As = math.pi / 4 * 25.4692 ** 2 * 3
a_hand = As * 420 / (0.85 * 20 * 250)
add("flexure a (singly)", "a = As fy / (0.85 f'c b), Whitney block",
    a_hand, r["a"])
add("flexure c = a/beta1", "Table 22.2.2.4.3, f'c 20 -> beta1 0.85",
    a_hand / 0.85, r["neutral_axis"])
add("flexure Mn (singly)", "Mn = As fy (d - a/2)",
    As * 420 * (500 - a_hand / 2) / 1e6, r["mn"])
add("flexure eps_t", "eps_t = 0.003 (d - c)/c",
    0.003 * (500 - a_hand / 0.85) / (a_hand / 0.85), r["epsilon_t"], rel=2e-3)

a9 = 1.00 * IN * IN
db9 = math.sqrt(4 * a9 / math.pi)
r2 = calculate_beam_moment(
    rebar_list=[{"d": 2.5 * IN, "diam": db9, "num": 3},
                {"d": 16.8 * IN, "diam": db9, "num": 6}],
    fc=4 * KSI, fy=60 * KSI, b=12 * IN, h=20 * IN, es=29000 * KSI)
fc_, fy_, b_, Es_ = 4 * KSI, 60 * KSI, 12 * IN, 29000 * KSI
As_, Asp_, d_, dp_ = 6 * a9, 3 * a9, 16.8 * IN, 2.5 * IN


def resid(c):
    a = 0.85 * c
    fsp = min(fy_, max(-fy_, 0.003 * (c - dp_) / c * Es_))
    ded = 0.85 * fc_ if dp_ < a else 0.0
    return (0.85 * fc_ * a * b_ + Asp_ * (fsp - ded)
            - As_ * min(fy_, max(-fy_, 0.003 * (d_ - c) / c * Es_)))


lo, hi = 1e-3, 20 * IN
for _ in range(400):
    m = (lo + hi) / 2
    lo, hi = (m, hi) if resid(m) < 0 else (lo, m)
c_hand = (lo + hi) / 2
add("flexure c (doubly, Cs deducted)",
    "bisection on Cc + A's(f's - 0.85f'c) - T = 0, W&M Ex 4-4 step 5 p170",
    c_hand, r2["neutral_axis"], rel=2e-3)
a_h = 0.85 * c_hand
fsp_h = min(fy_, max(-fy_, 0.003 * (c_hand - dp_) / c_hand * Es_))
Mn_h = (0.85 * fc_ * a_h * b_ * (d_ - a_h / 2)
        + Asp_ * (fsp_h - 0.85 * fc_) * (d_ - dp_)) / 1e6
add("flexure Mn (doubly)", "moments about the tension steel, independent axis",
    Mn_h, r2["mn"], rel=3e-3)

# ── shear ──────────────────────────────────────────────────────────────
vc = compute_concrete_shear_strength(fc=25, b=300, d=610)
add("Vc simple", "Vc = lam sqrt(f'c) bw d / 6, NSCP 422.5 / 318M-14",
    math.sqrt(25) / 6 * 300 * 610 / 1000, vc["vc_kn"])
vs = compute_steel_shear_strength(av=142, fyt=300, d=610, s=300)
add("Vs", "Vs = Av fyt d / s", 142 * 300 * 610 / 300 / 1000, vs["vs_kn"])

sp = compute_shear_spacing(fc=25, b=300, d=610, fyt=300,
                           vu_required=373e3 * 0.75, phi=0.75, av=258)
vs_req = 373e3 - math.sqrt(25) / 6 * 300 * 610
add("spacing Vs,req", "(Vu/phi) - Vc", vs_req, sp["vs_required"])
s_str = 258 * 300 * 610 / vs_req
s_avm = 258 / max(0.062 * math.sqrt(25) * 300 / 300, 0.35 * 300 / 300)
add("spacing s_avmin", "Av / max(0.062 sqrt(f'c) bw/fyt, 0.35 bw/fyt), T9.6.3.4",
    s_avm, sp["s_avmin"])
add("spacing governing s", "min(strength, Av,min, s_max) -- all three",
    min(s_str, s_avm, min(610 / 2, 600)), sp["spacing"])
add("spacing Vs,max", "(2/3) sqrt(f'c) bw d, 22.5.1.2",
    2 / 3 * math.sqrt(25) * 300 * 610, sp["vs_max"])

sd = shear_design(fc=28, fyv=400, phi=0.75, bw=350, h=600, cc=40, c=70, d=530,
                  vu=300, nu=500, s_chosen=150, n_legs=2, db_stirrup=10)
ag = 350 * 600
add("Vc with axial compression",
    "(1 + Nu/(14 Ag)) sqrt(f'c) bw d / 6, W&M Eq. (6-13aM) p282, no cap",
    (1 + 500e3 / (14 * ag)) * math.sqrt(28) / 6 * 350 * 530 / 1000, sd["vc"])
add("Av,min single constant", "0.062 sqrt(f'c) bw/fyt vs 0.35 bw/fyt, T9.6.3.4",
    max(0.062 * math.sqrt(28) * 350 / 400, 0.35 * 350 / 400), sd["av_min"])

# ── torsion ────────────────────────────────────────────────────────────
bw, h, dd, cc_, ds = 14 * IN, 24 * IN, 21.5 * IN, 1.5 * IN, 0.5 * IN
fc, fy = 4 * KSI, 60 * KSI
t = shear_torsion_design(fc=fc, fyv=fy, fy=fy, phi=0.75, bw=bw, h=h, cc=cc_,
                         c=h - dd, d=dd, vu=57.1 * KIP, tu=28 * KIPFT, nu=0,
                         s_chosen=125, n_legs=2, db_stirrup=ds, db_long=25.4)
aoh_h = (bw - 2 * cc_ - ds) * (h - 2 * cc_ - ds)
ph_h = 2 * ((bw - 2 * cc_ - ds) + (h - 2 * cc_ - ds))
acp, pcp = bw * h, 2 * (bw + h)
add("Aoh (stirrup centreline)", "R22.7.6.1 p465, outermost closed stirrup",
    aoh_h, t["aoh"])
add("ph", "perimeter of the stirrup centreline", ph_h, t["ph"])
add("phi*Tth", "0.083 lam sqrt(f'c) Acp^2/pcp, T22.7.4.1(a) p463",
    0.75 * 0.083 * math.sqrt(fc) * acp ** 2 / pcp / 1e6, t["tth"])
add("phi*Tcr", "0.33 lam sqrt(f'c) Acp^2/pcp, T22.7.5.1 p464",
    0.75 * 0.33 * math.sqrt(fc) * acp ** 2 / pcp / 1e6, t["tcr"])
add("At/s", "Eq. (22.7.6.1a) with Ao = 0.85 Aoh, theta 45: Tu/(phi 1.7 Aoh fyt)",
    28 * KIPFT * 1e6 / (0.75 * 1.7 * aoh_h * fy), t["at_s"])
add("Al", "Eq. (22.7.6.1b): (Tu/phi) ph / (1.7 Aoh fy) cot(theta)",
    28 * KIPFT * 1e6 * ph_h / (0.75 * 1.7 * aoh_h * fy), t["al"])
add("torsional s_max", "lesser of ph/8 and 300, 9.7.6.3.3 p160",
    min(ph_h / 8, 300.0), t["smax"])

tt = shear_torsion_design(fc=fc, fyv=fy, fy=fy, phi=0.75, bw=bw, h=h, cc=cc_,
                          c=h - dd, d=dd, vu=57.1 * KIP, tu=28 * KIPFT,
                          nu=-100.0, s_chosen=125, n_legs=2, db_stirrup=ds,
                          db_long=25.4)
add("phi*Tth under net tension",
    "T22.7.4.1(a) row (c): x sqrt(1 + Nu/(0.33 Ag lam sqrt(f'c))), Nu < 0",
    0.75 * 0.083 * math.sqrt(fc) * acp ** 2 / pcp / 1e6
    * math.sqrt(1 - 100e3 / (0.33 * bw * h * math.sqrt(fc))), tt["tth"])

# ── joint ──────────────────────────────────────────────────────────────
j = joint_shear_check(ve=81.8 * KIP, as1=4.36 * IN * IN, n_bars1=1,
                      as2=2.24 * IN * IN, n_bars2=1, fy=60 * KSI, fc=4 * KSI,
                      beam_width=24 * IN, joint_depth=24 * IN,
                      column_width=24 * IN, perpendicular_dist=0,
                      joint_config=1)
add("joint T1", "1.25 As fy, 18.8.2.1 p342",
    1.25 * 4.36 * IN * IN * 60 * KSI / 1000, j["t1"])
add("joint Vj", "T1 + T2 - Vcol, 18.8.4.1",
    (1.25 * (4.36 + 2.24) * IN * IN * 60 * KSI / 1000) - 81.8 * KIP,
    j["v_joint"])
add("joint Aj", "bj h with bj = min(bcol, b+h, b+2x), 15.5.2.2 p232",
    (24 * IN) ** 2, j["aj"])
add("joint phi*Vn", "0.85 x 1.7 lam sqrt(f'c) Aj, T418.8.4.1 + 21.2.4.4 p435",
    0.85 * 1.7 * math.sqrt(4 * KSI) * (24 * IN) ** 2 / 1000, j["phi_vn"])

j2 = joint_shear_check(ve=0, as1=2000, n_bars1=1, as2=1000, n_bars2=1,
                       fy=420, fc=28, beam_width=500, joint_depth=600,
                       column_width=400, perpendicular_dist=0, joint_config=1)
add("joint Aj capped at the column", "R15.5.2.2 p231-232, spandrel on a slim column",
    400.0 * 600.0, j2["aj"])



# ── the checks are the tests ───────────────────────────────────────────

@pytest.mark.parametrize("check", CHECKS, ids=[c["name"] for c in CHECKS])
def test_reported_value_matches_an_independent_rederivation(check):
    assert check["reported"] == pytest.approx(check["expected"],
                                              rel=check["tol"]), check["basis"]


def test_every_check_carries_a_basis():
    """A QAQC row without a stated basis is not a check."""
    for c in CHECKS:
        assert c["basis"]


def test_the_sheet_covers_all_four_modules():
    names = " ".join(c["name"] for c in CHECKS)
    for topic in ("flexure", "Vc", "Tth", "joint"):
        assert topic in names


if __name__ == "__main__":
    w = max(len(c["name"]) for c in CHECKS)
    print("=" * 100)
    print("QAQC -- independent re-derivation of every reported value")
    print("=" * 100)
    for c in CHECKS:
        print("%s  %-*s  expected %14.5g   reported %14.5g   dev %+7.4f%%"
              % ("PASS" if c["ok"] else "FAIL", w, c["name"],
                 c["expected"], c["reported"], c["dev"] * 100))
        print("      %s" % c["basis"])
    n = sum(1 for c in CHECKS if c["ok"])
    print("-" * 100)
    print("%d of %d checks pass." % (n, len(CHECKS)))
    print("A pass means the module is internally consistent with the provision")
    print("as printed. It is NOT a statement about any design.")
    raise SystemExit(0 if n == len(CHECKS) else 1)
