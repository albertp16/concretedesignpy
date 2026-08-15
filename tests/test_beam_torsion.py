# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam torsion tests -- concretedesignpy.calculators.beam_shear
==============================================================

``shear_torsion_design`` is the package's ONLY implementation of
Section 22.7. The second one, ``beam_torsion.torsion_design``, was
deleted: it computed phi*Tth in N.m and compared it against a Tu the web
form supplied in kN.m, so it answered "MAY NEGLECT" and "PASS" for a beam
4.8x over the threshold.

Governing edition: NSCP 2015 Section 422.7 (= ACI 318M-14).
Benchmark: Wight & MacGregor 7th ed. Printed pages are the printed folio.
W&M Chapter 7 examples are inch-pound; they are converted here, and the
conversion factors are stated at the top of the file.
"""

import math

import pytest

from concretedesignpy.calculators.beam_shear import shear_torsion_design

IN = 25.4                # mm per inch
KSI = 6.894757           # MPa per ksi
KIP = 4.44822            # kN per kip
KIPFT = 1.35582          # kN.m per kip-ft


def rel(actual, expected):
    return abs(actual / expected - 1.0)


def ex_7_2(**over):
    """W&M Example 7-2 (printed 358-362) converted to SI.

    Solid 14 x 24 in. cantilever, d = 21.5 in., cover 1.5 in.,
    No.4 stirrups (0.5 in.), f'c = 4 ksi, fy = 60 ksi,
    Tu = 28 k-ft = 37.96 kN.m, Vu = 57.1 kips = 254.0 kN.
    """
    kwargs = dict(
        fc=4 * KSI, fyv=60 * KSI, fy=60 * KSI, phi=0.75,
        bw=14 * IN, h=24 * IN, cc=1.5 * IN, c=(24 - 21.5) * IN, d=21.5 * IN,
        vu=57.1 * KIP, tu=28.0 * KIPFT, nu=0,
        s_chosen=125, n_legs=2, db_stirrup=0.5 * IN, db_long=25.4,
    )
    kwargs.update(over)
    return shear_torsion_design(**kwargs)


# ---------------------------------------------------------------- textbook

def test_wm_example_7_2_geometry():
    """W&M printed 360: Aoh = 215 in2, ph = 62 in.

    Aoh is bounded by the centreline of the outermost closed transverse
    reinforcement -- the STIRRUP, not the main bar. ACI 318-25M
    R22.7.6.1, printed 465. This function has always used db_stirrup;
    the deleted module used the main-bar diameter and was 7.1% low.
    """
    r = ex_7_2()
    assert rel(r["aoh"], 215 * IN * IN) < 0.01
    assert rel(r["ph"], 62 * IN) < 0.01


def test_wm_example_7_2_threshold_torsion():
    """W&M printed 359: phi*Tth = 5.87 k-ft = 7.96 kN.m, and Tu = 28 k-ft
    is 4.8x that, so torsion must be designed for.

    Basis: ACI 318-25M Table 22.7.4.1(a), printed 463 (read at the page):
    Tth = 0.083 lam sqrt(f'c) (Acp^2/pcp).
    """
    r = ex_7_2()
    assert rel(r["tth"], 5.87 * KIPFT) < 0.01
    assert r["torsion_action"] == "Design for Torsion"


def test_wm_example_7_2_transverse_steel():
    """W&M printed 361: At/s = 0.0204 in2/in per leg."""
    r = ex_7_2()
    assert rel(r["at_s"], 0.0204 * IN * IN / IN) < 0.01


def test_wm_example_7_2_longitudinal_steel():
    """W&M printed 361: Al = 1.27 in2 = 819 mm2.

    Regression pin on P1a/F-9. The divisor is 1.7 Aoh, not 2 Aoh:
    Eq. (22.7.6.1b) is Tn = (2 Ao Al fy / ph) tan(theta) and Section
    22.7.6.1.1 (printed 466) permits Ao = 0.85 Aoh. Using 2 Aoh made
    every Al exactly 0.85 of the required area -- 15.00% short. Before
    the fix this returned 694 mm2.
    """
    r = ex_7_2()
    assert rel(r["al"], 1.27 * IN * IN) < 0.01
    assert r["al"] == pytest.approx(816.2, abs=1.0)
    assert r["al"] > 800.0, "the 2*Aoh divisor has come back (would give 694)"


def test_wm_example_7_2_max_spacing():
    """W&M printed 353 and ACI 318-25M Section 9.7.6.3.3, printed 160:
    torsional stirrup spacing is the lesser of ph/8 and 300 mm. The
    deleted module used 305 (12 in. converted), 1.7% unconservative
    whenever the 305 governed."""
    r = ex_7_2()
    assert r["smax"] == pytest.approx(min(r["ph"] / 8.0, 300.0), abs=1.0)
    assert r["smax"] <= 300.0


def test_wm_example_7_3_hollow_geometry():
    """W&M Example 7-3, printed 363-367. Example 7-2 redesigned as a
    hollow 20 x 24 in. with 5 in. walls.
    Printed 364: Aoh = 16.5 x 20.5 = 338 in2, ph = 74.0 in.
    """
    r = ex_7_2(bw=20 * IN, h=24 * IN)
    assert rel(r["aoh"], 338 * IN * IN) < 0.01
    assert rel(r["ph"], 74 * IN) < 0.01


def test_wm_example_7_3_al_min_governs():
    """W&M printed 366: Al,min = 2.53 - 0.962 = 1.57 in2 governs over the
    Al required by Eq. (22.7.6.1b).

    TOLERANCE IS 2%, NOT 1%, AND HERE IS WHY. Section 9.6.4.3(a),
    printed 152, gives Al,min = 0.42 sqrt(f'c) Acp/fy in SI. W&M works in
    inch-pound, where the same provision reads 5 sqrt(f'c) Acp/fy. The SI
    0.42 is 1.16% above the exact conversion of the inch-pound 5, so this
    example cannot be reproduced inside 1% by a correct SI
    implementation. The gap is code rounding, not a defect -- exactly the
    trap the audit's standing notes warn about.
    """
    r = ex_7_2(bw=20 * IN, h=24 * IN)
    assert rel(r["al"], 1.57 * IN * IN) < 0.02
    # and it is Al,min doing the governing, not the Eq. (22.7.6.1b) area
    al_eq = 28.0 * KIPFT * 1e6 * r["ph"] / (1.7 * 0.75 * r["aoh"] * 60 * KSI)
    assert r["al"] > al_eq, "Al,min must govern in Example 7-3"


# -------------------------------------------------------------- regression

def test_axial_tension_lowers_the_torsion_thresholds():
    """Regression pin on P3. The Nu term of Tables 22.7.4.1(a) and
    22.7.5.1 row (c) was omitted -- Nu was an argument and was never used
    for torsion. Omitting it is conservative under compression and
    UNCONSERVATIVE under net tension.

    Basis: ACI 318-25M Table 22.7.4.1(a), printed 463; Table 22.7.5.1,
    printed 464. Both read at the rendered page. Nu positive for
    compression, negative for tension.
    """
    base = ex_7_2(nu=0)
    tension = ex_7_2(nu=-100.0)
    compression = ex_7_2(nu=+100.0)

    assert tension["tth"] < base["tth"] < compression["tth"]
    assert tension["tcr"] < base["tcr"] < compression["tcr"]

    ag = 14 * IN * 24 * IN
    factor = math.sqrt(1 - 100e3 / (0.33 * ag * math.sqrt(4 * KSI)))
    assert tension["tth"] == pytest.approx(base["tth"] * factor, rel=1e-3)


def test_heavy_net_tension_drives_the_threshold_to_zero():
    """The radicand is clamped at zero rather than raising a domain
    error, so overwhelming tension forces torsion to be designed for."""
    r = ex_7_2(nu=-100000.0)
    assert r["tth"] == 0.0
    assert r["tcr"] == 0.0
    assert r["torsion_action"] == "Design for Torsion"


def test_threshold_uses_the_printed_coefficients():
    """0.083 and 0.33, not 1/12 and 1/3. Tables 22.7.4.1(a) printed 463
    and 22.7.5.1 printed 464. The old 1/12 and 1/3 were +0.4% and +1.0%,
    both in the direction of neglecting torsion."""
    r = ex_7_2()
    bw, h, fc = 14 * IN, 24 * IN, 4 * KSI
    acp, pcp = bw * h, 2 * (bw + h)
    assert r["tth"] == pytest.approx(
        0.75 * 0.083 * math.sqrt(fc) * acp ** 2 / pcp / 1e6, rel=1e-3)
    assert r["tcr"] == pytest.approx(
        0.75 * 0.33 * math.sqrt(fc) * acp ** 2 / pcp / 1e6, rel=1e-3)


def test_vs_status_reaches_the_caller():
    """Regression pin on P3. shear_status was computed and discarded,
    overall_check being returned under that key, so
    "UNSAFE - Vs exceeds limit" could never reach a caller.

    The two are the same inequality, so they always agree -- this pins
    that the string is now reachable, not that a check was restored.
    """
    ok = ex_7_2()
    assert ok["vs_status"] in ("SAFE", "Provide Min. Rft")
    assert ok["shear_status"] == "SAFE"

    over = ex_7_2(vu=5000.0)
    assert over["vs_status"] == "UNSAFE - Vs exceeds limit"
    assert over["shear_status"] == "UNSAFE"


# ------------------------------------------------------------- unit sanity

def test_threshold_is_in_kilonewton_metres_not_newton_metres():
    """THE test for the class of bug this whole module exists because of.

    A 355.6 x 609.6 beam in f'c = 27.6 MPa has phi*Tth of order 8 kN.m.
    The deleted implementation returned 7957.84 for the same section --
    the same number 1000x too large, in N.m, compared against a kN.m
    demand. Anything outside 1..100 here is a unit error.
    """
    r = ex_7_2()
    assert 1.0 < r["tth"] < 100.0, "phi*Tth is off by an order of magnitude"
    assert 1.0 < r["tcr"] < 400.0
    assert r["tcr"] == pytest.approx(r["tth"] * 0.33 / 0.083, rel=1e-4)


def test_at_s_is_in_square_millimetres_per_millimetre():
    """At/s of order 0.5 mm2/mm. The deleted implementation returned
    0.000558 for the same section -- 1000x low."""
    r = ex_7_2()
    assert 0.05 < r["at_s"] < 5.0, "At/s is off by an order of magnitude"


def test_al_is_in_square_millimetres():
    """Al of order 800 mm2 for this section, not 0.8 and not 800 000."""
    r = ex_7_2()
    assert 100.0 < r["al"] < 5000.0


def test_torsion_scales_as_the_cube_of_a_uniform_size_change():
    """Tth ~ Acp^2/pcp scales as L^3. Doubling every dimension must
    multiply the threshold by 8, which no stray constant survives."""
    small = ex_7_2()
    big = ex_7_2(bw=28 * IN, h=48 * IN, d=43 * IN, cc=3.0 * IN,
                 c=(48 - 43) * IN)
    assert rel(big["tth"] / small["tth"], 8.0) < 0.01


# ------------------------------------------------- no silent acceptance

def test_a_beam_over_the_threshold_is_never_told_to_neglect_torsion():
    """The headline defect, asserted directly. Tu = 37.96 kN.m against a
    threshold of 7.96 kN.m is 4.8x over. The deleted module answered
    "MAY NEGLECT" and "status": "PASS" for exactly this section."""
    r = ex_7_2()
    assert r["torsion_action"] == "Design for Torsion"
    assert 28.0 * KIPFT > 4.5 * r["tth"], "the demand really is far over"
    assert r["at_s"] > 0.0, "torsion stirrups must be required"
    assert r["al"] > 0.0, "longitudinal torsion steel must be required"


def test_a_crushed_section_reports_unsafe():
    """Section 22.7.7.1 cross-sectional limit, printed 466. A section
    that cannot carry the combined shear and torsion stress must say so,
    not return a reinforcement schedule that reads as adequate."""
    r = ex_7_2(tu=1000.0, vu=2000.0)
    assert r["dim_check"] == "UNSAFE"


def test_negligible_torsion_returns_zero_not_a_spurious_area():
    """Below the threshold, At/s and Al are zero -- an explicit zero
    demand, not a small nonzero number from an unguarded formula."""
    r = ex_7_2(tu=0.5)
    assert r["torsion_action"] == "Neglect Torsion"
    assert r["at_s"] == 0.0
    assert r["al"] == 0.0
