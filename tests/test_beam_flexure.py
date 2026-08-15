# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam flexure tests -- concretedesignpy.calculators.beam_moment
==============================================================

Four kinds of test, per the rectification brief:

1. textbook pins   -- one per worked example, printed page in the docstring
2. regression pins -- the exact numbers in the P3 acceptance line
3. unit sanity     -- order of magnitude, so a 1000x error fails loudly
4. no silent acceptance -- an infeasible section must raise, never return

Governing edition: NSCP 2015 (= ACI 318M-14).
Benchmark: Wight & MacGregor, *Reinforced Concrete Mechanics and Design*,
7th ed. Printed pages are the printed folio (PDF - 1).
"""

import math

import pytest

from concretedesignpy.calculators.beam_moment import calculate_beam_moment

IN = 25.4                # mm per inch
KSI = 6.894757           # MPa per ksi
KIPFT = 1.35582          # kN.m per kip-ft


def rel(actual, expected):
    return abs(actual / expected - 1.0)


# ---------------------------------------------------------------- textbook

def test_wm_example_4_1m_singly_reinforced_si():
    """W&M Example 4-1M, printed 148.

    250 x 565, d = 500, As = 3 x 510 = 1530 mm2 (3-No.25M),
    f'c = 20 MPa, fy = 420 MPa.
    Book: a = 151 mm, c = 178 mm, eps_s = 0.00543, Mn = 273 kN.m.
    """
    r = calculate_beam_moment(
        rebar_list=[{"d": 500, "diam": 25.4692, "num": 3}],
        fc=20, fy=420, b=250, h=565,
    )
    assert rel(r["a"], 151.0) < 0.01
    assert rel(r["neutral_axis"], 178.0) < 0.01
    assert rel(r["epsilon_t"], 0.00543) < 0.01
    assert rel(r["mn"], 273.0) < 0.01
    assert r["classification"] == "tension-controlled"
    assert r["phi"] == pytest.approx(0.90)


def test_wm_example_4_4_doubly_reinforced():
    """W&M Example 4-4, printed 170-171.

    12 x 20 in., d = 16.8 in., d' = 2.5 in.,
    As = 6-No.9 = 6.00 in2, A's = 3-No.9 = 3.00 in2,
    f'c = 4 ksi, fy = 60 ksi, Es = 29 000 ksi.
    Book: c = 6.2 in., Mn = 428 k-ft.

    This is the example that fixes the compression-steel force law:
    step 5 (printed 170) writes Cs = A's (f's - 0.85 f'c), deducting the
    concrete the bar displaces inside the Whitney block.
    """
    a_no9 = 1.00 * IN * IN
    d_bar = math.sqrt(4 * a_no9 / math.pi)
    r = calculate_beam_moment(
        rebar_list=[{"d": 2.5 * IN, "diam": d_bar, "num": 3},
                    {"d": 16.8 * IN, "diam": d_bar, "num": 6}],
        fc=4 * KSI, fy=60 * KSI, b=12 * IN, h=20 * IN, es=29000 * KSI,
    )
    assert rel(r["mn"], 428.0 * KIPFT) < 0.01
    assert rel(r["neutral_axis"], 6.2 * IN) < 0.01


# -------------------------------------------------------------- regression

def test_displaced_concrete_deduction_is_applied_inside_the_block():
    """Regression pin on P3/F-1.

    A compression bar inside depth a carries A's (f's - 0.85 f'c), not
    A's f's. Basis: W&M Ex 4-4 step 5, printed 170.

    W&M Ex 4-4: before the fix c came out 5.99 in. (-3.35% against the
    book's 6.2 in.); after, 6.184 in. (-0.26%). Mn barely moves -- the
    error was in c, not Mn -- so c is what this pins.
    """
    a_no9 = 1.00 * IN * IN
    d_bar = math.sqrt(4 * a_no9 / math.pi)
    r = calculate_beam_moment(
        rebar_list=[{"d": 2.5 * IN, "diam": d_bar, "num": 3},
                    {"d": 16.8 * IN, "diam": d_bar, "num": 6}],
        fc=4 * KSI, fy=60 * KSI, b=12 * IN, h=20 * IN, es=29000 * KSI,
    )
    assert r["neutral_axis"] == pytest.approx(157.07, abs=0.5)

    top = min(r["rebar_forces"], key=lambda f: f["d"])
    assert top["d"] < r["a"], "this pin only means anything inside the block"
    # gross (undeducted) compressive stress at the top bar, negative;
    # taken from the bar's own reported strain, not from the rounded c
    fsp_gross = max(-60 * KSI, top["strain"] * 29000 * KSI)
    # reported stress is the gross stress relieved by 0.85 f'c
    assert top["stress"] == pytest.approx(fsp_gross + 0.85 * 4 * KSI, rel=1e-9)
    assert top["stress"] > fsp_gross, "deduction must reduce the bar force"


def test_deduction_not_applied_to_a_bar_below_the_stress_block():
    """A bar in compression but below depth a sits in cracked concrete
    and displaces nothing the Whitney block already counted, so no
    deduction applies to it."""
    r = calculate_beam_moment(
        rebar_list=[{"d": 60.0, "diam": 12.0, "num": 2},
                    {"d": 540.0, "diam": 25.0, "num": 4}],
        fc=28, fy=420, b=300, h=600,
    )
    a = r["a"]
    for f in r["rebar_forces"]:
        if f["stress"] < 0 and f["d"] >= a:
            gross = 0.003 * (r["neutral_axis"] - f["d"]) / r["neutral_axis"] * 200000
            assert f["stress"] == pytest.approx(max(-420.0, gross), rel=1e-6)


# ------------------------------------------------------------- unit sanity

def test_order_of_magnitude_singly_reinforced():
    """Hand check to one significant figure. A 1000x unit error fails.

    300 x 600, d = 535, As = 4-No.25 = 1963 mm2, f'c = 28, fy = 420.
      a  = As fy / (0.85 f'c b) = 1963*420/(0.85*28*300) = 115 mm
      Mn = As fy (d - a/2)      = 1963*420*(535-57.7)   = 393 kN.m
    """
    r = calculate_beam_moment(
        rebar_list=[{"d": 535, "diam": 25.0, "num": 4}],
        fc=28, fy=420, b=300, h=600,
    )
    assert 350.0 < r["mn"] < 440.0, "Mn is out by an order of magnitude"
    assert 100.0 < r["a"] < 130.0
    assert 0.0 < r["neutral_axis"] < 600.0
    assert 0.65 <= r["phi"] <= 0.90


@pytest.mark.parametrize("n_bars,expected", [(2, 209.0), (4, 393.0), (6, 555.0)])
def test_mn_tracks_the_steel_it_is_given(n_bars, expected):
    """Mn approximately As fy (d - a/2) at fixed geometry. Hand values:

        300 x 600, d = 535, No.25 bars (491 mm2), f'c = 28, fy = 420
        n=2: As=982   a= 57.7  Mn = 982*420*(535-28.9)  = 209 kN.m
        n=4: As=1963  a=115.5  Mn = 1963*420*(535-57.7) = 393 kN.m
        n=6: As=2945  a=173.2  Mn = 2945*420*(535-86.6) = 555 kN.m

    Tolerance is loose (15%) on purpose: this is an order-of-magnitude
    guard against a stray unit factor, not a precision check. The
    textbook pins above are the precision checks.
    """
    r = calculate_beam_moment(
        rebar_list=[{"d": 535, "diam": 25.0, "num": n_bars}],
        fc=28, fy=420, b=300, h=600)
    assert rel(r["mn"], expected) < 0.15


# ------------------------------------------------- no silent acceptance

def test_no_equilibrium_raises_rather_than_reporting_a_capacity():
    """When the compression side outruns every available tension force
    there is no equilibrium, and the solver must raise -- not return a
    number. A 10 mm top plate on a 300 x 600 with only 2-No.10 in
    tension is such a section.

    This module already had the right idiom (unlike three sibling repos
    in the same audit family); the test locks it in.
    """
    with pytest.raises(ValueError, match="equilibrium"):
        calculate_beam_moment(
            rebar_list=[{"d": 540.0, "diam": 10.0, "num": 2}],
            fc=21, fy=420, b=300, h=600, plate_top_thickness=10.0,
        )


def test_a_hopeless_section_never_reports_a_passing_verdict():
    """Whatever a near-degenerate section returns, it must not be a
    capacity that reads as adequate. Ten No.32 bars 30 mm from the top
    of a 250 x 400 converge on a token Mn; assert it stays token rather
    than silently reporting a full-depth capacity."""
    r = calculate_beam_moment(
        rebar_list=[{"d": 30.0, "diam": 32.0, "num": 10}],
        fc=21, fy=420, b=250, h=400)
    assert r["mn"] < 20.0
    assert r["neutral_axis"] < 40.0


def test_empty_rebar_list_raises():
    with pytest.raises(ValueError, match="rebar_list"):
        calculate_beam_moment(rebar_list=[], fc=28, fy=420, b=300, h=600)


@pytest.mark.parametrize("bad", [
    {"fc": 0}, {"fc": -28}, {"fy": 0}, {"b": 0}, {"h": -600},
])
def test_nonpositive_inputs_raise(bad):
    kwargs = dict(rebar_list=[{"d": 535, "diam": 25.0, "num": 4}],
                  fc=28, fy=420, b=300, h=600)
    kwargs.update(bad)
    with pytest.raises(ValueError):
        calculate_beam_moment(**kwargs)


def test_module_does_not_claim_to_check_as_min():
    """calculate_beam_moment computes capacity; it does not certify a
    section. Section 9.6.1.2 As,min is NOT checked, and the docstring
    must keep saying so."""
    doc = calculate_beam_moment.__module__
    import importlib
    mod = importlib.import_module(doc)
    assert "As,min" in mod.__doc__
    assert "9.6.1.2" in mod.__doc__
    # and the header must not have drifted back to claiming 318-19
    assert "ACI 318-19" not in mod.__doc__.split("Governing edition")[0]
