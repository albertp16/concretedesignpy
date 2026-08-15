# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Beam shear tests -- concretedesignpy.calculators.beam_shear
============================================================

Governing edition: NSCP 2015 Section 422.5 (= ACI 318M-14).
Benchmark: Wight & MacGregor 7th ed. Printed pages are the printed folio.
"""

import math

import pytest

from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength,
    compute_steel_shear_strength,
    compute_shear_spacing,
    shear_design,
)


def rel(actual, expected):
    return abs(actual / expected - 1.0)


# ---------------------------------------------------------------- textbook

def test_wm_example_6_1m_vc():
    """W&M Example 6-1M, printed 297-300.

    T-beam bw = 300, d = 610, f'c = 25. Book Vc = 153 kN, from
    Vc = sqrt(f'c) bw d / 6.
    """
    vc = compute_concrete_shear_strength(fc=25, b=300, d=610)
    assert rel(vc["vc_kn"], 153.0) < 0.01


def test_wm_example_6_1m_vs():
    """W&M Example 6-1M, printed 299. No.10M stirrups (Av = 142 mm2) at
    s = 300 give Vs = 86.6 kN, fyt = 300."""
    vs = compute_steel_shear_strength(av=142, fyt=300, d=610, s=300)
    assert rel(vs["vs_kn"], 86.6) < 0.01


def test_wm_example_6_1m_required_spacing():
    """W&M Example 6-1M, printed 300. No.13M stirrups (Av = 258 mm2) at
    Vu/phi = 373 kN require s = 215 mm."""
    r = compute_shear_spacing(fc=25, b=300, d=610, fyt=300,
                              vu_required=373e3 * 0.75, phi=0.75, av=258)
    assert r["status"] == "OK"
    assert rel(r["spacing"], 215.0) < 0.01


def test_wm_example_6_1m_through_shear_design():
    """The same case through shear_design(). Book Vs,max = 610 kN."""
    sd = shear_design(fc=25, fyv=300, phi=0.75, bw=300, h=675, cc=40, c=65,
                      d=610, vu=373 * 0.75, nu=0, s_chosen=200, n_legs=2,
                      db_stirrup=12.7)
    assert rel(sd["vc"], 153.0) < 0.01
    assert rel(sd["vs_max"], 610.0) < 0.01
    assert sd["shear_status"] == "SAFE"


# -------------------------------------------------------------- regression

def test_avmin_governs_on_the_strength_branch():
    """Regression pin on P2a/F-3.

    Av,min used to be applied only when Vu <= phi*Vc. On a wide, lightly
    loaded web with Vu just above phi*Vc, the strength requirement alone
    wants an absurd spacing and Av,min is what must govern.

    600 x 800 web, f'c 28, fyt 420, Av 157 mm2:
      strength alone  s ~ 15 750 mm
      s_max = d/2     = 400 mm
      Av,min          = max(0.062*sqrt(28)*600/420, 0.35*600/420)
                      = max(0.4688, 0.5) = 0.5 mm2/mm  ->  s = 314 mm
    Basis: ACI 318-25M Table 9.6.3.4, printed 151.
    """
    vc = math.sqrt(28) / 6 * 600 * 800
    r = compute_shear_spacing(fc=28, b=600, d=800, fyt=420,
                              vu_required=0.75 * vc + 1000, phi=0.75, av=157)
    assert r["status"] == "OK"
    assert r["av_min_per_s"] == pytest.approx(0.5, abs=1e-3)
    assert r["s_avmin"] == pytest.approx(314.0, abs=1.0)
    assert r["spacing"] == pytest.approx(314.0, abs=1.0)
    assert r["spacing"] < r["smax"], "Av,min, not s_max, must govern here"


def test_avmin_constant_is_0_062_everywhere():
    """Regression pin on P3/F-6. compute_shear_spacing() used 0.062 and
    the two design functions used 1/16 = 0.0625 -- two constants 0.8%
    apart for Table 9.6.3.4(a), printed 151. Now one."""
    # branch (a) governs branch (b) only above f'c = (0.35/0.062)^2
    # = 31.9 MPa, so pick f'c = 50 to put the constant under test.
    fc, bw, fyv = 50.0, 300.0, 420.0
    expected = 0.062 * math.sqrt(fc) * bw / fyv
    assert expected > 0.35 * bw / fyv, "test must exercise branch (a)"

    sd = shear_design(fc=fc, fyv=fyv, phi=0.75, bw=bw, h=675, cc=40, c=65,
                      d=610, vu=100, nu=0, s_chosen=200, n_legs=2,
                      db_stirrup=12.7)
    sp = compute_shear_spacing(fc=fc, b=bw, d=610, fyt=fyv,
                               vu_required=100e3, phi=0.75, av=258)
    assert sd["av_min"] == pytest.approx(expected, rel=1e-3)
    assert sp["av_min_per_s"] == pytest.approx(expected, rel=1e-3)
    assert sd["av_min"] == pytest.approx(sp["av_min_per_s"], rel=1e-6)


def test_axial_compression_vc_has_no_second_branch():
    """Regression pin on P3/F-5.

    The deleted cap 0.3 sqrt(f'c) bw d (1 + 0.3 Nu/Ag) was not supported
    by reference/ and never governed: its ratio to the governing branch
    is 1.8 (1 + 0.3q)/(1 + q/14), q = Nu/Ag, minimum 1.8 at Nu = 0.
    Vc must now equal W&M Eq. (6-13aM), printed 282, exactly.
    """
    fc, bw, h, d, nu_kn = 28.0, 350.0, 600.0, 530.0, 500.0
    ag = bw * h
    expected = (1 + nu_kn * 1000 / (14 * ag)) * math.sqrt(fc) / 6 * bw * d / 1000
    sd = shear_design(fc=fc, fyv=400, phi=0.75, bw=bw, h=h, cc=40, c=70, d=d,
                      vu=300, nu=nu_kn, s_chosen=150, n_legs=2, db_stirrup=10)
    assert sd["vc"] == pytest.approx(expected, rel=1e-3)
    assert sd["vc_note"] == "with axial compression"


def test_the_deleted_cap_provably_never_governed():
    """The claim behind F-5's downgrade, asserted rather than trusted."""
    for q in (0.0, 0.5, 1.0, 5.0, 20.0, 100.0):       # q = Nu/Ag, MPa
        ratio = 1.8 * (1 + 0.3 * q) / (1 + q / 14.0)
        assert ratio >= 1.8


# ------------------------------------------------------------- unit sanity

def test_vc_order_of_magnitude():
    """300 x 600, d = 535, f'c = 28: Vc = sqrt(28)/6 * 300 * 535
    = 141 600 N = 142 kN. A 1000x error fails loudly."""
    vc = compute_concrete_shear_strength(fc=28, b=300, d=535)
    assert 100.0 < vc["vc_kn"] < 200.0
    assert vc["vc"] == pytest.approx(vc["vc_kn"] * 1000, rel=1e-3)


def test_vs_order_of_magnitude():
    """Av = 157 mm2, fyt = 420, d = 535, s = 150:
    Vs = 157*420*535/150 = 235 000 N = 235 kN."""
    vs = compute_steel_shear_strength(av=157, fyt=420, d=535, s=150)
    assert 200.0 < vs["vs_kn"] < 270.0


def test_spacing_stays_in_millimetres():
    """A returned spacing must be a plausible stirrup pitch, not a
    number that has slipped a factor of 1000 in either direction."""
    r = compute_shear_spacing(fc=28, b=300, d=535, fyt=420,
                              vu_required=200e3, phi=0.75, av=157)
    assert r["status"] == "OK"
    assert 25.0 < r["spacing"] <= r["smax"] <= 600.0


# ------------------------------------------------- no silent acceptance

def test_overstressed_section_returns_no_spacing():
    """Regression pin on P2a/F-4. Above the Section 22.5.1.2 limit
    (printed 442) there is no stirrup spacing that works, so the
    function must not return one."""
    r = compute_shear_spacing(fc=28, b=300, d=500, fyt=420,
                              vu_required=0.75 * 3000e3, phi=0.75, av=157)
    assert r["spacing"] is None
    assert r["status"].startswith("UNSAFE")
    assert "22.5.1.2" in r["status"]
    assert r["vs_required"] > r["vs_max"]


def test_shear_design_reports_unsafe_rather_than_a_number():
    sd = shear_design(fc=28, fyv=420, phi=0.75, bw=300, h=600, cc=40, c=65,
                      d=535, vu=3000.0, nu=0, s_chosen=100, n_legs=2,
                      db_stirrup=12)
    assert sd["shear_status"] == "UNSAFE"
    assert sd["vu_max"] < 3000.0


def test_zero_spacing_raises():
    with pytest.raises(ValueError, match="spacing"):
        compute_steel_shear_strength(av=157, fyt=420, d=535, s=0)


@pytest.mark.parametrize("bad", [{"fc": 0}, {"b": 0}, {"d": -535}])
def test_nonpositive_geometry_raises(bad):
    kwargs = dict(fc=28, b=300, d=535)
    kwargs.update(bad)
    with pytest.raises(ValueError):
        compute_concrete_shear_strength(**kwargs)


def test_unknown_vc_type_raises():
    with pytest.raises(ValueError, match="vc_type"):
        compute_concrete_shear_strength(fc=28, b=300, d=535, vc_type="wishful")
