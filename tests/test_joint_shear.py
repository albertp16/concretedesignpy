# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
Joint shear tests -- concretedesignpy.calculators.joint_shear
==============================================================

Governing edition: NSCP 2015 Section 418.8.4 (= ACI 318M-14).
Benchmark: Wight & MacGregor 7th ed. Printed pages are the printed folio.
"""

import math

import pytest

from concretedesignpy.calculators.joint_shear import (
    joint_shear_check,
    NSCP_2015_TABLE_418_8_4_1,
    LEGAL_LAMBDA,
    PHI_SMF_JOINT,
)

IN = 25.4
KSI = 6.894757
KIP = 4.44822


def rel(actual, expected):
    return abs(actual / expected - 1.0)


def ex_19_3(**over):
    """W&M Example 19-3 (printed 1076-1078) converted to SI.

    Interior joint, 24 in. square column, 24 in. wide beams,
    f'c = 4 ksi, fy = 60 ksi, column shear Ve = 81.8 kips.
    As,top = 4.36 in2, As,bot = 2.24 in2.
    Beams flush with the column, so x = 0.
    """
    kwargs = dict(
        ve=81.8 * KIP,
        as1=4.36 * IN * IN, n_bars1=1,
        as2=2.24 * IN * IN, n_bars2=1,
        fy=60 * KSI, fc=4 * KSI,
        beam_width=24 * IN, joint_depth=24 * IN, column_width=24 * IN,
        perpendicular_dist=0, joint_config=1,
    )
    kwargs.update(over)
    return joint_shear_check(**kwargs)


# ---------------------------------------------------------------- textbook

def test_wm_example_19_3_demand():
    """W&M printed 1078: T1 = 327 kips, C2 = T2 = 168 kips,
    Vj = T1 + C2 - Vcol = 413 kips. The 1.25 fy overstrength is
    Section 18.8.2.1."""
    r = ex_19_3()
    assert rel(r["t1"], 327 * KIP) < 0.01
    assert rel(r["t2"], 168 * KIP) < 0.01
    assert rel(r["v_joint"], 413 * KIP) < 0.01


def test_wm_example_19_3_effective_joint_area():
    """W&M printed 1078: Aj = 576 in2 = 371 612 mm2, the full column
    section for a beam flush with a square column."""
    r = ex_19_3()
    assert rel(r["aj"], 576 * IN * IN) < 0.01
    assert rel(r["joint_width"], 24 * IN) < 0.01


def test_wm_example_19_3_capacity():
    """W&M printed 1078: phi*Vn = 619 kips.

    TOLERANCE IS 3%, NOT 1%. W&M works in psi with gamma = 20; NSCP 2015
    Table 418.8.4.1 and ACI 318M print 1.7, and 20/12 = 1.667. The 2.4%
    gap is the SI-versus-inch-pound rounding of the same coefficient, not
    an error -- see W&M Eq. (17-22b), printed 967.
    """
    r = ex_19_3()
    assert rel(r["phi_vn"], 619 * KIP) < 0.03
    assert r["phi_vn"] > r["v_joint"]
    assert r["status"] == "OK"


# -------------------------------------------------------------- regression

@pytest.mark.parametrize("beam_b,aj_wrong", [
    (450.0, 270000.0), (500.0, 300000.0),
    (600.0, 360000.0), (700.0, 420000.0),
])
def test_column_width_caps_aj_for_a_beam_wider_than_the_column(beam_b, aj_wrong):
    """Regression pin on P1b/F-10.

    Section 15.5.2.2, printed 232: "Effective joint width shall be the
    overall width of the column where the beam is wider than the column."
    R15.5.2.2, printed 231-232: "In no case is Aj greater than the column
    cross-sectional area."

    A spandrel on a 400 x 600 column. Aj must be 240 000 mm2 in every
    case; the second column is what the function used to return.
    """
    r = ex_19_3(beam_width=beam_b, joint_depth=600.0, column_width=400.0,
                fy=420, fc=28, ve=0,
                as1=2000, n_bars1=1, as2=1000, n_bars2=1)
    assert r["joint_width"] == 400.0
    assert r["aj"] == 240000.0
    assert r["aj"] < aj_wrong


def test_beam_narrower_than_the_column_is_unaffected():
    """The b + 2x limb still governs when it should. Section 15.5.2.2(b)
    reads "twice the perpendicular distance from longitudinal axis of
    beam to nearest side face of the column" = 2(b/2 + x) = b + 2x, with
    x measured beam face to column face as this function defines it.

    300 wide beam, x = 50, on a 600 wide column:
      column_width  600
      b + h         300 + 600 = 900
      b + 2x        300 + 100 = 400   <- governs
    """
    r = ex_19_3(beam_width=300.0, joint_depth=600.0, column_width=600.0,
                perpendicular_dist=50.0, fy=420, fc=28, ve=0,
                as1=1000, n_bars1=1, as2=500, n_bars2=1)
    assert r["joint_width"] == 400.0
    assert r["aj"] == 240000.0


def test_column_width_is_required_and_keyword_only():
    """Regression pin on P1b. Omitting it must fail loudly rather than
    silently reproducing the old, uncapped answer."""
    with pytest.raises(TypeError, match="column_width"):
        joint_shear_check(
            ve=0, as1=1000, n_bars1=1, as2=500, n_bars2=1,
            fy=420, fc=28, beam_width=500, joint_depth=600,
        )


def test_gamma_table_is_the_nscp_three_row_one():
    """Regression pin on P2b/F-11. Option A: the package targets
    NSCP 2015, so the three-row Table 418.8.4.1 is correct as written.
    ACI 318-25M Table 18.8.4.3, printed 342, is an eight-row successor
    (1.7/1.3/1.3/1.0/1.3/1.0/1.0/0.7) and moving to it is an edition
    change, not a fix."""
    assert NSCP_2015_TABLE_418_8_4_1 == {1: 1.7, 2: 1.2, 3: 1.0}
    for cfg, gamma in NSCP_2015_TABLE_418_8_4_1.items():
        r = ex_19_3(joint_config=cfg)
        assert r["gamma"] == gamma
        expected = 0.85 * gamma * math.sqrt(4 * KSI) * r["aj"] / 1000.0
        assert r["phi_vn"] == pytest.approx(expected, rel=1e-3)


# ------------------------------------------------------------- unit sanity

def test_forces_are_in_kilonewtons():
    """T1 = 1.25 * 4.36 in2 * 60 ksi = 327 kips = 1455 kN. Anything
    outside 100..10 000 is a unit error."""
    r = ex_19_3()
    assert 100.0 < r["t1"] < 10000.0
    assert 100.0 < r["phi_vn"] < 100000.0


def test_aj_is_in_square_millimetres():
    r = ex_19_3()
    assert 1e4 < r["aj"] < 1e7


def test_capacity_scales_with_the_square_root_of_fc():
    """Vn ~ sqrt(f'c). Quadrupling f'c must double phi*Vn."""
    lo = ex_19_3(fc=25.0)
    hi = ex_19_3(fc=100.0)
    assert rel(hi["phi_vn"] / lo["phi_vn"], 2.0) < 0.001


def test_capacity_is_linear_in_aj():
    a = ex_19_3(joint_depth=600.0, column_width=400.0, beam_width=400.0)
    b = ex_19_3(joint_depth=600.0, column_width=800.0, beam_width=800.0)
    assert rel(b["phi_vn"] / a["phi_vn"], 2.0) < 0.001


# ------------------------------------------------- no silent acceptance

def test_an_inadequate_joint_says_so():
    """A joint whose demand exceeds phi*Vn must return
    REINFORCEMENT REQUIRED, not OK."""
    r = ex_19_3(beam_width=300.0, joint_depth=300.0, column_width=300.0,
                fc=21.0)
    assert r["phi_vn"] < r["v_joint"]
    assert r["status"] == "REINFORCEMENT REQUIRED"


@pytest.mark.parametrize("bad_lamda", [0.0, 0.5, 0.85, 0.9, 1.2])
def test_illegal_lamda_raises(bad_lamda):
    """Regression pin on P2b. lamda was a free float. ACI 318-25M
    Table 15.5.2.1 footnote [1], printed 231: lambda shall be 0.75 for
    lightweight and 1.0 for normalweight. Two legal values."""
    with pytest.raises(ValueError, match="lamda"):
        ex_19_3(lamda=bad_lamda)


@pytest.mark.parametrize("good_lamda", LEGAL_LAMBDA)
def test_legal_lamda_accepted(good_lamda):
    r = ex_19_3(lamda=good_lamda)
    assert r["lamda"] == good_lamda


@pytest.mark.parametrize("bad_phi", [0.75, 0.9, 1.0])
def test_overriding_phi_raises(bad_phi):
    """Regression pin on P2b. ACI 318-25M Section 21.2.4.4, printed 435:
    "For beam-column joints of special moment frames ... phi for shear
    shall be 0.85." No caller discretion."""
    with pytest.raises(ValueError, match="phi"):
        ex_19_3(phi=bad_phi)


def test_phi_default_is_the_code_value():
    assert PHI_SMF_JOINT == 0.85
    assert ex_19_3()["phi"] == 0.85


@pytest.mark.parametrize("bad_config", [0, 4, 99, -1])
def test_unknown_joint_config_raises_instead_of_scoring_1_0(bad_config):
    """Regression pin on P2b. The lookup was dict.get(cfg, 1.0), so a
    typo silently scored the 'other joints' row -- the most favourable of
    the three that is not 1.7."""
    with pytest.raises(ValueError, match="joint_config"):
        ex_19_3(joint_config=bad_config)


def test_nonpositive_column_width_raises():
    with pytest.raises(ValueError, match="column_width"):
        ex_19_3(column_width=0)
