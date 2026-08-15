"""Parity tests for the VENDORED jacketing engine.

``concretedesignpy/calculators/column_jacket.py`` is a byte-for-byte copy of the
engine maintained in the APEC RC Column Jacketing service, which in turn vendors
it from the workspace draft.  Vendoring breaks single-source-of-truth by
construction, so these tests exist to make any drift loud.

Two kinds of test here, and the distinction matters when you add one:

*   **Independent recomputation** — the reference value is derived from the code
    equation inside the test, not pasted from a previous run.  These catch the
    engine changing behaviour.  They are the real protection; prefer them.
*   **Regression pins** — a recorded value with a tight tolerance, used only
    where independent recomputation would mean reimplementing a layer
    integration.  These catch accidental change but would happily pin a bug, so
    they are labelled.  If one fails, work out which behaviour changed before
    updating the number.

Reference case: the worked example in TN-RET-001 Project Applications —
existing 400x400, f'c 21 MPa, 8-D20 Grade 275; jacket 100 mm to 600x600,
f'c 28 MPa, 12-D20 Grade 415, D12 ties at 100, plain perimeter hoop.
"""

import math

import pytest

from concretedesignpy.calculators.column_jacket import (
    Bar,
    Concrete,
    JacketedColumn,
    Steel,
    TieSet,
    capacity_at_P,
    rect_perimeter_bars,
)
from concretedesignpy.calculators.column_jacket_design import ENGINE_SHA256

BE = HE = 400.0
T = 100.0
FC_E, FC_J = 21.0, 28.0
FY_E, FY_J = 275.0, 415.0
TIE_DIA, TIE_S, TIE_FY = 12.0, 100.0, 275.0
COVER_TIE = 40.0
DB = 20.0


@pytest.fixture(scope="module")
def col():
    bars_e = rect_perimeter_bars(BE, HE, 50.0, 3, 3, DB, Steel(FY_E),
                                 "existing", T)
    bars_j = rect_perimeter_bars(BE + 2 * T, HE + 2 * T,
                                 COVER_TIE + TIE_DIA + DB / 2, 4, 4, DB,
                                 Steel(FY_J), "jacket", 0.0)
    leg = math.pi * TIE_DIA ** 2 / 4 * 2
    ties = TieSet(Av_x=leg, Av_z=leg, s=TIE_S, fyt=TIE_FY, db=TIE_DIA)
    return JacketedColumn(BE, HE, T, Concrete(FC_E), Concrete(FC_J),
                          bars_e, bars_j, ties, cover_j=COVER_TIE,
                          n_sup_b=2, n_sup_h=2, db_long=DB)


# ------------------------------------------------------------------ geometry


def test_bar_counts_follow_the_perimeter_rule(col):
    """n per face on a perimeter, corners counted once, gives 4(n-1)."""
    assert len(col.bars_e) == 4 * (3 - 1) == 8
    assert len(col.bars_j) == 4 * (4 - 1) == 12


def test_jacket_bars_are_uniformly_spaced(col):
    """Equally spaced perimeter bars must have equal gaps -- an unequal set
    means the layout generator has drifted."""
    ys = sorted({round(b.y, 4) for b in col.bars_j})
    gaps = [b - a for a, b in zip(ys, ys[1:])]
    assert all(g == pytest.approx(gaps[0], abs=1e-3) for g in gaps), \
        "jacket bar rows not uniformly spaced: {}".format(ys)


# --------------------------------------------------- independent recomputation


def test_Po_matches_independent_recomputation(col):
    """ACI 318-19 Eq. (22.4.2.2), applied per concrete region.

    Derived here from the equation, not copied from a run.
    """
    Ag_e = BE * HE
    Ag_j = (BE + 2 * T) * (HE + 2 * T) - Ag_e
    expected = (0.85 * FC_E * (Ag_e - col.Ast_e) + FY_E * col.Ast_e
                + 0.85 * FC_J * (Ag_j - col.Ast_j) + FY_J * col.Ast_j)
    assert col.Po() == pytest.approx(expected, rel=1e-12)
    # and the published hand-check value, to 0.01%
    assert col.Po() / 1e3 == pytest.approx(9737.1, rel=1e-4)


def test_phiPn_max_uses_0_80_and_0_65_for_ties(col):
    """Table 22.4.2.1 (0.80 P_o for ties) and Table 21.2.2 (phi = 0.65)."""
    assert col.alpha_max == 0.80
    assert col.phi_comp == 0.65
    assert col.phiPn_max() == pytest.approx(0.65 * 0.80 * col.Po(), rel=1e-12)


def test_beta1_follows_table_22_2_2_4_3():
    assert Concrete(21.0).beta1 == pytest.approx(0.85)
    assert Concrete(28.0).beta1 == pytest.approx(0.85)
    assert Concrete(45.0).beta1 == pytest.approx(0.85 - 0.05 * (45 - 28) / 7)
    assert Concrete(90.0).beta1 == pytest.approx(0.65)          # floor


def test_shear_matches_independent_recomputation(col):
    """ACI 318-19 Table 22.5.5.1(a) SI, with the 22.5.5.1.2 cap on the axial
    term."""
    Nu = 3200e3
    sh = col.shear_capacity(Nu=Nu)
    fc, lam, bw = min(FC_E, FC_J), 1.0, BE + 2 * T
    d = max(b.y for b in col.bars)
    Ag = bw * (HE + 2 * T)

    axial = min(Nu / (6 * Ag), 0.05 * fc)                       # 22.5.5.1.2
    assert axial == pytest.approx(0.05 * fc), "the Nu/(6Ag) cap should bind here"
    Vc = (0.17 * lam * math.sqrt(fc) + axial) * bw * d
    assert Vc <= 0.42 * lam * math.sqrt(fc) * bw * d, "the 0.42 cap should NOT bind"
    Vs = (math.pi * TIE_DIA ** 2 / 4 * 2) * TIE_FY * d / TIE_S

    assert sh["d"] == pytest.approx(d)
    assert sh["Vc"] == pytest.approx(Vc, rel=1e-12)
    assert sh["Vs"] == pytest.approx(Vs, rel=1e-12)
    assert sh["phiVn"] == pytest.approx(0.75 * (Vc + Vs), rel=1e-12)


def test_confinement_matches_independent_recomputation(col):
    """Mander via fib B14 Eq. (6-35), (6-5), (6-36)."""
    cf = col.confined("existing")
    bc = hc = (BE + 2 * T) - 2 * (COVER_TIE + TIE_DIA / 2)
    leg = math.pi * TIE_DIA ** 2 / 4 * 2
    rho_s = (leg * hc + leg * bc) / (TIE_S * bc * hc)
    assert col.rho_s_ties == pytest.approx(rho_s, rel=1e-12)

    fl = 0.5 * col.mander_ke() * rho_s * TIE_FY
    r = fl / FC_E
    fcc = FC_E * (2.254 * math.sqrt(1 + 7.94 * r) - 2 * r - 1.254)
    eps_cu = 0.004 + 1.4 * rho_s * TIE_FY * 0.09 / fcc

    assert cf.fl == pytest.approx(fl, rel=1e-12)
    assert cf.fcc == pytest.approx(fcc, rel=1e-12)
    assert cf.eps_cu == pytest.approx(eps_cu, rel=1e-12)


def test_mander_ke_matches_independent_recomputation(col):
    """k_e from the tie and bar geometry, not assumed."""
    bc = hc = (BE + 2 * T) - 2 * (COVER_TIE + TIE_DIA / 2)
    d_bar = COVER_TIE + TIE_DIA + DB / 2
    span = (BE + 2 * T) - 2 * d_bar
    w = span / (2 - 1) - DB                       # 2 restrained bars per face
    sum_w2 = 2 * 1 * w ** 2 + 2 * 1 * w ** 2      # four faces, one gap each
    s_clear = TIE_S - TIE_DIA
    rho_cc = col.Ast_tot / (bc * hc)
    ke = ((1 - sum_w2 / (6 * bc * hc)) * (1 - s_clear / (2 * bc))
          * (1 - s_clear / (2 * hc)) / (1 - rho_cc))
    assert col.mander_ke() == pytest.approx(ke, rel=1e-12)
    assert col.mander_ke() < 0.55, "a plain perimeter hoop must not approach 0.80"


def test_confinement_raises_eps_cu_above_unconfined(col):
    """The jacket ties buy back the ultimate strain the preload consumes."""
    assert col.confined("existing").eps_cu > 0.0038
    assert col.confined("existing").strength_gain > 1.0


# ----------------------------------------------- the three method decisions


def test_default_stress_block_is_the_single_block(col):
    """The default must be the conservative single block, not per_material."""
    hi = JacketedColumn(BE, HE, T, Concrete(FC_E), Concrete(45.0),
                        col.bars_e, col.bars_j, col.ties_j, cover_j=COVER_TIE,
                        n_sup_b=2, n_sup_h=2, db_long=DB)
    default = max(hi.interaction_code(n_points=60)["phiM"])
    single = max(hi.interaction_code(stress_block="single", n_points=60)["phiM"])
    assert default == pytest.approx(single, rel=1e-12)


def test_single_block_is_strictly_below_per_material(col):
    """Guards against BOTH branches being broken to the same rule, which the
    default-check above would miss.  Verified upstream by fault injection --
    do not merge these two tests."""
    hi = JacketedColumn(BE, HE, T, Concrete(FC_E), Concrete(45.0),
                        col.bars_e, col.bars_j, col.ties_j, cover_j=COVER_TIE,
                        n_sup_b=2, n_sup_h=2, db_long=DB)
    single = max(hi.interaction_code(stress_block="single", n_points=60)["phiM"])
    permat = max(hi.interaction_code(stress_block="per_material",
                                     n_points=60)["phiM"])
    assert single < permat * (1 - 1e-6), (
        "per_material must remain strictly higher (it is the unconservative rule)")


def test_shear_d_comes_from_the_extreme_tension_bar(col):
    """Not a guessed cover + tie + half-bar stack: the guess disagreed with the
    real layout by ~10 mm and desynchronised shear from the interaction
    diagram."""
    assert col.shear_capacity(Nu=0.0)["d"] == pytest.approx(
        max(b.y for b in col.bars))


# ------------------------------------------------------- behavioural guards


def test_capacity_outside_axial_range_is_nan_not_clamped(col):
    """Clamping would report a capacity the section does not have."""
    ia = col.interaction_code(n_points=60)
    assert math.isnan(capacity_at_P(ia["phiP"], ia["phiM"], 9e12))
    assert math.isnan(capacity_at_P(ia["phiP"], ia["phiM"], -9e12))


def test_degenerate_case_zero_jacket_equals_bare_column(col):
    """t = 0 with one material must reproduce the bare-column path exactly."""
    shifted = [Bar(b.y - T, b.z, b.area, b.steel, b.zone) for b in col.bars_e]
    # A bare t=0 section raises by design -- n_faces == 0 is a real error
    # everywhere except here, so the un-jacketed column has its own constructor.
    bare = JacketedColumn.bare_column(
        b_e=BE, h_e=HE, conc_e=Concrete(FC_E), conc_j=Concrete(FC_E),
        bars_e=shifted, bars_j=[], ties_j=col.ties_j,
        cover_j=COVER_TIE, db_long=DB)
    a = col.interaction_code(existing_only=True, n_points=60)
    b = bare.interaction_code(n_points=60)
    assert (capacity_at_P(a["phiP"], a["phiM"], 1000e3)
            == pytest.approx(capacity_at_P(b["phiP"], b["phiM"], 1000e3),
                             rel=1e-9))


def test_a_bare_section_must_use_its_own_constructor():
    """n_faces == 0 is a real error everywhere else."""
    ties = TieSet(Av_x=226.0, Av_z=226.0, s=100.0, fyt=275.0, db=12.0)
    with pytest.raises(ValueError):
        JacketedColumn(BE, HE, 0.0, Concrete(FC_E), Concrete(FC_J),
                       [], [], ties, cover_j=COVER_TIE)


def test_interface_reports_both_routes_and_picks_neither(col):
    """Crediting the 1.8 MPa cohesion term is an Engineer-of-Record decision,
    not the engine's."""
    it = col.interface_check(C_jacket=2.0e6, dP_jacket=1.4e6, Vu=300e3,
                             L_shear_span=1400, L_transfer=2800,
                             rho_v=0.001, fy_dowel=415.0)
    for key in ("phi_v_n_dowel", "phi_v_n_bond", "DCR_dowel", "DCR_bond",
                "ok_dowel_only", "ok_with_bond"):
        assert key in it
    assert "ok" not in it, \
        "the engine must not collapse the two routes into one verdict"


def test_continuity_removes_the_axial_interface_term(col):
    kw = dict(C_jacket=2.0e6, dP_jacket=1.4e6, Vu=300e3, L_shear_span=1400,
              L_transfer=2800, rho_v=0.001, fy_dowel=415.0)
    disc = col.interface_check(continuity="discontinuous", **kw)
    cont = col.interface_check(continuity="continuous", **kw)
    assert disc["v_axial"] > 0
    assert cont["v_axial"] == 0.0
    assert cont["v_u"] < disc["v_u"]


def test_not_roughened_surface_gets_no_cohesion_credit(col):
    """Table 16.4.4.2: the 1.8 MPa cohesion term needs a roughened interface,
    and ACI 318-19 Table 22.9.4.2 makes every mu a multiple of lambda."""
    kw = dict(C_jacket=2.0e6, dP_jacket=1.4e6, Vu=300e3, L_shear_span=1400,
              L_transfer=2800, rho_v=0.001, fy_dowel=415.0)
    rough = col.interface_check(surface="roughened_6mm", **kw)
    plain = col.interface_check(surface="not_roughened", **kw)
    assert plain["mu"] == pytest.approx(0.6 * 1.0)
    assert rough["mu"] == pytest.approx(1.0 * 1.0)
    assert plain["v_n_bond"] == 0.55, "no cohesion credit without roughening"
    assert rough["v_n_bond"] > plain["v_n_bond"]


def test_shear_friction_cap_branches_on_surface_condition(col):
    """ACI 318-19 Table 22.9.4.4: only monolithic or fully roughened take the
    least of 0.2f'c, 3.3 + 0.08f'c and 11.0 MPa.  Every other case takes the
    lesser of 0.2f'c and 5.5.  A single unbranched set credited the roughened
    ceiling to a not-roughened interface -- +40% at f'c = 55 MPa, which is the
    only strength range where the two differ at all.
    """
    leg = math.pi * TIE_DIA ** 2 / 4 * 2
    ties = TieSet(Av_x=leg, Av_z=leg, s=TIE_S, fyt=TIE_FY, db=TIE_DIA)
    hi = JacketedColumn(BE, HE, T, Concrete(55.0), Concrete(55.0),
                        col.bars_e, col.bars_j, ties, cover_j=COVER_TIE)
    kw = dict(C_jacket=2.0e6, dP_jacket=1.4e6, Vu=300e3, L_shear_span=1400,
              L_transfer=2800, rho_v=0.001, fy_dowel=415.0)
    rough = hi.interface_check(surface="roughened_6mm", **kw)
    plain = hi.interface_check(surface="not_roughened", **kw)
    assert rough["sf_cap"] == pytest.approx(min(0.2 * 55, 3.3 + 0.08 * 55, 11.0))
    assert plain["sf_cap"] == pytest.approx(min(0.2 * 55, 5.5))
    assert plain["sf_cap"] < rough["sf_cap"]


def test_preload_strain_scales_with_the_creep_coefficient(col):
    """Mechanics, not a code equation -- and it is why an unshored jacket does
    not deliver the monolithic capacity."""
    short = col.preload_strain(900e3, creep_factor=0.0)
    crept = col.preload_strain(900e3, creep_factor=2.0)
    assert crept == pytest.approx(3.0 * short, rel=1e-12)
    assert col.preload_strain(0.0, 2.0) == 0.0


def test_unshored_jacket_takes_less_than_a_monolithic_split(col):
    """Only the increment after casting is shared; P0 stays in the core."""
    ss = col.service_stress_split(P0=900e3, P_total=2240e3, creep_factor=2.0)
    assert ss["jacket_share"] < ss["jacket_share_monolithic"]
    assert ss["core_overstress"] > 1.0
    shored = col.service_stress_split(P0=0.0, P_total=2240e3, creep_factor=2.0)
    assert shored["core_overstress"] == pytest.approx(1.0, abs=1e-9)


def test_engine_stays_a_plain_library():
    """No web framework and no relative imports, so the workspace's own
    verification harness still applies to this copy unchanged."""
    from pathlib import Path
    import concretedesignpy.calculators.column_jacket as engine

    src = Path(engine.__file__).read_text()
    for forbidden in ("fastapi", "pydantic", "starlette", "flask",
                      "from ..", "from ."):
        assert forbidden not in src, \
            "vendored engine must not reference {!r}".format(forbidden)


def test_engine_hash_is_recorded():
    """Not a correctness check -- it proves the hash is computable, so a
    drifted vendor copy shows up in every design result."""
    assert len(ENGINE_SHA256) == 64
    assert ENGINE_SHA256 == ENGINE_SHA256.lower()


# -------------------------------------------------------------- model bounds


def test_two_adjacent_faces_is_refused_by_the_engine():
    """Unsymmetric about both principal axes -- the capacity is a biaxial
    SURFACE, not a curve.  TN-RET-002 TODO C2, BLOCKING."""
    ties = TieSet(2 * math.pi * 36, 2 * math.pi * 36, 100.0, 275.0, 12.0)
    col = JacketedColumn(400, 400, 100.0, Concrete(21.0), Concrete(28.0),
                         [], [], ties, cover_j=40.0,
                         t_top=100.0, t_bot=0.0, t_left=100.0, t_right=0.0)
    col.bars_j = rect_perimeter_bars(col.B, col.H, 62.0, 4, 4, 20.0,
                                     Steel(415.0), "jacket", 0.0, 0.0,
                                     faces=col.faces)
    assert col.doubly_unsymmetric
    with pytest.raises(NotImplementedError):
        col.interaction_code()


def test_an_unsymmetric_section_refuses_when_only_one_sense_exists():
    """Both diagrams or neither.

    J3-U is unsymmetric top-to-bottom, so the weaker sense governs -- but
    mirroring it puts both concretes at the extreme compression fibre, where
    the single stress block has no beta1.  Returning the one sense that happens
    to compute would report the section as checked when only half its capacity
    envelope exists.
    """
    ties = TieSet(2 * math.pi * 36, 2 * math.pi * 36, 100.0, 275.0, 12.0)
    col = JacketedColumn(400, 400, 100.0, Concrete(21.0), Concrete(28.0),
                         [], [], ties, cover_j=40.0,
                         t_top=0.0, t_bot=100.0, t_left=100.0, t_right=100.0)
    col.bars_j = rect_perimeter_bars(col.B, col.H, 62.0, 4, 4, 20.0,
                                     Steel(415.0), "jacket", 0.0, 0.0,
                                     faces=col.faces)
    col.bars_e = rect_perimeter_bars(400, 400, 50.0, 3, 3, 20.0, Steel(275.0),
                                     "existing", col.t_top, col.zc_e)
    for sense in ("as_built", "reversed"):
        with pytest.raises(NotImplementedError):
            col.interaction_code(sense=sense)


@pytest.mark.parametrize("method", [
    "core_dims", "rho_s_ties", "mander_ke", "stiffness_summary",
])
def test_four_sided_only_methods_refuse_a_partial_jacket(method):
    """Answering with four-sided behaviour is the silent extrapolation this
    module exists to prevent."""
    ties = TieSet(2 * math.pi * 36, 2 * math.pi * 36, 100.0, 275.0, 12.0)
    col = JacketedColumn(400, 400, 100.0, Concrete(21.0), Concrete(28.0),
                         [], [], ties, cover_j=40.0,
                         t_top=100.0, t_bot=100.0, t_left=100.0, t_right=0.0)
    with pytest.raises(NotImplementedError):
        attr = getattr(col, method)
        attr() if callable(attr) else attr


def test_monolithic_coefficients_refuse_a_one_sided_jacket():
    """Lampropoulos et al. studied four- and three-sided jackets only; the
    three-sided values are an UPPER BOUND for a smaller jacket, not an
    estimate, and must not be borrowed."""
    ties = TieSet(2 * math.pi * 36, 2 * math.pi * 36, 100.0, 275.0, 12.0)
    col = JacketedColumn(400, 400, 100.0, Concrete(21.0), Concrete(28.0),
                         [], [], ties, cover_j=40.0,
                         t_top=100.0, t_bot=0.0, t_left=0.0, t_right=0.0)
    mono = col.monolithic_coefficients(3200e3)
    assert mono["applicable"] is False
    assert mono["n_faces"] == 1
    assert mono["coefficients"] is None


def test_monolithic_coefficients_clamp_and_say_so():
    """The paper publishes no values outside the tabulated span, so
    extrapolating past an anchor would be inventing data."""
    bars_e = rect_perimeter_bars(BE, HE, 50.0, 3, 3, DB, Steel(FY_E),
                                 "existing", T)
    bars_j = rect_perimeter_bars(BE + 2 * T, HE + 2 * T, 62.0, 4, 4, DB,
                                 Steel(FY_J), "jacket", 0.0)
    leg = math.pi * TIE_DIA ** 2 / 4 * 2
    ties = TieSet(Av_x=leg, Av_z=leg, s=TIE_S, fyt=TIE_FY, db=TIE_DIA)
    col = JacketedColumn(BE, HE, T, Concrete(FC_E), Concrete(FC_J),
                         bars_e, bars_j, ties, cover_j=COVER_TIE)
    mono = col.monolithic_coefficients(50e3)      # nu far below the 0.1 anchor
    assert mono["applicable"] is True
    assert mono["out_of_range"] is True
    # clamped at the anchor, never extrapolated past it
    assert mono["coefficients"]["K_F"] <= 0.90
