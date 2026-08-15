"""Contract tests for the jacketing design boundary.

These check the contract, not the arithmetic — the arithmetic is covered by
``test_column_jacket_engine.py``.  What matters here is that the module refuses
nonsense, converts units exactly once, and never returns a number without the
provenance and advisories that make it auditable.

The reference case throughout is the TN-RET-001 worked example.
"""

import copy
import math

import pytest

from concretedesignpy.calculators.column_jacket_design import (
    build_jacketed_column,
    column_jacket_design,
    validate_jacket_request,
)


def reference_inputs(**overrides):
    p = {
        "existing": {"width": 400, "depth": 400, "fc": 21, "fy": 275,
                     "bar_dia": 20, "bars_per_face_width": 3,
                     "bars_per_face_depth": 3, "cover_to_bar_centre": 50},
        "jacket": {"thickness": 100, "fc": 28, "fy": 415, "bar_dia": 20,
                   "bars_per_face_width": 4, "bars_per_face_depth": 4,
                   "clear_cover_to_tie": 40, "tie_dia": 12, "tie_spacing": 100,
                   "tie_fy": 275, "tie_legs_each_way": 2,
                   "bars_restrained_per_face": 2},
        "demand": {"Pu": 3200, "Mu": 520, "Vu": 300, "clear_height": 2800},
        "construction": {"P0_at_casting": 900, "P_service": 2240,
                         "creep_coefficient": 2.0,
                         "continuity": "discontinuous", "dowel_dia": 16,
                         "dowels_per_row": 2, "dowel_row_spacing": 250,
                         "dowel_fy": 415},
    }
    p = copy.deepcopy(p)
    for k, v in overrides.items():
        p[k] = {**p[k], **v} if isinstance(v, dict) else v
    return p


@pytest.fixture(scope="module")
def result():
    return column_jacket_design(**reference_inputs())


def codes(res):
    return {a["code"] for a in res["advisories"]}


# ------------------------------------------------------------------ contract


def test_every_result_is_traceable(result):
    """A number without provenance is not auditable, which defeats the point."""
    p = result["provenance"]
    assert len(p["engine_sha256"]) == 64
    assert p["engine_version"].startswith("TN-RET-001-Rev01")
    assert "ACI 318-19" in p["code_basis"]
    assert "Engineer of Record" in p["disclaimer"]
    assert result["request_echo"]["existing"]["fc"] == 21, \
        "the validated inputs must be echoed so the result is self-contained"


def test_units_are_kN_and_kNm_at_the_boundary(result):
    """The engine works in N; this module must expose kN.  9737 kN, not 9.7e6.

    The failure is silent -- both values look plausible in a report -- so it is
    pinned explicitly.
    """
    assert result["axial"]["Po_jacketed"] == pytest.approx(9737.1, rel=1e-3)
    assert result["axial"]["phiPn_max_jacketed"] == pytest.approx(5063.3, rel=1e-3)
    assert result["shear"]["phiVn"] == pytest.approx(693.8, rel=1e-3)
    assert result["interaction"]["phiMn_at_Pu_jacketed"] == pytest.approx(
        535.25, rel=1e-3)


def test_defaults_are_applied_and_echoed():
    """An omitted construction block must not silently change the model."""
    inputs = reference_inputs()
    del inputs["construction"]
    req = validate_jacket_request(inputs["existing"], inputs["jacket"],
                                  inputs["demand"], None)
    c = req["construction"]
    assert c["continuity"] == "discontinuous"      # the conservative default
    assert c["creep_coefficient"] == 2.0
    assert c["P0_at_casting"] == 0.0
    assert req["jacket"]["sides"] == "four"
    assert req["jacket"]["spiral"] is False


def test_a_misspelt_field_is_refused_not_ignored():
    """Silently ignoring an unknown key would apply a default the caller
    believed they had overridden."""
    with pytest.raises(ValueError, match="unknown field"):
        column_jacket_design(**reference_inputs(
            jacket={"tie_spacing_mm": 100}))


def test_jacket_increases_capacity(result):
    a = result["axial"]
    assert a["Po_jacketed"] > a["Po_existing"]
    assert a["strength_gain"] == pytest.approx(2.78, rel=1e-2)


def test_existing_column_fails_and_jacketed_passes(result):
    assert result["interaction"]["utilisation_existing"] is None, \
        "existing column is beyond capacity at Pu -> None, not a clamped number"
    assert result["interaction"]["utilisation_jacketed"] < 1.0
    assert result["adequate"] is True
    assert result["governing_checks"] == []


def test_interaction_curves_are_returned(result):
    for key in ("design_existing", "design_jacketed", "nominal_jacketed"):
        pts = result["interaction"][key]
        assert len(pts) > 20
        assert all({"P", "M"} == set(p) for p in pts)
        assert all(isinstance(p["P"], float) for p in pts)


def test_confinement_ke_is_computed_not_assumed(result):
    """A plain perimeter hoop restrains only the four corner bars.  Assuming
    0.80 without crossties overstates f'cc by roughly 12%."""
    c = result["confinement"]
    assert c["ke_computed"] == pytest.approx(0.396, abs=0.01)
    assert c["ke_computed"] < 0.55
    assert c["eps_cu_core"] > c["eps_cu_unconfined"]


def test_interface_reports_both_routes(result):
    i = result["interface"]
    assert i["DCR_dowels_only"] > 1.0, "dowels alone should not close this design"
    assert i["DCR_with_bond"] < 1.0
    assert i["relies_on_bond"] is True
    assert i["aci562_interface_reinforcement_required"] is True


def test_geometry_matches_the_section_that_was_analysed(result):
    """The drawing payload and the shear check must see the same bars.

    ``d`` for shear is ``max(bar.y)``.  If the geometry block disagreed with
    it, the figure and the shear result would describe different sections.
    """
    g = result["geometry"]
    assert g["B"] == result["section"]["B"] and g["H"] == result["section"]["H"]
    assert len(g["bars"]) == (result["section"]["n_bars_existing"]
                              + result["section"]["n_bars_jacket"])
    assert {b["zone"] for b in g["bars"]} == {"existing", "jacket"}
    assert max(b["y"] for b in g["bars"]) == pytest.approx(
        result["shear"]["d_used"], rel=1e-9)
    ex = result["request_echo"]["existing"]
    assert g["interface_perimeter"] == pytest.approx(
        2 * (ex["width"] + ex["depth"]))


def test_build_jacketed_column_agrees_with_the_report(result):
    """The exposed section object must be the one the report was built from."""
    inputs = reference_inputs()
    col = build_jacketed_column(inputs["existing"], inputs["jacket"])
    assert col.B == result["section"]["B"]
    assert col.Po() / 1e3 == pytest.approx(result["axial"]["Po_jacketed"])
    assert len(col.bars) == len(result["geometry"]["bars"])


# ----------------------------------------------------------------- advisories


def test_advisories_are_present_and_typed(result):
    for expected in ("STIFFNESS_FEEDBACK", "JACKET_DISCONTINUOUS",
                     "RELIES_ON_INTERFACE_BOND", "SHRINKAGE_NOT_INCLUDED",
                     "UNSHORED_JACKET", "NO_CROSSTIES",
                     "FOUNDATION_NOT_CHECKED", "NSCP_CLAUSE_NUMBERING"):
        assert expected in codes(result), "missing advisory {}".format(expected)
    assert all(a["severity"] in ("info", "warning", "critical")
               for a in result["advisories"])
    assert all(a["message"] for a in result["advisories"])


def test_continuity_changes_the_interface_demand_and_advisory():
    disc = column_jacket_design(**reference_inputs())
    cont = column_jacket_design(**reference_inputs(
        construction={"continuity": "continuous"}))
    assert cont["interface"]["v_axial"] == 0.0
    assert cont["interface"]["v_u"] < disc["interface"]["v_u"]
    assert "JACKET_DISCONTINUOUS" not in codes(cont)


def test_shoring_removes_the_preload_advisory():
    shored = column_jacket_design(**reference_inputs(
        construction={"P0_at_casting": 0}))
    assert "UNSHORED_JACKET" not in codes(shored)
    assert shored["preload"]["core_overstress"] == pytest.approx(1.0, abs=1e-9)


def test_crossties_raise_ke_and_drop_the_advisory():
    tied = column_jacket_design(**reference_inputs(
        jacket={"bars_restrained_per_face": 4}))
    base = column_jacket_design(**reference_inputs())
    assert tied["confinement"]["ke_computed"] > base["confinement"]["ke_computed"]
    assert tied["confinement"]["fcc_core"] > base["confinement"]["fcc_core"]
    assert "NO_CROSSTIES" not in codes(tied)


def test_shear_critical_case_is_flagged():
    """Wide tie spacing drops Vs; the failure mode must flip and say so."""
    r = column_jacket_design(**reference_inputs(
        jacket={"tie_spacing": 400, "tie_dia": 10}))
    assert r["shear"]["failure_mode"] == "shear-critical"
    assert "SHEAR_CRITICAL" in codes(r)
    assert "shear-critical failure mode" in r["governing_checks"]
    assert r["adequate"] is False


# ----------------------------------------------------------------- validation


@pytest.mark.parametrize("patch", [
    {"existing": {"fc": 0}},
    {"existing": {"fc": 500}},
    {"existing": {"width": 10}},
    {"existing": {"bars_per_face_width": 1}},
    {"existing": {"bar_dia": 200}},
    {"jacket": {"thickness": 5}},
    {"jacket": {"tie_spacing": 0}},
    {"jacket": {"fy": 5000}},
    {"jacket": {"tie_dia": 2}},
    {"jacket": {"sides": "five"}},
    {"demand": {"clear_height": 10}},
    {"demand": {"Mu": -5}},
    {"construction": {"creep_coefficient": 99}},
    {"construction": {"continuity": "maybe"}},
    {"construction": {"dowel_row_spacing": 0}},
])
def test_nonsense_input_is_rejected(patch):
    """The engine would return a confident number for an impossible section.
    Validation is where that gets refused."""
    with pytest.raises(ValueError):
        column_jacket_design(**reference_inputs(**patch))


def test_missing_required_field_is_rejected():
    inputs = reference_inputs()
    del inputs["jacket"]["tie_spacing"]
    with pytest.raises(ValueError, match="jacket.tie_spacing is required"):
        column_jacket_design(**inputs)


@pytest.mark.parametrize("bad", [float("nan"), float("inf"), -float("inf")])
def test_nan_and_inf_demand_are_refused(bad):
    """A garbage Pu flips comparisons silently -- including the two mandatory
    ACI 562 interface requirements, which are `>` tests."""
    with pytest.raises(ValueError, match="finite"):
        column_jacket_design(**reference_inputs(demand={"Pu": bad}))


def test_cover_that_collapses_the_lever_arm_is_rejected():
    """A 400 mm column with 160 mm cover leaves 80 mm between outer bar rows --
    under 5 bar diameters.  The engine would integrate it and return a
    confident, meaningless flexural capacity."""
    with pytest.raises(ValueError, match="lever arm"):
        column_jacket_design(**reference_inputs(
            existing={"cover_to_bar_centre": 160}))


def test_jacket_too_thin_for_its_own_reinforcement_is_rejected():
    with pytest.raises(ValueError, match="cannot house"):
        column_jacket_design(**reference_inputs(
            jacket={"thickness": 55, "bar_dia": 25, "tie_dia": 12}))


def test_P0_greater_than_Pu_is_rejected():
    with pytest.raises(ValueError, match="P0_at_casting exceeds Pu"):
        column_jacket_design(**reference_inputs(
            construction={"P0_at_casting": 9999}))


def test_a_net_tension_demand_is_currently_unreachable():
    """Records an INHERITED bound, and is not an endorsement of it.

    ``P0_at_casting`` is >= 0 and the sanity check refuses ``P0 > Pu``, so a
    negative (net tension) ``Pu`` is refused for every legal P0 -- including a
    fully shored column at P0 = 0, which is a physically sensible input. The
    ``Demand`` bound admits ``Pu > -1e6``, so the two rules disagree about
    whether net tension is in scope.

    The refusal is the conservative direction, so it is left as it stands:
    widening a bound is a decision for the Engineer of Record, not a fix to
    make in passing. If net tension is to be supported, this is the bound to
    revisit -- ``P0_at_casting > max(Pu, 0)`` is the likely intent.
    """
    with pytest.raises(ValueError, match="P0_at_casting exceeds Pu"):
        column_jacket_design(**reference_inputs(
            demand={"Pu": -500, "Mu": 100},
            construction={"P0_at_casting": 0}))


def test_two_adjacent_faces_is_refused_not_approximated():
    """Unsymmetric about both axes -- its capacity is a biaxial SURFACE.
    Returning the uniaxial curve would silently ignore M_z."""
    with pytest.raises(ValueError) as exc:
        column_jacket_design(**reference_inputs(
            jacket={"sides": "two_adjacent"}))
    assert "biaxial" in str(exc.value) and "C2" in str(exc.value)


def test_four_sided_is_the_default_and_still_computes():
    inputs = reference_inputs()
    assert "sides" not in inputs["jacket"], \
        "the worked example should rely on the default"
    assert column_jacket_design(**inputs)["adequate"] is True
    inputs["jacket"]["sides"] = "four"
    assert column_jacket_design(**inputs)["adequate"] is True


# ------------------------------------------------------------ partial jackets


@pytest.mark.parametrize("sides", ["one", "two", "three"])
def test_partial_jackets_compute_but_declare_what_is_missing(sides):
    """One, two opposite and three sides compute their geometry, axial strength
    and P-M interaction.

    What must NOT happen is a partial jacket quietly receiving four-sided
    confinement, shear, interface, stiffness and detailing.  The engine refuses
    those, the result records the refusals, and ``adequate`` stays false -- a
    section whose shear was never checked has not been shown adequate.
    """
    d = column_jacket_design(**reference_inputs(jacket={"sides": sides}))

    # the geometry is real, and smaller than the four-sided section
    assert d["geometry"]["B"] <= 600 and d["geometry"]["H"] <= 600
    assert d["axial"]["phiPn_max_jacketed"] < 5063
    assert len(d["geometry"]["bars"]) < 20

    # and the gaps are declared, not filled in
    assert d["unavailable"], "{}-sided declared nothing unavailable".format(sides)
    assert d["shear"] is None and d["confinement"] is None
    assert d["interface"] is None and d["stiffness"] is None
    assert d["detailing"] == []
    assert d["adequate"] is False, "unchecked is not the same as passed"
    assert "PARTIAL_JACKET_INCOMPLETE" in codes(d)


def test_partial_jacket_advisory_does_not_claim_what_is_not_computed():
    """The one message that exists to disclose omissions must not conceal them.

    ``induced_eccentricity`` exists in the engine and has no call site here, so
    P*e is never added to M_u.  On TN-RET-002's one-sided case the reported DCR
    is 0.82 where the note's own table gives 1.09.
    """
    d = column_jacket_design(**reference_inputs(jacket={"sides": "one"}))
    adv = next(a for a in d["advisories"]
               if a["code"] == "PARTIAL_JACKET_INCOMPLETE")
    assert "does NOT add it to M_u" in adv["message"]
    assert "as-built sense" in adv["message"]
    assert "0.82" in adv["message"] and "1.09" in adv["message"], \
        "the consequence must be quantified"


def test_a_one_sided_jacket_gets_no_monolithic_coefficients():
    """The three-sided values are an UPPER BOUND for a smaller jacket, not an
    estimate, and must not be borrowed."""
    d = column_jacket_design(**reference_inputs(jacket={"sides": "one"}))
    assert d["monolithic"]["applicable"] is False
    assert d["monolithic"]["n_faces"] == 1
    assert "NO_MONOLITHIC_COEFFICIENTS" in codes(d)


# -------------------------------------------------- monolithic coefficients


def test_monolithic_coefficient_is_reported_and_not_applied(result):
    """A jacketed column does not reach the strength of a monolithic section.

    Both results must be present and neither may replace the other.  Whether
    K_F multiplies phi*M_n or M_n is TN-RET-001 TODO C8/C9 and BLOCKING.
    """
    m = result["monolithic"]
    assert m["applicable"] and m["K_F"] == pytest.approx(0.7429, abs=5e-4)

    # both, side by side -- the uncorrected value is NOT overwritten
    assert m["phiMn_monolithic"] == pytest.approx(535.25, rel=1e-3)
    assert m["phiMn_x_KF"] == pytest.approx(397.6, rel=1e-3)
    assert m["utilisation_monolithic"] == pytest.approx(0.9715, rel=1e-3)
    assert m["utilisation_x_KF"] == pytest.approx(1.308, rel=1e-3)

    # the uncorrected numbers elsewhere in the report are untouched
    assert result["interaction"]["phiMn_at_Pu_jacketed"] == pytest.approx(
        535.25, rel=1e-3)

    # and the caller is told, loudly, that a correction exists
    crit = next(a for a in result["advisories"]
                if a["code"] == "MONOLITHIC_COEFFICIENT_NOT_APPLIED")
    assert crit["severity"] == "critical"
    assert "1.31" in crit["message"] or "1.308" in crit["message"]


def test_monolithic_carries_the_direction_of_conservatism_rule(result):
    """K_F is a strength REDUCTION, so it is safe only where a smaller strength
    is the safe answer.

    Applying it to a capacity-design quantity lowers V_e, which looks
    conservative and is not.  Vp_at_hinging is reported WITHOUT K_F.
    """
    assert "K_NOT_FOR_CAPACITY_DESIGN" in codes(result)
    adv = next(a for a in result["advisories"]
               if a["code"] == "K_NOT_FOR_CAPACITY_DESIGN")
    assert adv["severity"] == "critical"
    assert result["monolithic"]["do_not_apply_to"], \
        "the exclusions must be machine-readable"
    kf = result["monolithic"]["K_F"]
    assert result["shear"]["Vp_at_hinging"] == pytest.approx(652.5, rel=1e-2)
    assert result["shear"]["Vp_at_hinging"] != pytest.approx(652.5 * kf,
                                                             rel=1e-3)


def test_deformation_penalty_is_reported_as_procedure_dependent(result):
    """K_du/K_dy costs ductility while plastic rotation barely moves, so an
    ASCE 41 linear procedure absorbs the whole penalty and a nonlinear one
    absorbs almost none."""
    assert result["monolithic"]["ductility_factor"] == pytest.approx(
        0.8205, abs=5e-4)
    assert "DEFORMATION_PENALTY_IS_PROCEDURE_DEPENDENT" in codes(result)


def test_an_overloaded_column_reports_no_capacity_instead_of_crashing():
    """P_u beyond the section's axial range has no phi*M_n at all.

    That must surface as "no capacity" -- None, a failing DCR and the
    advisories intact -- not as a clamped number and not as a formatting crash
    in the advisory that exists to quantify the correction.
    """
    r = column_jacket_design(**reference_inputs(demand={"Pu": 9000}))
    assert r["interaction"]["phiMn_at_Pu_jacketed"] is None
    assert r["interaction"]["utilisation_jacketed"] is None
    assert r["adequate"] is False
    assert "P-M interaction" in r["governing_checks"]

    adv = next(a for a in r["advisories"]
               if a["code"] == "MONOLITHIC_COEFFICIENT_NOT_APPLIED")
    assert "outside the section's axial range" in adv["message"]
    assert r["monolithic"]["K_F"] is not None, \
        "the coefficient itself is still defined and must still be reported"


def test_large_jacket_caution_fires_past_r_of_two():
    """Past r = 2 the tabulated behaviour turns non-monotonic."""
    base = column_jacket_design(**reference_inputs())
    assert base["monolithic"]["large_jacket_caution"] is False   # r = 1.25
    big = column_jacket_design(**reference_inputs(
        jacket={"thickness": 200}))                              # r = 2.25
    assert big["monolithic"]["r"] > 2.0
    assert big["monolithic"]["large_jacket_caution"] is True
    assert "LARGE_JACKET_COEFFICIENTS_NON_MONOTONIC" in codes(big)


# --------------------------------------------------------------------- QAQC


def test_qaqc_all_checks_pass_on_the_reference_case(result):
    """The QAQC block re-derives reported values along independent arithmetic
    paths. On the reference case every one must agree."""
    q = result["qaqc"]
    assert q["all_pass"] is True
    assert q["n_pass"] == len(q["checks"])
    names = {c["name"] for c in q["checks"]}
    for expected in ("Nominal axial strength Po", "phi*Mn at Pu",
                     "Design shear strength phi*Vn", "Interface demand v_u",
                     "Confined strength f'cc"):
        assert expected in names, "missing QAQC check {}".format(expected)


@pytest.mark.parametrize("patch", [
    {},                                           # reference
    {"jacket": {"spiral": True}},                 # phi/alpha switch
    {"jacket": {"tie_spacing": 400, "tie_dia": 10}},  # shear-critical
    {"existing": {"fc": 90, "fy": 600}},          # fy above the 550 cap
    {"demand": {"Pu": 9000, "Mu": 520, "Vu": 300,
                "clear_height": 2800}},           # beyond axial capacity
])
def test_qaqc_survives_the_variants(patch):
    """Independent recomputation must track the report across the parameter
    space -- including the fy > 550 cap and a Pu with no capacity (both the
    reported value and the recomputation must be None there, and None == None
    is a pass, not a crash)."""
    r = column_jacket_design(**reference_inputs(**patch))
    assert r["qaqc"]["all_pass"] is True, [
        c["name"] for c in r["qaqc"]["checks"] if not c["pass"]]


def test_qaqc_adapts_to_a_partial_jacket():
    """Checks whose subject was refused (shear, interface, confinement, and
    the four-sided cage arithmetic) must be absent, not failed."""
    r = column_jacket_design(**reference_inputs(jacket={"sides": "three"}))
    q = r["qaqc"]
    assert q["all_pass"] is True
    names = {c["name"] for c in q["checks"]}
    assert "Design shear strength phi*Vn" not in names
    assert "Interface demand v_u" not in names
    assert "Nominal axial strength Po" not in names, \
        "the four-sided bar-cage recomputation does not apply to a partial cage"
    assert "Overall width B" in names and "phi*Mn at Pu" in names


def test_qaqc_catches_a_corrupted_report():
    """The block must actually be able to fail: perturb a reported value the
    way a units bug would and the matching check must flip to FAIL."""
    import importlib

    m = importlib.import_module(
        "concretedesignpy.calculators.column_jacket_design")
    r = column_jacket_design(**reference_inputs())
    # simulate the classic silent failure: a value scaled by 10 (the same
    # class of defect as N-mm -> kN-m dividing by 1e7 instead of 1e6)
    sec = dict(r["section"])
    ax = dict(r["axial"])
    ax["Po_jacketed"] = ax["Po_jacketed"] / 10.0
    q = m._qaqc_block(r["request_echo"], sec, ax, r["interaction"],
                      r["shear"], r["confinement"], r["interface"],
                      r["request_echo"]["demand"])
    assert q["all_pass"] is False
    failed = {c["name"] for c in q["checks"] if not c["pass"]}
    assert "Nominal axial strength Po" in failed


# ---------------------------------------------------------------- regression


@pytest.mark.parametrize("path,expected,rel", [
    (("interaction", "phiMn_at_Pu_jacketed"), 535.25, 1e-3),   # kN.m
    (("interface", "v_u"), 1.3778, 1e-3),                      # MPa
    (("stiffness", "EI_ratio"), 5.6910, 1e-4),
    (("axial", "Po_jacketed"), 9737.1, 1e-4),                  # kN
])
def test_regression_pins(path, expected, rel, result):
    """REGRESSION PINS, not independent checks.  They would happily pin a bug.
    If one fails, work out which behaviour changed before updating the number.

    Every value here is one PUBLISHED by the upstream jacketing service, so
    they also confirm the port did not move the arithmetic.  Do not add a pin
    by pasting whatever this module currently returns -- that records a run,
    not a reference.  Confined strength is deliberately absent: it is covered
    properly by ``test_confinement_matches_independent_recomputation``, which
    re-derives f'cc from the fib B14 equation.
    """
    got = result[path[0]][path[1]]
    assert got == pytest.approx(expected, rel=rel)


def test_result_is_json_serialisable(result):
    """The report is meant to cross a process boundary (web route, file, log).
    A numpy scalar or a NaN would serialise wrong or not at all."""
    import json

    text = json.dumps(result, allow_nan=False)
    assert "NaN" not in text and "Infinity" not in text
    assert json.loads(text)["adequate"] is True


def test_no_reported_number_is_a_nan(result):
    def walk(node):
        if isinstance(node, dict):
            for v in node.values():
                walk(v)
        elif isinstance(node, list):
            for v in node:
                walk(v)
        elif isinstance(node, float):
            assert not math.isnan(node) and not math.isinf(node), \
                "a NaN or inf reached the report"

    walk(result)
