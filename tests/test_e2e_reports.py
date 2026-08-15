# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
End-to-end tests for the three calculation sheets
=================================================

These go through the whole stack the browser goes through: a real Flask
app, a real JSON POST, real routing, real serialisation. They are not unit
tests of the report builders -- ``tests/test_beam_*.py`` cover the
calculators and ``tests/js/calcsheet.test.js`` covers the renderer. What is
tested HERE is everything between those two, which is exactly where the two
defects this repo has actually shipped both lived:

* ``beam_torsion.html`` read ``result.t_threshold`` while the module
  returned ``torsion_threshold``, so the UI rendered nothing;
* ``joint-shear.html`` read ``result.vu_joint`` and ``result.vn_joint``,
  neither of which is a key its module returns.

Both are invisible to a unit test on either side. A contract test between
them is the only thing that catches them, so the payload keys the JS reads
are asserted against the payload keys the server sends.
"""

import json
import math
import re
from pathlib import Path

import pytest

from concretedesignpy.webapp import create_app

IN = 25.4
KSI = 6.894757
KIP = 4.44822
KIPFT = 1.35582

REPO = Path(__file__).resolve().parents[1]


@pytest.fixture
def client():
    app = create_app()
    app.config["TESTING"] = True
    with app.test_client() as c:
        yield c


# ── payloads that mirror what each form actually submits ──────────────

FLEXURE = dict(rebar_list=[{"d": 500, "diam": 25.4692, "num": 3}],
               fc=20, fy=420, b=250, h=565, es=200000, mu_demand=200)

SHEAR = dict(fc=25, fyv=300, bw=300, h=675, d=610, vu=280,
             s_chosen=200, n_legs=2, db_stirrup=12.7, cc=40, c=65,
             nu=0, phi=0.75)

JOINT = dict(ve=81.8 * KIP, as1=4.36 * IN * IN, n_bars1=1,
             as2=2.24 * IN * IN, n_bars2=1, fy=60 * KSI, fc=4 * KSI,
             beam_width=24 * IN, joint_depth=24 * IN, column_width=24 * IN,
             perpendicular_dist=0, joint_config=1, lamda=1.0)

ENDPOINTS = [
    ("/api/beam/flexure-report", FLEXURE, "Beam Flexural Capacity"),
    ("/api/beam/shear-report", SHEAR, "Beam Shear Design"),
    ("/api/joint/shear-report", JOINT, "Beam-Column Joint Shear"),
]


def post(client, url, payload):
    r = client.post(url, json=payload)
    assert r.status_code == 200, r.get_data(as_text=True)
    body = r.get_json()
    assert body["status"] == "success", body
    return body["result"]


# ══════════════════════════════════════════════════════════════════════
# 1 — the round trip
# ══════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("url,payload,title", ENDPOINTS)
def test_report_round_trips_over_http(client, url, payload, title):
    d = post(client, url, payload)
    assert d["title"] == title
    assert d["basis"]
    assert d["sections"]
    assert d["summary"]
    assert d["qaqc"]["checks"]
    assert d["provenance"]["report_version"]


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_response_is_strict_json_with_no_nan(client, url, payload, _t):
    """A NaN serialises as the bare token NaN, which is not valid JSON and
    which every downstream parser then disagrees about. The report builder
    maps non-finite values to null; this asserts none slips through."""
    raw = client.post(url, json=payload).get_data(as_text=True)
    # The bare tokens, unquoted -- a quoted "Infinity" is a legal JSON
    # string, an unquoted one is Python's json module emitting something no
    # other parser accepts.
    assert not re.search(r"[:\[,]\s*(NaN|-?Infinity)\b", raw)
    json.loads(raw)  # strict mode below would accept them; this would not
    json.loads(raw, parse_constant=lambda c: (_ for _ in ()).throw(
        AssertionError("non-finite constant %r in the response" % c)))


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_every_qaqc_check_passes_end_to_end(client, url, payload, _t):
    q = post(client, url, payload)["qaqc"]
    failed = [c["name"] for c in q["checks"] if not c["pass"]]
    assert not failed, "QAQC failures: %s" % failed
    assert q["all_pass"] is True
    assert q["n_pass"] == q["n_total"] == len(q["checks"])


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_every_qaqc_check_states_its_method(client, url, payload, _t):
    """A QAQC row without a stated independent method is decoration."""
    for c in post(client, url, payload)["qaqc"]["checks"]:
        assert c["name"] and c["method"], c
        assert len(c["method"].split()) >= 3, (
            "%r is a formula, not a stated method" % c["method"])


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_every_step_carries_a_reference(client, url, payload, _t):
    d = post(client, url, payload)
    for sec in d["sections"]:
        assert sec["heading"]
        assert sec["steps"], sec["heading"]
        for s in sec["steps"]:
            assert set(s) >= {"ref", "desc", "eq", "result", "status"}
            assert s["ref"], s
            assert s["desc"], s
            assert s["result"] is not None, s


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_advisories_and_unavailable_are_never_empty(client, url, payload, _t):
    """Every one of these sheets has real omissions. A sheet reporting none
    has lost them, which is worse than having them."""
    d = post(client, url, payload)
    assert d["advisories"], "a sheet with no advisories has lost them"
    assert d["unavailable"], "a sheet claiming to check everything is lying"
    for a in d["advisories"]:
        assert a["code"] and a["text"]
        assert a["severity"] in ("critical", "warning", "info")
    for u in d["unavailable"]:
        assert u["check"] and u["why"] and u["clause"]


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_complete_is_false_while_checks_are_unperformed(client, url, payload, _t):
    """`adequate` and `complete` are two fields on purpose. Collapsing them
    is what lets a sheet say 'adequate' about a section whose As,min was
    never looked at."""
    d = post(client, url, payload)
    assert d["complete"] is False
    assert d["adequate"] is True
    assert d["adequate"] != d["complete"]


# ══════════════════════════════════════════════════════════════════════
# 2 — the numbers survive the round trip
# ══════════════════════════════════════════════════════════════════════

def test_flexure_sheet_reproduces_wm_example_4_1m(client):
    """W&M Example 4-1M, printed 148: Mn = 273 kN.m, a = 151, c = 178.
    The same pin as tests/test_beam_flexure.py, asserted through HTTP so a
    serialisation or routing change cannot quietly move it."""
    d = post(client, "/api/beam/flexure-report", FLEXURE)
    raw = d["raw"]
    assert abs(raw["mn"] / 273.0 - 1) < 0.01
    assert abs(raw["a"] / 151.0 - 1) < 0.01
    assert abs(raw["neutral_axis"] / 178.0 - 1) < 0.01
    assert raw["phi"] == pytest.approx(0.90)


def test_shear_sheet_reproduces_wm_example_6_1m(client):
    """W&M Example 6-1M, printed 297-300: Vc = 153 kN."""
    d = post(client, "/api/beam/shear-report", SHEAR)
    assert abs(d["raw"]["vc"] / 153.0 - 1) < 0.01


def test_joint_sheet_reproduces_wm_example_19_3(client):
    """W&M Example 19-3, printed 1076-1078: Vj = 413 kips, Aj = 576 in2."""
    d = post(client, "/api/joint/shear-report", JOINT)
    raw = d["raw"]
    assert abs(raw["v_joint"] / (413 * KIP) - 1) < 0.01
    assert abs(raw["aj"] / (576 * IN * IN) - 1) < 0.01


def test_joint_sheet_shows_the_column_cap_doing_its_work(client):
    """A 500 mm spandrel on a 400 x 600 column. Aj must be the column
    section, and the sheet must SHOW the min() that produced it -- a fix
    the reader cannot see on the sheet is a fix nobody checks."""
    p = dict(JOINT, beam_width=500.0, joint_depth=600.0, column_width=400.0,
             fy=420, fc=28, ve=0, as1=2000, n_bars1=1, as2=1000, n_bars2=1)
    d = post(client, "/api/joint/shear-report", p)
    assert d["raw"]["aj"] == 240000.0
    assert d["raw"]["joint_width"] == 400.0
    text = json.dumps(d)
    assert "R15.5.2.2" in text
    assert "min(b_{col}" in text or "\\min(b_{col}" in text


# ══════════════════════════════════════════════════════════════════════
# 3 — verdict states
# ══════════════════════════════════════════════════════════════════════

def test_flexure_without_a_demand_draws_no_verdict(client):
    """No demand means no verdict. Not a pass, not a fail, and above all
    not a D/C of zero."""
    p = dict(FLEXURE)
    p.pop("mu_demand")
    d = post(client, "/api/beam/flexure-report", p)
    assert d["adequate"] is None
    assert d["has_demand"] is False
    assert d["governing_checks"] == []
    assert d["summary"][-1]["value"] == "not evaluated"


def test_flexure_with_a_blank_demand_string_is_treated_as_no_demand(client):
    """The form sends "" for an empty optional field, not null."""
    d = post(client, "/api/beam/flexure-report", dict(FLEXURE, mu_demand=""))
    assert d["adequate"] is None


def test_flexure_over_capacity_fails_and_names_the_clause(client):
    d = post(client, "/api/beam/flexure-report", dict(FLEXURE, mu_demand=400))
    assert d["adequate"] is False
    assert len(d["governing_checks"]) == 1
    g = d["governing_checks"][0]
    assert g["dc"] > 1.0
    assert "9.5.1.1" in g["clause"]


def test_overstressed_shear_reports_every_failing_check(client):
    """Above the 22.5.1.2 limit the section is inadequate. The sheet must
    say so, name the clause, and not offer a spacing that reads usable."""
    d = post(client, "/api/beam/shear-report", dict(SHEAR, vu=900))
    assert d["adequate"] is False
    checks = [g["check"] for g in d["governing_checks"]]
    assert "Section 22.5.1.2 cross-sectional limit" in checks
    assert any("22.5.1.2" in g["clause"] for g in d["governing_checks"])
    assert d["raw"]["shear_status"] == "UNSAFE"


def test_joint_over_capacity_fails(client):
    p = dict(JOINT, beam_width=300.0, joint_depth=300.0, column_width=300.0,
             fc=21.0)
    d = post(client, "/api/joint/shear-report", p)
    assert d["adequate"] is False
    assert d["raw"]["status"] == "REINFORCEMENT REQUIRED"


def test_axial_tension_reaches_the_shear_sheet_as_a_critical_advisory(client):
    d = post(client, "/api/beam/shear-report", dict(SHEAR, nu=-150))
    codes = {a["code"]: a["severity"] for a in d["advisories"]}
    assert codes.get("V-TENSION") == "critical"
    assert d["raw"]["vc"] < post(client, "/api/beam/shear-report",
                                 SHEAR)["raw"]["vc"]


def test_compression_controlled_section_raises_the_prohibition_advisory(client):
    """Section 9.3.3.1 prohibits a compression-controlled beam outright.
    This sheet applies the phi penalty and does not apply the prohibition,
    so it must say that in the loudest band it has."""
    p = dict(FLEXURE, rebar_list=[{"d": 500, "diam": 32.0, "num": 8}],
             mu_demand=100)
    d = post(client, "/api/beam/flexure-report", p)
    if d["raw"]["classification"] != "tension-controlled":
        codes = {a["code"]: a["severity"] for a in d["advisories"]}
        assert codes.get("F-BRITTLE") == "critical"
        assert d["advisories"][0]["code"] == "F-BRITTLE"


# ══════════════════════════════════════════════════════════════════════
# 4 — the joint sheet cannot be talked into an illegal input
# ══════════════════════════════════════════════════════════════════════

@pytest.mark.parametrize("bad,frag", [
    ({"lamda": 0.9}, "lamda"),
    ({"joint_config": 4}, "joint_config"),
    ({"column_width": 0}, "column_width"),
])
def test_joint_report_rejects_illegal_input(client, bad, frag):
    r = client.post("/api/joint/shear-report", json=dict(JOINT, **bad))
    assert r.status_code == 400
    assert frag in r.get_json()["message"]


def test_joint_report_ignores_a_phi_in_the_payload(client):
    """Section 21.2.4.4 fixes phi at 0.85. A payload carrying phi = 0.9
    must not move the answer -- the route does not read the field."""
    a = post(client, "/api/joint/shear-report", JOINT)
    b = post(client, "/api/joint/shear-report", dict(JOINT, phi=0.9))
    assert a["raw"]["phi_vn"] == b["raw"]["phi_vn"]
    assert b["raw"]["phi"] == 0.85


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_missing_required_field_is_a_400_not_a_500(client, url, payload, _t):
    for key in list(payload)[:3]:
        p = {k: v for k, v in payload.items() if k != key}
        r = client.post(url, json=p)
        assert r.status_code == 400, "%s without %s" % (url, key)
        assert r.get_json()["status"] == "error"


@pytest.mark.parametrize("url,payload,_t", ENDPOINTS)
def test_a_non_numeric_field_is_a_400_not_a_500(client, url, payload, _t):
    key = "fc"
    r = client.post(url, json=dict(payload, **{key: "banana"}))
    assert r.status_code == 400
    assert r.get_json()["status"] == "error"


# ══════════════════════════════════════════════════════════════════════
# 5 — the page/payload contract
#
# The two UI defects this repo shipped were both a template reading a key
# its module never returned. Neither a Python test nor a Jest test can see
# that; only a contract test across the boundary can.
# ══════════════════════════════════════════════════════════════════════

CALC_PAGES = [
    ("/calculator/beam-moment", "/api/beam/flexure-report"),
    ("/calculator/beam-shear", "/api/beam/shear-report"),
    ("/calculator/joint-shear", "/api/joint/shear-report"),
]


@pytest.mark.parametrize("page,endpoint", CALC_PAGES)
def test_calculator_page_loads_the_sheet_assets(client, page, endpoint):
    html = client.get(page).get_data(as_text=True)
    assert "css/calcsheet.css" in html
    assert "js/calcsheet.js" in html
    assert 'id="calc-sheet"' in html
    assert "buildSheetPayload" in html
    assert endpoint in html


@pytest.mark.parametrize("path", ["/static/css/calcsheet.css",
                                  "/static/js/calcsheet.js"])
def test_sheet_assets_are_served(client, path):
    assert client.get(path).status_code == 200


def test_renderer_reads_only_keys_the_server_sends(client):
    """The contract test. Every top-level `d.<key>` the renderer reads must
    exist in every report the server produces.

    This is the test that would have caught `result.t_threshold` and
    `result.vu_joint` on the day they were written.
    """
    js = (REPO / "concretedesignpy/webapp/static/js/calcsheet.js").read_text()
    read = set(re.findall(r"\bd\.([a-z_]+)\b", js))
    # `d` is also the parameter name in helpers that take a sub-object;
    # these are the report-level reads.
    assert read, "no key reads found -- has the renderer been rewritten?"

    for url, payload, _t in ENDPOINTS:
        d = post(client, url, payload)
        missing = sorted(k for k in read if k not in d)
        assert not missing, "%s renders keys the server never sends: %s" % (
            url, missing)


def test_navbar_version_tracks_the_package(client):
    """The badge was the literal string 'v0.7 | NSCP 2015 / ACI 318-19' and
    went stale. It is injected now, so this pins that it stays injected."""
    import concretedesignpy
    html = client.get("/").get_data(as_text=True)
    assert "v%s" % concretedesignpy.__version__ in html
    assert "v0.7 |" not in html
    # and the edition label must not claim 318-19 for the whole package
    # (the badge moved from the old navbar to the sidebar foot in the
    # solver-workbench redesign; the pin follows the element)
    m = re.search(r'foot-version[^>]*>([^<]+)<', html)
    assert m is not None, "version badge element missing"
    assert "318-19" not in m.group(1)


def test_js_fixture_shape_matches_the_server(client):
    """The Jest fixtures must not drift from the real payload, or the Jest
    suite passes while the page is broken."""
    fixture = (REPO / "tests/js/fixtures.js").read_text()
    fixture_keys = set(re.findall(r"^\s{8}([a-z_]+):", fixture, re.M))
    assert fixture_keys, "fixture shape not detected"

    d = post(client, "/api/beam/flexure-report", FLEXURE)
    unknown = sorted(k for k in fixture_keys if k not in d)
    assert not unknown, "fixture invents keys the server does not send: %s" % unknown

    # and the fixture must cover the keys the renderer depends on
    for required in ("sections", "summary", "qaqc", "advisories",
                     "unavailable", "governing_checks", "provenance",
                     "adequate", "complete", "title", "basis"):
        assert required in fixture_keys, required
        assert required in d, required


# ══════════════════════════════════════════════════════════════════════
# 6 — the printed sheet
# ══════════════════════════════════════════════════════════════════════

def test_stylesheet_declares_a4():
    css = (REPO / "concretedesignpy/webapp/static/css/calcsheet.css").read_text()
    assert "@page" in css
    assert "size: A4 portrait" in css
    assert "print-color-adjust: exact" in css, (
        "without this Chrome prints the PASS chips white")


def test_stylesheet_hides_the_chrome_when_printing():
    css = (REPO / "concretedesignpy/webapp/static/css/calcsheet.css").read_text()
    printed = css.split("@media print")[1]
    for chrome in (".navbar", ".calc-form", ".cs-actions", ".footer"):
        assert chrome in printed, "%s would print" % chrome


def test_stylesheet_controls_page_breaks():
    css = (REPO / "concretedesignpy/webapp/static/css/calcsheet.css").read_text()
    printed = css.split("@media print")[1]
    assert "break-inside: avoid" in printed
    assert "break-after: avoid" in printed
    assert "table-header-group" in printed, (
        "a QAQC table spanning a page break needs its header repeated")


def test_screen_width_equals_the_a4_content_width():
    """182mm is 210mm less two 14mm margins. If the two ever disagree, the
    on-screen sheet stops being a preview of the printed one."""
    css = (REPO / "concretedesignpy/webapp/static/css/calcsheet.css").read_text()
    assert "--sheet-w: 182mm" in css
    margin = re.search(r"@page\s*\{[^}]*margin:\s*(\d+)mm", css)
    assert margin, "no @page margin found"
    assert 210 - 2 * int(margin.group(1)) == 182
