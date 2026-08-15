"""Web layer tests for the column jacketing calculator.

These check the plumbing between the Flask app and the calculation module, and
the two page properties that are not cosmetic: the page must post to the
server for every number it shows, and it must not be able to hide an advisory.
"""

import copy
import json
import re
from pathlib import Path

import pytest

from concretedesignpy.webapp.app import create_app

TEMPLATE = (Path(__file__).resolve().parent.parent / "concretedesignpy" /
            "webapp" / "templates" / "calculators" / "column-jacket.html")


def reference_payload(**overrides):
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
def client():
    return create_app().test_client()


@pytest.fixture(scope="module")
def page_source():
    return TEMPLATE.read_text()


# ------------------------------------------------------------------ routing


def test_the_calculator_is_listed_on_the_landing_page(client):
    """The home grid is a hard-coded list -- a page not in it is unreachable."""
    body = client.get("/").get_data(as_text=True)
    assert "Column Concrete Jacketing" in body
    assert "/calculator/column-jacket" in body


def test_the_page_renders(client):
    r = client.get("/calculator/column-jacket")
    assert r.status_code == 200
    body = r.get_data(as_text=True)
    assert "RC Column Concrete Jacketing" in body
    assert 'id="calc-form"' in body


def test_design_endpoint_returns_the_success_envelope(client):
    r = client.post("/api/column-jacket/design", json=reference_payload())
    assert r.status_code == 200
    body = r.get_json()
    assert body["status"] == "success"
    res = body["result"]
    assert res["interaction"]["phiMn_at_Pu_jacketed"] == pytest.approx(535.25, rel=1e-3)
    assert res["axial"]["Po_jacketed"] == pytest.approx(9737.1, rel=1e-3)
    assert res["adequate"] is True


def test_response_is_json_serialisable_over_the_wire(client):
    """A NaN would serialise to invalid JSON and break the browser parse."""
    raw = client.post("/api/column-jacket/design",
                      json=reference_payload()).get_data(as_text=True)
    assert "NaN" not in raw and "Infinity" not in raw
    json.loads(raw)


def test_construction_block_is_optional(client):
    payload = reference_payload()
    del payload["construction"]
    r = client.post("/api/column-jacket/design", json=payload)
    assert r.status_code == 200
    echo = r.get_json()["result"]["request_echo"]["construction"]
    assert echo["continuity"] == "discontinuous"


# --------------------------------------------------------------- refusals


@pytest.mark.parametrize("patch,fragment", [
    ({"existing": {"cover_to_bar_centre": 160}}, "lever arm"),
    ({"jacket": {"thickness": 55, "bar_dia": 25}}, "cannot house"),
    ({"jacket": {"sides": "two_adjacent"}}, "biaxial"),
    ({"construction": {"P0_at_casting": 9999}}, "P0_at_casting exceeds Pu"),
    ({"existing": {"fc": 500}}, "at most"),
])
def test_a_refused_input_returns_400_with_the_exact_reason(client, patch, fragment):
    """The validation layer's message names the field and why the bound exists.
    Collapsing it to "invalid input" at the HTTP boundary would throw away the
    one thing that makes the refusal actionable."""
    r = client.post("/api/column-jacket/design", json=reference_payload(**patch))
    assert r.status_code == 400
    body = r.get_json()
    assert body["status"] == "error"
    assert fragment in body["message"]


def test_a_missing_input_group_is_a_400_not_a_500(client):
    payload = reference_payload()
    del payload["jacket"]
    r = client.post("/api/column-jacket/design", json=payload)
    assert r.status_code == 400
    assert r.get_json()["status"] == "error"


def test_a_partial_jacket_returns_200_and_declares_the_gaps(client):
    """Unavailable is not an error -- the geometry and P-M are real."""
    r = client.post("/api/column-jacket/design",
                    json=reference_payload(jacket={"sides": "three"}))
    assert r.status_code == 200
    res = r.get_json()["result"]
    assert res["adequate"] is False
    assert res["unavailable"]
    assert res["shear"] is None


# ------------------------------------------------- page is a thin client


def test_the_page_posts_to_the_design_endpoint(page_source):
    assert "/api/column-jacket/design" in page_source


def test_the_page_carries_no_engine_constants(page_source):
    """Every number shown must come from the server. A constant computed in the
    browser is one no test of the calculation layer can see, so it can drift
    from the engine indefinitely without anything failing.

    A code coefficient may be QUOTED inside a displayed formula -- the report
    prints `P_o = 0.85 f'_c (A_g - A_st) + f_y A_st` as documentation of what
    the server did. So the LaTeX literals are stripped first, and the check
    then applies to the executable remainder: a constant surviving that strip
    is one the page is using, not showing.
    """
    script = page_source.split("<script>")[-1]
    executable = re.sub(r"\\\\\[.*?\\\\\]", "", script, flags=re.S)
    # ...and rgba()/rgb() alphas are colours, not coefficients.
    executable = re.sub(r"rgba?\([^)]*\)", "", executable)
    assert "0.85" in script and "0.85" not in executable, \
        "0.85 should appear only inside a displayed formula or a colour"
    for forbidden in ("0.85", "0.003", "beta1 =", "phi =", "0.65", "0.75",
                      "Math.sqrt", "Math.pow", "eval(", "Function("):
        assert forbidden not in executable, \
            "the page must not compute {!r} itself".format(forbidden)


def test_the_page_does_no_arithmetic_on_returned_values(page_source):
    """Layout arithmetic for the two figures is allowed and is confined to the
    draw functions; anything else would be the page second-guessing the
    engine."""
    script = page_source.split("<script>")[-1]
    allowed = ["g.B / 2", "g.tie_cage_b / 2", "(g.H - g.tie_cage_h) / 2",
               "(g.H + g.tie_cage_h) / 2", "b.dia / 2", "-halfB - PAD",
               "halfB + PAD", "g.H + PAD", "b.z - r", "b.z + r",
               "b.y - r", "b.y + r"]
    for a in allowed:
        script = script.replace(a, "")
    for forbidden in (" * 1000", " / 1000", " * 0.85", " / 100 ", " * 1e"):
        assert forbidden not in script, \
            "unexpected arithmetic {!r} on a server value".format(forbidden)


def test_plotly_strings_use_unicode_not_html_entities(page_source):
    """Plotly does not decode named HTML entities, so `&phi;` renders as the
    literal seven characters in a chart title, axis label or legend. It DOES
    honour <sub>/<sup>. Caught by looking at the page, not by any assertion --
    hence this one.

    Entities remain correct in the figure captions, which are written with
    innerHTML -- so only the region that builds the traces and layout, up to
    the Plotly.newPlot call, is scanned.
    """
    for fn in ("function drawPM", "function drawSection"):
        plot_args = page_source.split(fn)[1].split("Plotly.newPlot")[0]
        for entity in ("&phi;", "&mdash;", "&times;", "&sub;", "&ldquo;"):
            assert entity not in plot_args, \
                "{} passes {!r} to Plotly; use the Unicode character".format(fn, entity)


def test_advisories_cannot_be_hidden(page_source):
    """The panel renders one entry per advisory, with no filtering, capping or
    severity gate. `adequate` covers the computed checks only, and a page that
    showed the DCRs alone would be quietly downgrading the answer."""
    panel = page_source.split("Advisories (")[1].split("</div>")[0]
    assert "d.advisories.map(advisory)" in page_source
    for forbidden in (".filter(", ".slice(", "severity ===", "> 3", "[0]"):
        assert forbidden not in panel, \
            "the advisories panel must not {!r}".format(forbidden)


def test_the_hover_reveal_hides_nothing(page_source):
    """Severity, code and the opening sentence are always painted; the rest is
    display:none until hover or focus. The full message stays in the DOM."""
    assert "a.message.slice(0, cut + 1)" in page_source
    assert "a.message.slice(cut + 1)" in page_source
    assert 'class="cj-rest"' in page_source
    # Hidden by `display` only -- never sliced out of the DOM -- and revealed
    # on hover AND focus, so a keyboard user is not cut off from it.
    rest_rule = re.search(r"\.cj-adv \.cj-rest \{([^}]*)\}", page_source).group(1)
    assert "display: none" in rest_rule
    assert re.search(r"\.cj-adv:hover \.cj-rest,\s*\.cj-adv:focus-within \.cj-rest",
                     page_source)


def test_print_does_not_drop_the_advisories(page_source):
    """Screen hides the remainder behind hover only because hover exists there.
    Paper has no hover, so the print rules must reveal every advisory in full
    rather than carrying a disclosure that something was omitted."""
    print_block = page_source.split("@media print {")[1]
    reveal = re.search(r"\.cj-adv \.cj-rest \{([^}]*)\}", print_block).group(1)
    assert "display: block !important" in reveal
    assert ".cj-adv" not in print_block.split("display: none !important")[0] \
        .split(".cj-adv .cj-rest")[0].replace(
            ".report-section, .cj-adv, .results-panel", ""), \
        "no print rule may hide an advisory"


def test_every_number_field_accepts_decimals(page_source):
    """input[type=number] defaults to step=1, which silently rejects 21.5 MPa.
    Integer counts are the only fields allowed to omit step."""
    for m in re.finditer(r'<input type="number"[^>]*>', page_source):
        tag = m.group(0)
        if 'data-int="1"' in tag:
            assert 'step=' not in tag, "integer field should not set step: " + tag
        else:
            assert 'step="any"' in tag, "real-valued field needs step=any: " + tag


def test_the_form_covers_every_input_the_module_accepts(page_source):
    """A field the module accepts but the page omits silently takes a default
    the engineer never saw."""
    # NB: the package re-exports a FUNCTION named column_jacket_design, which
    # shadows the module of the same name on attribute access -- so both
    # `from ... import column_jacket_design` and `import ....column_jacket_design
    # as m` hand back the function. import_module goes to sys.modules and gets
    # the module itself. (`from ....column_jacket_design import X` is
    # unaffected, which is why the public API still works.)
    import importlib

    m = importlib.import_module(
        "concretedesignpy.calculators.column_jacket_design")

    specs = {"existing": m._EXISTING_SPEC, "jacket": m._JACKET_SPEC,
             "demand": m._DEMAND_SPEC, "construction": m._CONSTRUCTION_SPEC}
    for group, spec in specs.items():
        for field in spec:
            needle = 'data-group="{}" data-key="{}"'.format(group, field)
            assert needle in page_source, \
                "form is missing {}.{}".format(group, field)
