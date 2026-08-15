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


def test_response_carries_a_passing_qaqc_block(client):
    """The report's QAQC section is server-computed. Its independent
    recomputations must agree with the report on the reference case."""
    r = client.post("/api/column-jacket/design", json=reference_payload())
    q = r.get_json()["result"]["qaqc"]
    assert q["all_pass"] is True
    assert q["n_pass"] == len(q["checks"]) >= 14
    for c in q["checks"]:
        assert c["name"] and c["method"], "a QAQC row without its method is not auditable"


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
    """Layout arithmetic for the two figures is allowed and is FENCED between
    the [figure-layout] markers; outside that region the page may format and
    compare server values but never derive one from another. Arithmetic in the
    report builder would be the page second-guessing the engine."""
    script = page_source.split("<script>")[-1]
    assert "// [figure-layout:start]" in script, "the fence markers must exist"
    assert "// [figure-layout:end]" in script
    pre, rest = script.split("// [figure-layout:start]", 1)
    _, post = rest.split("// [figure-layout:end]", 1)
    remainder = pre + post
    # displayed formulae are annotation, not computation
    remainder = re.sub(r"\\\\\[.*?\\\\\]", "", remainder, flags=re.S)
    for forbidden in (" * 1000", " / 1000", " * 0.85", " / 100 ", " * 1e",
                      " / 2", " * 2", "Math.sqrt", "Math.pow"):
        assert forbidden not in remainder, (
            "arithmetic {!r} on a server value outside the "
            "figure-layout fence".format(forbidden))


def test_qaqc_table_renders_server_computed_pairs_only(page_source):
    """Every expected/computed pair in the QAQC table comes from the server's
    qaqc block. The page compares nothing itself -- c.pass is rendered, never
    derived, so the QAQC verdicts stay inside the tested calculation layer."""
    assert "q.checks.forEach" in page_source
    body = page_source.split("function buildQaqc")[1].split("function buildSummary")[0]
    assert "c.expected" in body and "c.computed" in body and "c.pass" in body
    for forbidden in ("Math.abs", "toFixed(") :
        # formatting via num() is fine; comparison or arithmetic is not
        assert forbidden not in body.replace("num(c.expected", "").replace(
            "num(c.computed", ""), \
            "buildQaqc must not process the pair itself: {}".format(forbidden)
    assert "qq-banner-fail" in body, \
        "a failed QAQC run must produce a loud banner, not a quiet row"


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
    """The Internal Review panel renders one entry per advisory, with no
    filtering, capping or severity gate. They are separated from the client
    report by direction -- separation is not permission to lose any."""
    body = page_source.split("function buildInternalReview")[1] \
                      .split("function buildReport")[0]
    assert "d.advisories.map(advisory)" in body
    for forbidden in (".filter(", ".slice(", "severity ===", "> 3", "[0]"):
        assert forbidden not in body, \
            "the internal review panel must not {!r}".format(forbidden)


# ----------------------------------- client report / internal review split


def test_the_client_report_carries_no_advisory_content(page_source):
    """User direction, 2026-08-15 (revised after review of the disclosure
    block on screen): the submitted report carries NO advisory content at all
    -- not the text, not the codes, not the count.

    buildReport must therefore reference neither the advisory renderer nor the
    disclosure summary; both belong to buildInternalReview.
    """
    report = page_source.split("function buildReport")[1] \
                        .split("// ------------------------------------------------------------------ submit")[0]
    assert "d.advisories" not in report, \
        "the client report must not touch d.advisories in any form"
    assert "buildDisclosure" not in report, \
        "the disclosure summary belongs to the internal review, not the report"
    assert "advisory(" not in report


def test_the_disclosure_summary_lives_in_the_internal_review(page_source):
    """The count and the full code list are not lost -- they head the Internal
    Review panel, so the Engineer of Record still sees in one line how many
    findings the calculation raised and how many are critical.

    Built from the response, never hand-written: a hand-maintained list would
    drift the moment an advisory is added.
    """
    ir = page_source.split("function buildInternalReview")[1] \
                    .split("function buildReport")[0]
    assert "buildDisclosure(d)" in ir

    body = page_source.split("function buildDisclosure")[1] \
                      .split("// ------------------------------------------------- internal review")[0]
    assert "d.advisories.map(" in body, "the code list must be built from the response"
    assert "a.code" in body and "a.severity" in body
    assert "d.advisories.length" in body, "the count must be reported"
    for forbidden in (".filter(", ".slice("):
        assert forbidden not in body, \
            "the summary must cover every advisory, not a subset ({})".format(forbidden)
    # the exclusion is recorded as DIRECTED, so the decision stays visible to
    # whoever reads the internal copy
    assert "at the direction of the" in body and "Engineer of Record" in body


def test_the_internal_review_is_outside_the_client_report(page_source):
    """The client report is printed by hiding the internal review, so the
    review must not be a descendant of the document being submitted."""
    report_div = page_source.index('id="results-report"')
    review_div = page_source.index('id="internal-review"')
    assert review_div > report_div
    # the report container is self-closing before the review begins
    between = page_source[report_div:review_div]
    assert between.count("<div") - between.count("</div>") <= 1, \
        "internal-review appears to be nested inside results-report"


def test_print_keys_the_two_documents_off_a_body_class(page_source):
    """Default print is the CLIENT report; the internal review prints only
    under body.print-internal. Anything else risks advisories reaching a
    client submission through a stray selector."""
    print_block = page_source.split("@media print {")[1]
    assert ".internal-review { display: none; }" in print_block, \
        "the client print must hide the internal review"
    assert "body.print-internal .internal-review { display: block !important; }" in print_block
    # the advisory remainder is revealed on paper only in the internal document
    assert "body.print-internal .cj-adv .cj-rest { display: block !important; }" in print_block


def test_the_client_verdict_does_not_point_at_missing_advisories(page_source):
    """The screen verdict says "read every advisory below"; the client copy
    cannot, because they are not there. Both wordings must exist and the print
    rules must swap them."""
    verdict = page_source.split("function renderVerdict")[1].split("function drawPM")[0]
    assert 'class="screen-only"' in verdict and 'class="client-only"' in verdict
    # the client wording must not overclaim, and must not reference advisories
    # the document does not contain
    assert "computed code checks only" in verdict
    assert "advisor" not in verdict.split('class="client-only"')[1].split("</span>")[0].lower()
    print_block = page_source.split("@media print {")[1]
    assert ".screen-only { display: none !important; }" in print_block
    assert ".client-only { display: inline !important; }" in print_block


def test_the_internal_print_mode_cannot_leak_into_the_next_client_report(page_source):
    """A cancelled print dialog must not leave the page armed to print
    advisories into the following client report."""
    handler = page_source.split("btn-print-internal').addEventListener")[1]
    assert "finally" in handler and "classList.remove('print-internal')" in handler
    client = page_source.split("getElementById('btn-print').addEventListener")[1] \
                        .split("btn-print-internal")[0]
    assert "classList.remove('print-internal')" in client, \
        "the client print must clear the internal flag before printing"


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
