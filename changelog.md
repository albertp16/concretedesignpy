# Changelog - concretedesignpy

## Version 0.7.3 | August 15, 2026

### Changed
- **The client report now carries no advisory content at all** — not the text,
  not the codes, not the count (user direction, 2026-08-15, given after
  reviewing the v0.7.2 disclosure block on screen). The report runs Input
  Summary → Calculation Sheet → Detailing → QAQC → Results Summary →
  Provenance, and the word "advisory" does not appear in it.
- The disclosure summary (count, critical count, full code list, and the note
  that the exclusion was directed) moves to the head of the **Internal Review**
  panel, so the Engineer of Record still sees in one line what the calculation
  raised. It is still built from the response, never hand-written.
- The client verdict no longer references advisories, which would have been a
  reference to something the document does not contain. It is now a pure scope
  statement: the report covers the computed code checks only, is not a
  statement that the retrofit is safe or complete, and is subject to the
  review and judgement of the Engineer of Record. The verdict must not
  overclaim even when the caveats are filed elsewhere.

### Decision record
- v0.7.2 omitted the advisory text but disclosed the codes; v0.7.3 removes the
  disclosure too. The concern — that an omission which leaves no trace on the
  signed sheet is indistinguishable from an answer that never had the finding
  — was raised and the direction was reaffirmed. That is the Engineer of
  Record's call to make; this entry is the record of it.

---

## Version 0.7.2 | August 15, 2026

### Changed
- **Advisories separated from the client report** (user direction, 2026-08-15;
  the report is submitted to a client). The report no longer reproduces the
  advisory text. It carries a **disclosure block** instead — the count, the
  critical count, and every advisory code with its severity, built from the
  response and never hand-written — stating that the omission was directed.
  Both halves are pinned by tests: reverting the omission needs a new
  decision, and losing the disclosure would turn a directed omission into a
  silent one, which is indistinguishable from an answer that never had the
  finding.
- Advisories now render in a separate **Internal Review** panel, marked as
  excluded from the client submission and deliberately placed outside the
  report container so no selector can sweep it into the client print.
- **Two print modes** — *Print Client Report* (default) and *Print Internal
  Review*, keyed off a `body.print-internal` class. The internal handler
  clears the class in a `finally`, so a cancelled dialog cannot leave the page
  armed to print advisories into the next client report.
- The verdict's closing sentence is now doubled: on screen it points at the
  Internal Review panel; on the client copy it states that advisories are
  issued separately. A verdict that pointed at advisories the document does
  not contain would be a dangling reference on the sheet that gets signed.

---

## Version 0.7.1 | August 15, 2026

### Added
- **Calculation-sheet report for Column Jacketing** — the page report now
  follows the APEC three-column REFERENCES | CALCULATIONS | RESULTS format:
  a 31-row Input Summary (Symbol | Input | Description | Value, every field
  echoed so the report is self-contained), eight numbered calculation
  sections with symbolic and substituted equations, a Results Summary table,
  and the detailing/advisories/provenance blocks
- **Server-side QAQC block** — `column_jacket_design()` results now carry a
  `qaqc` section: 16 independent recomputations (9 for a partial jacket)
  that re-derive reported values from the raw inputs along separate
  arithmetic paths — geometry, bar counts and areas, Po and phi*Pn,max
  longhand, phi*Mn at Pu by re-interpolating the returned curve, the full
  shear stack, v_u summation, and Mander f'cc from the reported f_l. The
  page renders the table; it compares nothing itself. A deliberate-corruption
  test proves the block can actually fail
- **Professional figures** — the section drawing gains true dimension lines
  with witness ticks (B, H), a leader for t, centrelines and a labelled core;
  the P-M chart gains the demand annotation and drops the raw nominal
  envelope, which reached Po and pure tension and crushed the design curves
  into a corner of the axes

### Improved
- QAQC comparisons coerce numpy booleans so results stay JSON-serialisable

---

## Version 0.7.0 | August 15, 2026

### Added
- **RC Column Concrete Jacketing** — new modules `calculators/column_jacket.py`
  (engine) and `calculators/column_jacket_design.py` (design boundary), per
  TN-RET-001, ACI 318-19, ACI 562-16 and ASCE 41-17:
  - P-M interaction of the composite section using a **single** ACI stress
    block with β₁ from the extreme compression fibre. Giving each concrete its
    own β₁c is physically inverted and unconservative; it remains available as
    `stress_block='per_material'` for comparison only
  - Mander confinement with **k_e computed from the actual tie and bar
    geometry**, not assumed. A plain perimeter hoop gives k_e ≈ 0.40, not the
    0.80 assumed by habit, which overstates f'cc by roughly 12%
  - One-way shear with `d` taken from the **actual extreme tension bar**, so
    the shear check and the interaction diagram agree about where the steel is
  - Interface shear transfer reporting **both** ACI capacity routes (shear
    friction and the Table 16.4.4.2 bond route) and picking neither
  - Unshored-jacket preload split, stiffness feedback, detailing checks, axial
    strength, and induced-eccentricity geometry
  - Lampropoulos et al. (2012) monolithic coefficients **reported and never
    applied**, with both the corrected and uncorrected results side by side
  - Partial (one/two-opposite/three-sided) jackets compute geometry, axial
    strength and P-M interaction, and declare the checks that have no model
    instead of answering with four-sided behaviour. Two *adjacent* faces are
    refused outright — that section's capacity is a biaxial surface
  - Every result carries `provenance` (including a SHA-256 of the engine
    source) and an `advisories` list of traps the arithmetic cannot rule out
- **Column Jacketing calculator page** — `/calculator/column-jacket`, backed by
  `POST /api/column-jacket/design`. All 35 inputs, P-M interaction chart
  (existing vs jacketed with the demand point), a dimensioned section figure
  drawn from the engine's own bar coordinates, a MathJax calculation report,
  detailing table, and the full advisories panel. Ported from the APEC
  jacketing service's browser calculator, keeping its two load-bearing
  properties: the page holds no equations or engine constants — every number
  and both figures come from the server — and advisories render in full,
  never filtered, capped or truncated (the hover reveal only unpaints the
  remainder, and print reveals everything)
- **Test suite** — first `tests/` directory in the repo, 109 tests covering the
  jacketing engine, its design boundary, and the web layer. Independent
  recomputations are separated from labelled regression pins

### Changed
- `setup.py` version synced to `__init__.__version__` (it had been left at
  0.5.0 while the package reported 0.6.0)
- Navbar version badge was hardcoded `v0.3`; now `v0.7`
- README: retrofitting/strengthening features, library usage, testing section,
  and the ACI 562-16 / fib B14 / Lampropoulos references

---

## Version 0.4.0 | April 5, 2026

### Added
- **Advanced Moment-Curvature Analysis** — Hognestad parabolic concrete model with fiber-layer approach, axial load support, 60-point smooth curve
- **4-Panel Section Analysis** — Strain/stress charts with full section visualization
- **Strength Reduction Factor Plot** — Phi plot matching ACI 318-19 Table 21.2.2 with projected design values
- **Confinement Effectiveness Charts** — 19-ratio curves with projections, Fig. 5 style (Mander et al.)
- **Biaxial Confinement** — Full biaxial formulas in Mander report
- **Rebar Input Improvement** — Bar diameter and quantity inputs replace direct As entry, auto-compute area
- **Complete Mander Report** — 13-step transparent report with all intermediate calculations, paper references
- **Excel Export** — Stirrup spacing design with Excel export (beam shear)
- **Print Report** — Print-ready output with engineering notation

### Improved
- Moment-curvature page rewritten with Plotly charts, section visualization, QAQC report
- Mander report enhanced with paper references, placeholders, biaxial formulas, summary
- Confinement chart: inverted y-axis, fixed title/label overlap, 19 ratio curves
- Beam shear page rewritten with MathJax report, Plotly charts, clean layout
- RC section diagram: proper stirrup legs, force arrows, Vc vs Nu/Ag chart
- Split layout (inputs left, output right) for better UX
- Combined shear & torsion design with ACI 318M-14 references

### Fixed
- Phi plot: fixed x-axis range matching reference design
- Confinement chart: inverted y-axis, fixed title/label overlap
- RC section: top bars at all legs, black fill, proper rendering
- Vc vs Nu/Ag chart matching ACI 318M-14 Fig. R22.5.6.1 exactly
- d/c dimensions, Tu curved arrow, Tu input fixes

### Changed
- Replaced FLECHA attribution with Hognestad (1951) paper reference
- Updated README: removed Railway deploy section, improved documentation

---

## Version 0.3 | 2025

### Added
- Column P-M interaction diagram with Plotly
- Joint shear verification for special moment frames
- Development length hook geometry
- Alternative inertia calculator
- Mander confined concrete model
- Moment-curvature 6-point analysis
- Flask web application with Plotly charts and MathJax reports
- Railway/Heroku deployment support

---

## Version 0.2 | 2025

### Added
- Beam shear and torsion calculators
- Beam deflection calculator
- Column flexural strength ratio
- Web interface with route-based API

---

## Version 0.1.1 | March 3, 2025

### Added
- **RC Beam Solver** — Computation of Beam Capacity based on NSCP 2015 with validation checks
- **Midas RC Pushover MGT Generator** — Initial development (partial implementation)
