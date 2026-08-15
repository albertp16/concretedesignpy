# concretedesignpy — rectification prompt (Claude Code)

> Paste everything below the rule into Claude Code from the repo root
> `/Users/albertpamonag/Desktop/Projects/concretedesignpy`.
> Companion to `vault/_Clause Verification Prompt.md`, which is the *audit* prompt.
> This one is the *rectification* prompt: the audit is done, the findings are in
> `Software/14`, and the job is to land them without inventing anything new.

---

You are working in `concretedesignpy` (Python, Flask + calculators) at
`/Users/albertpamonag/Desktop/Projects/concretedesignpy`.

The **Structural Mind vault** at
`/Users/albertpamonag/Desktop/Programs/structural-mind/Structural Mind/` is the
source of truth for every technical claim you make in this task.

## 0 — Read these first, in this order, before touching code

1. `<vault>/CLAUDE.md` — the vault's own operating rules. **They bind you here too.**
2. `<vault>/vault/Software/14 - concretedesignpy Review — Beam Flexure, Shear, Torsion & Joint Shear.md`
   — the audit you are implementing. §8 is the rectification order. **Do not re-derive
   its findings; implement them.**
3. `<vault>/vault/Software/concretedesignpy_verify.py` — the harness that produced every
   number in that note. **Run it before you change anything and capture the baseline.**
4. `<vault>/vault/reference/README.md` — the page-offset register. You need it to quote a page.
5. Skim `<vault>/vault/Software/12 - rc-beam Solver Deep Review…` and
   `<vault>/vault/Software/13 - PBDPy Clause Audit…` — the same three failure modes
   (silent zeros, edition drift, duplicated implementations) recur across all three repos.

## 1 — Hard rules

**1.1 Two sources, two jobs — never mix them.**

| Source | For | NOT for |
|---|---|---|
| `<vault>/vault/reference/` | **All technical basis** — provisions, equations, tables, coefficients, acceptance criteria | project facts |
| the repo, its git history, its docstrings, its web form, Slack | **What the code currently does and why** | **any technical basis whatsoever** |

A comment in `beam_shear.py` saying `# ACI 318M-14 Cl. 22.7.6.1a` is **a claim to be
checked, not a citation.** Several of them are wrong — that is the whole point of the audit.

**1.2 Cite from `reference/` only.** Every provision you invoke carries:

> **Basis:** `[file]`, Section/Clause `[n]`, p. `[printed page]`

using the **printed** page (apply the `reference/README.md` offset —
`ACI 318-25M` printed = PDF − 1; `Wight MacGregor 7th` printed = PDF − 1; NSCP Ch. 4
folios are chapter-relative, no linear offset).

**1.3 If `reference/` cannot support it, say exactly:**
> The reference folder does not contain enough basis to support this answer.

Then **stop and ask**. Never invent a clause, a page, a coefficient or a limit.

**1.4 Image-only files may not be quoted from a text search.** NSCP 2015 Vol. I, ACI 369.1-22,
ACI 562-25, AISC 358 and everything under `_Scanned/` have no text layer — render the page
at 300–400 dpi and *read* it, or mark the page unverified.

**1.5 Separate code requirements from software guidance.** A Flask route, a Jinja template,
a function signature and "delete this module" are **implementation decisions**. Mark them as
such. Only the provisions cited from `reference/` are code.

**1.6 Never widen scope.** These modules are **out of scope and must not be touched**:
`column_jacket*`, `frp_*`, `column_biaxial/flexural/interaction`, `moment_curvature`,
`development_length`, `beam_deflection`, `mander`, `rebar_layout`, `mgt/`, `stress-strain/`.
If a change to an in-scope module would break one of these, stop and report it.

## 2 — Governing edition — settle this before writing a line

Every docstring in the package claims **ACI 318-19**. Every implementation is
**ACI 318-14 / NSCP 2015**. That is the right code for Philippine work and the wrong label.

**Decide once, explicitly, and record the decision in `CHANGELOG` and in the note:**

- **Option A (recommended, default):** the package targets **NSCP 2015 (≡ ACI 318-14 numbering)**.
  Fix every docstring to say so. Then the pre-2019 φ, the 318-14 `Vc`, and the 3-row joint
  table are all **correct as written** and only the genuine defects below get fixed.
- **Option B:** the package targets **ACI 318-19/-25**. Then φ, `Vc` (size-effect `λs`) and the
  joint table all have to move as well, and every existing result changes.

**Do not silently mix them.** If you are unsure, take Option A and flag Option B as a
follow-up — that is the conservative reading of the governing PH code.

## 3 — The work, in dependency order

Land these as **separate commits**, each with its own test. Do not batch.

### P0 — Delete the duplicate torsion implementation ⛔⛔ *blocking*

`concretedesignpy/calculators/beam_torsion.py` is mis-scaled by 1000× and returns
`"MAY NEGLECT"` and `"status": "PASS"` for a beam at 4.8× the threshold.

- `beam_torsion.py:82` — `phi_tth` divides by `10.0 * 100.0` where N·mm → kN·m needs `1e6`.
- `beam_torsion.py:119` — `at_per_s` carries the same error.
- `beam_torsion.py:103–109` — the §22.7.7.1 crushing terms are mis-scaled **by two different
  factors** (shear ~100× low, torsion ~10⁵× low), so `combined_stress` cannot reach the limit.
- `templates/calculators/beam-torsion.html:37` labels the field `Tu (kN-m)`, which is what
  makes the comparison unit-inconsistent on every call.

**Do this:** delete `torsion_design()`, re-point `POST /api/beam/torsion`
(`webapp/routes/beam.py:399–421`) at `beam_shear.shear_torsion_design()`, and adapt the
template. **Do not patch the five scale factors individually** — they are individually
plausible and that is exactly how they survived.

⚠ Also fix the template's dead keys: it reads `result.t_threshold` and `result.torsion_status`;
the module returns `torsion_threshold` and `status`. The UI has been showing nothing.

**Acceptance:** W&M Ex 7-2 (printed 358–362), converted to SI, must return
`φT_th ≈ 7.96 kN·m`, `At/s ≈ 0.5183 mm²/mm`, and *"design for torsion"* at `Tu = 37.96 kN·m`.

### P1a — `Aℓ` is 15.00 % short ⛔ *unconservative, one line*

`beam_shear.py:370` — `al = tu * 1e6 * ph / (2 * phi * aoh * fy)`

**ACI 318-25M Eq. (22.7.6.1b), printed 465** gives `Tn = (2 Ao Aℓ fy / ph) tan θ` with
**`Ao = 0.85 Aoh`** (§22.7.6.1.1) and θ = 45°, so the divisor is **`1.7 Aoh`**, not `2 Aoh`.
Identical to **W&M Eq. (7-36), printed 349**. Ratio is exactly 1.7/2 = 0.85 for every section.

**Acceptance:** W&M Ex 7-2 → `Aℓ = 816 mm²` (currently 694). Add a regression pin at that value.

### P1b — the joint has no column-width cap on `Aj` ⛔ *unconservative*

`joint_shear.py:71–74` computes `bj = min(beam_width + joint_depth, beam_width + 2·perpendicular_dist)`
and stops. `joint_shear_check` takes **no column-width argument at all**.

**ACI 318-25M §15.5.2.2, printed 231–232**, verbatim: *"Effective joint width shall be the
overall width of the column where the beam is wider than the column. Where the column is
wider than the beam, effective joint width shall not exceed the lesser of (a) and (b)…"* and
**R15.5.2.2, printed 231:** *"In no case is `Aj` greater than the column cross-sectional area."*

**Do this:** add a required `column_width` parameter; `bj = min(column_width, b + h, b + 2x)`.
Keep `perpendicular_dist` measured **beam face → column face** (its current docstring
definition) — that form is algebraically identical to ACI's `2x` measured from the beam axis,
and the existing behaviour for a beam narrower than the column is correct. Say so in the
docstring so the next reader does not "fix" it back.

**Acceptance:** a 500 mm beam on a 400 × 600 column must give `Aj = 240 000 mm²`, not 300 000.
W&M Ex 19-3 (printed 1076–1078) must still give `Aj = 371 612 mm²`, `Vj = 1838 kN`.

### P2a — shear: two guards that exist in siblings but not in `compute_shear_spacing()`

`beam_shear.py:130–181`. When `Vu > φVc` (`:160–161`) the spacing comes from `Vs` alone and
`Av,min` is never applied; and there is no `Vs ≤ ⅔√f'c bw d` test at all — `shear_design()`
and `shear_torsion_design()` both have one.

- `Av,min` — **ACI 318-25M Table 9.6.3.4, printed 151** (`0.062√f'c bw/fyt` and `0.35 bw/fyt`).
- `Vn ≤ Vc + 0.66√f'c bw d` — **ACI 318-25M §22.5.1.2, printed 442**.

**Acceptance:** a wide lightly-loaded web must return the `Av,min` spacing, not the `Vs` one;
an over-stressed section must raise or return an explicit `UNSAFE`, never a number.

### P2b — joint: express the real table

`joint_shear.py:77` — `{1: 1.7, 2: 1.2, 3: 1.0}`.

Under **Option A** this **is** NSCP 2015 Table 418.8.4.1 (printed 4-118) and is correct —
leave the values, but rename the enum so it says NSCP, and document that it is a 3-row table.

Under **Option B**, **ACI 318-25M Table 18.8.4.3, printed 342** is an **eight-row** table keyed
on *column continuity (§15.5.2.3) × beam continuity (§15.5.2.4) × transverse-beam confinement
(§15.5.2.5)*: **1.7 / 1.3 / 1.3 / 1.0 / 1.3 / 1.0 / 1.0 / 0.7**. ⚠ That table is legible in the
text layer but **read the rendered page before transcribing it** — the extraction interleaves rows.

Either way: constrain `lamda` to {0.75, 1.0} (Table 18.8.4.3 footnote [1]) and stop letting a
caller override `phi` — **§21.2.4.4, printed 435** fixes it at 0.85 for special-moment-frame joints.

### P3 — the small ones

| Item | Where | Basis | Printed |
|---|---|---|---|
| φ law — import the strain-based one already in `provision.py` *(Option B only)* | `beam_moment.py:53–61` | ACI 318-25M Table 21.2.2 | 430 |
| `Cs = A′s(f′s − 0.85 f′c)` — deduct displaced concrete | `beam_moment.py:206–214` | W&M Ex 4-4 step 5 | 170 |
| `Aoh`/`ph` on the **stirrup** diameter, not the main bar | `beam_torsion.py:66–69` | ACI 318-25M §22.7.6.1 commentary | 465 |
| torsional stirrup `s_max` → **300 mm**, not 305 | `beam_torsion.py:131` | ACI 318-25M §9.7.6.3.3 | 160 |
| one `Av,min` constant — 0.062, not `1/16` in two places | `beam_shear.py:163` vs `:388`, `:530` | ACI 318-25M Table 9.6.3.4 | 151 |
| the §22.5.6.1 axial cap expression | `beam_shear.py:250`, `:502` | ACI 318-25M §22.5.6.1 | — |
| `Nu` term of the torsion thresholds | `beam_shear.py:326–328` | ACI 318-25M Table 22.7.4.1(a), Table 22.7.5.1 | 463, 464 |
| `shear_status` computed then discarded | `beam_shear.py:274–282` | — *(software)* | — |

⚠ **Two of these are cosmetic and you must say so rather than dress them up:** the displaced-concrete
term moves `Mn` by a **median +0.05 %, max +1.05 %** over 4 000 sections, and the §22.5.6.1 cap
**never governs** for any `Nu ≥ 0` because `min()` always selects the other branch. Fix them for
correctness; do not report them as safety findings.

### P4 — docstrings

Every module header claims ACI 318-19. Make them match §2's decision. Also state plainly what
each function does **not** check — `calculate_beam_moment` computes capacity and verifies
neither `As,min` (§9.6.1.2, printed 149), nor `ρ ≤ 0.025`, nor §18.6.3.1/.2.

## 4 — Tests are part of the work, not a follow-up

`tests/` currently holds **three files, all `column_jacket`**. There is **no test anywhere** on
`beam_moment`, `beam_shear`, `beam_torsion` or `joint_shear` — which is why a 1000× error shipped.

Add `tests/test_beam_flexure.py`, `tests/test_beam_shear.py`, `tests/test_beam_torsion.py`,
`tests/test_joint_shear.py` containing:

1. **Textbook pins** — one test per worked example, asserting the printed value within 1 %,
   with the printed page in the test's docstring:

   | Topic | Example | Printed | Pin |
   |---|---|---|---|
   | flexure, singly, SI | W&M 4-1M | 148 | `Mn = 273 kN·m`, `a = 151`, `c = 178` |
   | flexure, **doubly** | W&M 4-4 | 170–171 | `Mn = 428 k-ft`, `c = 6.2 in.` |
   | shear stirrup design, SI | W&M 6-1M | 297–300 | `Vc = 153 kN`, `s = 215 mm`, `Vs = 86.6 kN` |
   | torsion + shear, solid | W&M 7-2 | 358–362 | `φT_th = 5.87 k-ft`, `At/s = 0.0204 in²/in`, `Aℓ = 1.27 in²` |
   | torsion + shear, **hollow** | W&M 7-3 | 363–367 | `Aoh = 338 in²`, `ph = 74 in`, `Aℓ,min governs at 1.57 in²` |
   | **interior joint** | W&M 19-3 | 1076–1078 | `Vj = 413 kips`, `Aj = 576 in²`, `φVn = 619 kips` |

2. **Regression pins on the four fixes** — the exact numbers in the Acceptance lines above.
3. **Unit-sanity tests** — the class of bug P0 is. For each public function, feed a section
   whose answer you can hand-compute to one significant figure and assert the *order of
   magnitude*. A 1000× error must fail loudly.
4. **No-silent-acceptance tests** — assert that a section with zero capacity, zero stirrups or
   an infeasible layout **raises or returns an explicit failure**, never a passing verdict.
   `beam_moment` already does this correctly; use it as the model.

Run `pytest` and report the count. **Do not mark an item done with a failing or skipped test.**

## 5 — Adversarial pass — mandatory before you report

Before writing the summary, go back over your own work and try to break it:

1. **Re-read every printed page you cited**, in the file, at the page. Not from memory, not from
   this prompt. If a page number here is wrong, say so — this document is not authority either.
2. **Quantify every finding before you rank it.** A change that "looks unconservative" gets a
   number attached — a sweep, a worked case — or it is reported as unquantified.
3. **Check whether a wrong expression ever governs.** Two findings in the original audit were
   downgraded this way and one was retracted outright.
4. **Check the parameter definition before you accuse the caller.** The retracted finding came
   from feeding a function a quantity its own docstring defined differently.
5. **Diff the behaviour, not just the source.** Run the harness before and after and put the
   before/after table in the report.

## 6 — Output

1. **Answer** — what changed, per commit.
2. **Basis & Citations** — table: Point | Reference file | Clause/Section | Printed page | Relevance.
3. **Verification** — "Verified from reference folder" / "Partially verified…" / "Not verified…",
   per item, plus the pytest result and the before/after harness table.
4. **Summary** — what is closed, what is still open, what you refused to do and why.

Then write back to the vault:

- append a dated entry to `<vault>/vault/_Knowledge Log.md` (newest at top, matching the
  existing `<!-- log-YYYY-MM-DD-slug -->` format);
- amend `<vault>/vault/Software/14 - concretedesignpy Review…` §8 to mark each item closed,
  **without rewriting the findings** — an audit note keeps its history;
- keep `<vault>/vault/Software/00 - Software Index.md` and `<vault>/vault/00 - Index.md` in sync;
- mirror the fix list to `concretedesignpy/CLAUSES.md` in the repo, the same way PBDPy does.

## 7 — Standing traps in this vault

- **`ACI 318-19` on file is INCH-POUND.** SI coefficients are in **Appendix C, printed 589**.
  Use **`ACI 318-25M`** (printed = PDF − 1) for anything numeric.
- **W&M 7th is written to ACI 318-14**, so it **cannot** catch an edition finding — it prints
  the old φ rule on p. 169. It is a benchmark for arithmetic, not for provisions.
- **W&M Ch. 6/7/19 examples are inch-pound** except 6-1M. Convert, and say you converted.
- **W&M rounds where the code does not** — Eq. (6-21M) prints `(1/3) bw/fyt` where ACI 318M
  prints **0.35**; the joint γ is **20/15/12 psi** = 1.67/1.25/1.00 SI against NSCP's **1.7/1.2/1.0**.
  A 2 % gap against the textbook is not automatically a bug.
- **Albert wants backend modules, tests and a CLI — not HTML calculators.** Do not add UI
  beyond re-pointing the route P0 requires.
- **The vault may be edited by another session.** Re-read a vault file immediately before you
  write it, and merge additively. Never force-overwrite.
