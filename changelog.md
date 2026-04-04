# Changelog - concretedesignpy

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
