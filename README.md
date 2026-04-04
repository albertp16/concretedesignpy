# ConcreteDesignPy

[![PyPI version](https://img.shields.io/pypi/v/concretedesignpy)](https://pypi.org/project/concretedesignpy/)
[![Python](https://img.shields.io/pypi/pyversions/concretedesignpy)](https://pypi.org/project/concretedesignpy/)
[![License: MIT](https://img.shields.io/badge/License-MIT-blue.svg)](https://opensource.org/licenses/MIT)
[![GitHub release](https://img.shields.io/github/v/release/albertp16/concretedesignpy)](https://github.com/albertp16/concretedesignpy/releases)
[![GitHub issues](https://img.shields.io/github/issues/albertp16/concretedesignpy)](https://github.com/albertp16/concretedesignpy/issues)

A Python library and web application for reinforced concrete structural design per **NSCP 2015** and **ACI 318-19**.

> **Disclaimer:** This is an internal tool for Albert Pamonag Engineering Consultancy (APEC). All results must be verified by a licensed professional engineer. The developers assume no liability for errors or misuse.

## Features

### Beam Design
- **Flexural Capacity** — Strain compatibility method with iterative neutral axis solver, strength reduction factor (phi) plot per ACI 318-19 Table 21.2.2
- **Shear Strength** — Concrete and steel shear capacity, stirrup spacing design with Excel export (ACI 318M-14)
- **Torsion Design** — Combined shear-torsion interaction checks
- **Deflection** — Immediate and long-term deflection with effective moment of inertia

### Column Design
- **P-M Interaction** — Strain compatibility analysis with Plotly interaction diagrams, load combination overlay
- **Minimum Flexural Strength** — Column-to-beam strength ratio verification

### Analysis Tools
- **Moment-Curvature** — Dual-mode M-phi analysis:
  - *Quick (6-Point)* — Closed-form key points (cracking, yield, ultimate) per ACI 318-19
  - *Advanced (Incremental)* — Hognestad parabolic concrete model with fiber-layer approach, axial load support, 60-point smooth curve
- **Mander Confined Concrete** — Full Mander, Priestley & Park (1988) model with 13-step transparent report, biaxial confinement (Fig. 5), Plotly stress-strain curves, and confinement effectiveness charts
- **Joint Shear** — Verification for special moment frames (Section 422.7)

### Detailing
- **Development Length** — Hook geometry for deformed bars and stirrups (Section 425)
- **Alternative Inertia** — Effective moment of inertia per Section 424.2.3.5

### Web Application
- Flask-based interactive calculators
- Plotly charts for all visualizations (M-phi diagrams, P-M interaction, stress-strain curves, section analysis)
- MathJax-rendered calculation reports with full formula substitution
- Print-ready output with engineering notation

## Installation

```bash
pip install -e .
```

### Dependencies

- Python 3.8+
- Flask
- NumPy
- Gunicorn (production)
- openpyxl

## Running Locally

```bash
python run.py
```

Then open [http://localhost:5000](http://localhost:5000) in your browser.

## Project Structure

```
concretedesignpy/
├── calculators/          # Backend calculation modules
│   ├── beam_moment.py    # Flexural capacity (neutral axis solver)
│   ├── beam_shear.py     # Concrete & steel shear strength
│   ├── beam_torsion.py   # Torsion design checks
│   ├── beam_deflection.py
│   ├── column_interaction.py  # P-M interaction diagram
│   ├── column_flexural.py     # Min flexural strength ratio
│   ├── joint_shear.py
│   ├── mander.py              # Confined concrete model
│   ├── moment_curvature.py    # M-phi analysis
│   ├── development_length.py  # Hook geometry
│   ├── alternative_inertia.py
│   ├── rebar_layout.py        # Section rebar coordinates
│   └── diagrams.py            # SVG diagram generators
├── webapp/
│   ├── app.py            # Flask application factory
│   ├── routes/           # API endpoints (beam, column, section, joint, mander)
│   ├── templates/        # Jinja2 HTML templates with MathJax
│   └── static/           # CSS & JavaScript
├── Procfile              # Railway/Heroku deployment
├── railway.json          # Railway config
├── requirements.txt      # Python dependencies
└── setup.py              # Package setup
```

## API Endpoints

See the [Wiki — API Reference](https://github.com/albertp16/concretedesignpy/wiki/API-Reference) for the full list of 17 POST endpoints.

## Code Standards & References

- **NSCP 2015** — National Structural Code of the Philippines
- **ACI 318-19 / ACI 318M-14** — American Concrete Institute
- **ACI SP-17** — Design Handbook
- **Mander, Priestley & Park (1988)** — Theoretical stress-strain model for confined concrete, *J. Structural Engineering*, ASCE, Vol. 114, No. 8
- **Hognestad (1951)** — A study of combined bending and axial load in reinforced concrete members, University of Illinois Bulletin No. 399

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature-name`)
3. Commit your changes
4. Submit a pull request

## License

MIT License. See [LICENSE](LICENSE) for details.
