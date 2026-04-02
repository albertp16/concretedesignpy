# ConcreteDesignPy

A Python library and web application for reinforced concrete structural design per **NSCP 2015** and **ACI 318-19**.

> **Disclaimer:** This is an internal tool for Albert Pamonag Engineering Consultancy (APEC). All results must be verified by a licensed professional engineer. The developers assume no liability for errors or misuse.

## Features

- **Beam Design** — Moment capacity (strain compatibility), shear strength, torsion design, deflection checks
- **Column Design** — P-M interaction diagrams, rebar layout generation, minimum flexural strength checks
- **Joint Shear** — Verification for special moment frames (Section 422.7)
- **Mander Model** — Confined concrete strength and ultimate strain
- **Moment-Curvature** — 6-point M-φ analysis
- **Development Length** — Hook geometry for deformed bars and stirrups (Section 425)
- **Alternative Inertia** — Effective moment of inertia per Section 424.2.3.5
- **Web Application** — Flask-based interactive calculators with inline SVG diagrams and MathJax LaTeX reports

## Installation

```bash
pip install -e .
```

### Dependencies

- Python 3.8+
- Flask
- NumPy

## Running the Web App

```bash
python -m concretedesignpy.webapp.app
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
│   ├── moment_curvature.py    # M-φ analysis
│   ├── development_length.py  # Hook geometry
│   ├── alternative_inertia.py
│   ├── rebar_layout.py        # Section rebar coordinates
│   └── diagrams.py            # SVG diagram generators
├── webapp/
│   ├── app.py            # Flask application
│   ├── routes/           # API endpoints (beam, column, section, joint, mander)
│   ├── templates/        # Jinja2 HTML templates with MathJax
│   └── static/           # CSS & JavaScript
```

## API Endpoints

| Endpoint | Description |
|---|---|
| `POST /api/beam/moment` | Beam moment capacity |
| `POST /api/beam/shear` | Shear strength or spacing |
| `POST /api/beam/torsion` | Torsion design |
| `POST /api/beam/deflection` | Deflection check |
| `POST /api/column/interaction` | P-M interaction diagram |
| `POST /api/column/flexural` | Min flexural strength |
| `POST /api/joint/shear` | Joint shear verification |
| `POST /api/mander/confined` | Mander confined concrete |
| `POST /api/section/development-length` | Hook geometry |
| `POST /api/section/moment-curvature` | M-φ analysis |
| `POST /api/section/alternative-inertia` | Alternative Ie |

## Code Standards

- NSCP 2015 (National Structural Code of the Philippines)
- ACI 318-19 (American Concrete Institute)
- ACI SP-17 (Design Handbook)
- Mander, Priestley & Park (1988) — Confined concrete model

## Contributing

1. Fork the repository
2. Create a feature branch (`git checkout -b feature-name`)
3. Commit your changes
4. Submit a pull request

## License

MIT License. See [LICENSE](LICENSE) for details.
