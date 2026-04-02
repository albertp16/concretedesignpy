# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Main page routes."""

from flask import Blueprint, render_template

main_bp = Blueprint("main", __name__)


@main_bp.route("/")
def index():
    """Landing page with calculator directory."""
    calculators = [
        {"id": "beam-moment", "name": "Beam Moment Capacity",
         "desc": "Flexural capacity with neutral axis solver",
         "category": "Beam"},
        {"id": "beam-shear", "name": "Beam Shear & Torsion",
         "desc": "Combined shear & torsion design per ACI 318M-14",
         "category": "Beam"},
        {"id": "beam-torsion", "name": "Beam Torsion Design",
         "desc": "Torsion checks per ACI/NSCP",
         "category": "Beam"},
        {"id": "beam-deflection", "name": "Beam Deflection",
         "desc": "Short-term and long-term deflection",
         "category": "Beam"},
        {"id": "column-interaction", "name": "Column Interaction & Rebar Layout",
         "desc": "P-M interaction diagram with section layout",
         "category": "Column"},
        {"id": "column-flexural", "name": "Column Flexural Check",
         "desc": "Minimum flexural strength ratio",
         "category": "Column"},
        {"id": "joint-shear", "name": "Joint Shear",
         "desc": "Joint shear for special moment frames",
         "category": "Joint"},
        {"id": "mander", "name": "Mander Confined Concrete",
         "desc": "Confined concrete strength and strain",
         "category": "Material"},
        {"id": "development-length", "name": "Development Length",
         "desc": "Hook geometry per NSCP 2015 Section 425",
         "category": "Detailing"},
        {"id": "moment-curvature", "name": "Moment-Curvature",
         "desc": "6-point moment-curvature relationship",
         "category": "Analysis"},
        {"id": "alternative-inertia", "name": "Alternative Inertia",
         "desc": "Effective moment of inertia per NSCP 2015",
         "category": "Analysis"},
    ]
    return render_template("index.html", calculators=calculators)


@main_bp.route("/calculator/<calc_id>")
def calculator(calc_id):
    """Render a specific calculator page."""
    return render_template(f"calculators/{calc_id}.html")
