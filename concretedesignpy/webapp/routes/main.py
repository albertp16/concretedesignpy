# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Main page routes."""

from flask import Blueprint, abort, render_template

from concretedesignpy.webapp import nav

main_bp = Blueprint("main", __name__)


@main_bp.route("/")
def index():
    """Landing page with calculator directory."""
    return render_template(
        "index.html",
        groups=nav.calculators_by_category(),
        calculator_count=len(nav.CALCULATORS),
    )


@main_bp.route("/calculator/<calc_id>")
def calculator(calc_id):
    """Render a specific calculator page."""
    calc = nav.get_calculator(calc_id)
    if calc is None:
        abort(404)
    return render_template(f"calculators/{calc_id}.html", calc=calc)
