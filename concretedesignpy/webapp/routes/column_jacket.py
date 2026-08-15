# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""RC column concrete jacketing API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.column_jacket_design import column_jacket_design

column_jacket_bp = Blueprint("column_jacket", __name__)


@column_jacket_bp.route("/design", methods=["POST"])
def design():
    """Full jacketing design report for an existing rectangular column.

    Unlike the other routes, this one does NOT coerce each field with
    float()/int() here. ``column_jacket_design`` carries its own validation
    layer -- bounds, the bar-cage lever-arm check, the jacket-thickness check
    and the model-scope refusals -- and that layer is the safety boundary. Its
    ValueError message is the exact refusal and is passed through verbatim, so
    the caller is told which field is out of range and why, rather than
    "invalid input".
    """
    data = request.get_json()
    try:
        result = column_jacket_design(
            existing=data["existing"],
            jacket=data["jacket"],
            demand=data["demand"],
            construction=data.get("construction"),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
