# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""FRP strengthening API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.frp_beam import frp_flexural_strengthening

frp_bp = Blueprint("frp", __name__)


@frp_bp.route("/beam", methods=["POST"])
def frp_beam():
    """FRP flexural strengthening of RC beam."""
    data = request.get_json()
    try:
        result = frp_flexural_strengthening(
            b=float(data["b"]),
            d=float(data["d"]),
            h=float(data["h"]),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            as_steel=float(data["as_steel"]),
            mu_required=float(data["mu_required"]),
            tf=float(data["tf"]),
            ffu_star=float(data["ffu_star"]),
            efu_star=float(data["efu_star"]),
            ef_frp=float(data["ef_frp"]),
            n_plies=int(data["n_plies"]),
            exposure=data.get("exposure", "interior"),
            fiber=data.get("fiber", "carbon"),
            mdl=float(data.get("mdl", 0)),
            k_initial=float(data.get("k_initial", 0.334)),
            psi_f=float(data.get("psi_f", 0.85)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
