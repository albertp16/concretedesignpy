# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Beam calculation API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.beam_moment import calculate_beam_moment
from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength,
    compute_steel_shear_strength,
    compute_shear_spacing,
    shear_torsion_design,
)
from concretedesignpy.calculators.beam_torsion import torsion_design
from concretedesignpy.calculators.beam_deflection import deflection_computation
from concretedesignpy.calculators.diagrams import (
    svg_beam_cross_section,
    svg_shear_diagram,
    svg_torsion_section,
)

beam_bp = Blueprint("beam", __name__)


@beam_bp.route("/moment", methods=["POST"])
def beam_moment():
    """Calculate beam moment capacity."""
    data = request.get_json()
    try:
        result = calculate_beam_moment(
            rebar_list=data["rebar_list"],
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            b=float(data["b"]),
            h=float(data["h"]),
            es=float(data.get("es", 200000)),
        )
        # Generate cross-section SVG
        result["svg"] = svg_beam_cross_section(
            b=float(data["b"]),
            h=float(data["h"]),
            rebar_forces=result["rebar_forces"],
            c=result["neutral_axis"],
            a=result["a"],
        )
        # Round rebar forces for display
        for rf in result["rebar_forces"]:
            rf["strain"] = round(rf["strain"], 6)
            rf["stress"] = round(rf["stress"], 2)
            rf["force"] = round(rf["force"] / 1000, 2)  # Convert to kN
            rf["area"] = round(rf["area"], 2)
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/shear", methods=["POST"])
def beam_shear():
    """Calculate beam shear capacity or required spacing."""
    data = request.get_json()
    try:
        mode = data.get("mode", "capacity")
        phi = float(data.get("phi", 0.75))

        if mode == "capacity":
            vc_result = compute_concrete_shear_strength(
                fc=float(data["fc"]),
                b=float(data["b"]),
                d=float(data["d"]),
                lamda=float(data.get("lamda", 1.0)),
                vc_type=data.get("vc_type", "simple"),
            )
            vs_result = compute_steel_shear_strength(
                av=float(data["av"]),
                fyt=float(data["fyt"]),
                d=float(data["d"]),
                s=float(data["s"]),
            )
            vn = vc_result["vc"] + vs_result["vs"]
            vu = phi * vn
            result = {
                **vc_result,
                **vs_result,
                "vn": round(vn, 2),
                "vn_kn": round(vn / 1000, 2),
                "vu": round(vu, 2),
                "vu_kn": round(vu / 1000, 2),
            }
            result["svg"] = svg_shear_diagram(
                vc_result["vc_kn"], vs_result["vs_kn"],
                round(vu / 1000, 2), phi,
            )
        else:
            result = compute_shear_spacing(
                fc=float(data["fc"]),
                b=float(data["b"]),
                d=float(data["d"]),
                fyt=float(data["fyt"]),
                vu_required=float(data["vu_required"]),
                phi=float(data.get("phi", 0.75)),
                av=float(data["av"]),
                lamda=float(data.get("lamda", 1.0)),
                vc_type=data.get("vc_type", "simple"),
            )

        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/shear-torsion", methods=["POST"])
def beam_shear_torsion():
    """Combined shear and torsion design per ACI 318M-14."""
    data = request.get_json()
    try:
        result = shear_torsion_design(
            fc=float(data["fc"]),
            fyv=float(data["fyv"]),
            fy=float(data["fy"]),
            phi=float(data.get("phi", 0.75)),
            bw=float(data["bw"]),
            h=float(data["h"]),
            cc=float(data["cc"]),
            c=float(data["c"]),
            d=float(data["d"]),
            vu=float(data["vu"]),
            tu=float(data.get("tu", 0)),
            nu=float(data.get("nu", 0)),
            s_chosen=float(data.get("s_chosen", 150)),
            n_legs=int(data.get("n_legs", 4)),
            db_stirrup=float(data.get("db_stirrup", 10)),
            db_long=float(data.get("db_long", 12)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/torsion", methods=["POST"])
def beam_torsion():
    """Perform torsion design checks."""
    data = request.get_json()
    try:
        result = torsion_design(
            width=float(data["width"]),
            height=float(data["height"]),
            cover=float(data["cover"]),
            db=float(data["db"]),
            tf=float(data.get("tf", 0)),
            beff=float(data.get("beff", data["width"])),
            phi_torsion=float(data.get("phi_torsion", 0.75)),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            tu=float(data["tu"]),
            vc=float(data["vc"]),
            ds=float(data["ds"]),
            smax_shear=float(data["smax_shear"]),
            s_actual=float(data["s_actual"]),
            av=float(data["av"]),
            s=float(data["s"]),
        )
        g = result["geometry"]
        result["svg"] = svg_torsion_section(
            float(data["width"]), float(data["height"]),
            g["x1"], g["y1"], g["aoh"], g["ph"],
            float(data["cover"]), float(data["db"]),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@beam_bp.route("/deflection", methods=["POST"])
def beam_deflection():
    """Calculate beam deflection."""
    data = request.get_json()
    try:
        result = deflection_computation(
            b=float(data["b"]),
            h=float(data["h"]),
            d=float(data["d"]),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            clearspan=float(data["clearspan"]),
            as_tension=float(data["as_tension"]),
            n_bars_comp=int(data.get("n_bars_comp", 0)),
            db_comp=float(data.get("db_comp", 0)),
            es=float(data.get("es", 200000)),
            ma=float(data["ma"]) if data.get("ma") else None,
            point_load=float(data.get("point_load", 0)),
            uniform_load=float(data.get("uniform_load", 0)),
            beam_type=data.get("beam_type", "simply_supported"),
            member_type=data.get("member_type", "floor"),
            sustained_duration=data.get("sustained_duration", "5_years_or_more"),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
