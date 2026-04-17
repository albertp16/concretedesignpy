# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""Section, detailing, and analysis API routes."""

from flask import Blueprint, request, jsonify

from concretedesignpy.calculators.development_length import (
    bend_and_hook_deformed_bars,
    bend_and_hook_stirrups,
)
from concretedesignpy.calculators.moment_curvature import (
    moment_curvature_analysis,
    moment_curvature_advanced,
    moment_curvature_opensees,
    moment_rotation_from_mphi,
    plastic_hinge_length,
    asce41_backbone_beam,
)
from concretedesignpy.calculators.mander import confined_stress_strain
from concretedesignpy.calculators.rebar_layout import section_generate_rect
from concretedesignpy.calculators.alternative_inertia import (
    alternative_inertia_column_wall,
    alternative_inertia_beam_slab,
)
from concretedesignpy.calculators.diagrams import svg_hook_geometry

section_bp = Blueprint("section", __name__)


@section_bp.route("/development-length", methods=["POST"])
def development_length():
    """Hook geometry calculation."""
    data = request.get_json()
    try:
        bar_type = data.get("bar_type", "deformed")
        if bar_type == "deformed":
            result = bend_and_hook_deformed_bars(
                bar_size=float(data["bar_size"]),
                angle=int(data["angle"]),
            )
        else:
            result = bend_and_hook_stirrups(
                bar_size=float(data["bar_size"]),
                angle=int(data["angle"]),
            )
        # Add SVG diagram
        result["svg"] = svg_hook_geometry(
            result["bend_diameter"], result["lext"],
            result["bar_size"], result["angle"], result["ldh"],
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/moment-curvature", methods=["POST"])
def moment_curvature():
    """Moment-curvature analysis."""
    data = request.get_json()
    try:
        result = moment_curvature_analysis(
            b=float(data["b"]),
            h=float(data["h"]),
            d=float(data["d"]),
            fc=float(data["fc"]),
            fy=float(data["fy"]),
            as_tension=float(data["as_tension"]),
            es=float(data.get("es", 200000)),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


def _parse_section_vertices(data):
    """Parse an optional list of (x, y) vertices from the request payload."""
    raw = data.get("section_vertices")
    if not raw:
        return None
    verts = [(float(v[0]), float(v[1])) for v in raw]
    if len(verts) < 3:
        return None
    return verts


@section_bp.route("/moment-curvature-advanced", methods=["POST"])
def moment_curvature_adv():
    """Advanced moment-curvature with selectable concrete model."""
    data = request.get_json()
    try:
        b_val = float(data["b"])
        h_val = float(data["h"])
        fc_val = float(data["fc"])

        # Compression steel (optional)
        d_prime = float(data.get("d_prime", 0))
        as_compression = float(data.get("as_compression", 0))

        # Concrete model selection
        concrete_model = data.get("concrete_model", "hognestad")
        mander_params = None
        if concrete_model == "mander":
            mander_params = confined_stress_strain(
                fc=fc_val,
                fy_transverse=float(data.get("fy_transverse", 275)),
                b=b_val,
                h=h_val,
                cover=float(data.get("cover", 40)),
                db_main=float(data.get("db_tension", 25)),
                db_tie=float(data.get("db_tie", 10)),
                n_bars_x=int(data.get("n_bars_x", 4)),
                n_bars_y=int(data.get("n_bars_y", 4)),
                s_tie=float(data.get("s_tie", 100)),
            )

        section_vertices = _parse_section_vertices(data)

        result = moment_curvature_advanced(
            b=b_val,
            h=h_val,
            d=float(data["d"]),
            fc=fc_val,
            fy=float(data["fy"]),
            as_tension=float(data["as_tension"]),
            es=float(data.get("es", 200000)),
            axial_load=float(data.get("axial_load", 0)),
            d_prime=d_prime,
            as_compression=as_compression,
            concrete_model=concrete_model,
            mander_params=mander_params,
            section_vertices=section_vertices,
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/moment-curvature-opensees", methods=["POST"])
def moment_curvature_ops():
    """Moment-curvature via OpenSeesPy zero-length section element."""
    data = request.get_json()
    try:
        b_val = float(data["b"])
        h_val = float(data["h"])
        fc_val = float(data["fc"])

        concrete_model = data.get("concrete_model", "hognestad")
        mander_params = None
        if concrete_model == "mander":
            mander_params = confined_stress_strain(
                fc=fc_val,
                fy_transverse=float(data.get("fy_transverse", 275)),
                b=b_val, h=h_val,
                cover=float(data.get("cover", 40)),
                db_main=float(data.get("db_tension", 25)),
                db_tie=float(data.get("db_tie", 10)),
                n_bars_x=int(data.get("n_bars_x", 4)),
                n_bars_y=int(data.get("n_bars_y", 4)),
                s_tie=float(data.get("s_tie", 100)),
            )

        n_tension = int(data.get("n_tension", 4))
        n_comp = int(data.get("n_comp", 0))
        cover_val = float(data.get("cover", 40))
        d_prime = float(data.get("d_prime", 0))
        as_comp = float(data.get("as_compression", 0))

        result = moment_curvature_opensees(
            b=b_val, h=h_val,
            d=float(data["d"]),
            fc=fc_val,
            fy=float(data["fy"]),
            as_tension=float(data["as_tension"]),
            es=float(data.get("es", 200000)),
            axial_load=float(data.get("axial_load", 0)),
            n_increments=100,
            d_prime=d_prime,
            as_compression=as_comp,
            cover=cover_val,
            n_bars_tension=n_tension,
            n_bars_comp=n_comp,
            concrete_model=concrete_model,
            mander_params=mander_params,
        )
        return jsonify({"status": "success", "result": result})
    except ImportError as e:
        return jsonify({"status": "error", "message": str(e)}), 400
    except (KeyError, ValueError, TypeError, RuntimeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/moment-rotation", methods=["POST"])
def moment_rotation():
    """Moment-rotation analysis (M-theta from M-phi)."""
    data = request.get_json()
    try:
        b_val = float(data["b"])
        h_val = float(data["h"])
        fc_val = float(data["fc"])
        fy_val = float(data["fy"])
        d_prime = float(data.get("d_prime", 0))
        as_compression = float(data.get("as_compression", 0))

        # Concrete model
        concrete_model = data.get("concrete_model", "hognestad")
        mander_params = None
        if concrete_model == "mander":
            mander_params = confined_stress_strain(
                fc=fc_val,
                fy_transverse=float(data.get("fy_transverse", 275)),
                b=b_val, h=h_val,
                cover=float(data.get("cover", 40)),
                db_main=float(data.get("db_tension", 25)),
                db_tie=float(data.get("db_tie", 10)),
                n_bars_x=int(data.get("n_bars_x", 4)),
                n_bars_y=int(data.get("n_bars_y", 4)),
                s_tie=float(data.get("s_tie", 100)),
            )

        section_vertices = _parse_section_vertices(data)

        # Run M-phi first
        mphi_result = moment_curvature_advanced(
            b=b_val, h=h_val,
            d=float(data["d"]),
            fc=fc_val, fy=fy_val,
            as_tension=float(data["as_tension"]),
            es=float(data.get("es", 200000)),
            axial_load=float(data.get("axial_load", 0)),
            d_prime=d_prime,
            as_compression=as_compression,
            concrete_model=concrete_model,
            mander_params=mander_params,
            section_vertices=section_vertices,
        )

        # Plastic hinge length
        lp_method = data.get("lp_method", "manual")
        if lp_method == "manual":
            lp = float(data.get("lp", 300))
        else:
            lp = plastic_hinge_length(
                method=lp_method,
                h=h_val,
                fy=fy_val,
                db=float(data.get("db_tension", 25)),
                shear_span=float(data.get("shear_span", 0)),
            )

        # Convert M-phi to M-theta
        mtheta_result = moment_rotation_from_mphi(mphi_result, lp)

        # ASCE 41 / FEMA 356 backbone curve
        backbone_result = None
        backbone_code = data.get("backbone_code", "none")
        if backbone_code in ("asce41", "fema356"):
            duct = mtheta_result.get("ductility_rotation")
            if duct:
                backbone_result = asce41_backbone_beam(
                    my_knm=duct["moment_yield_knm"],
                    theta_y_rad=duct["theta_y_rad"],
                    fc=fc_val, fy=fy_val, b=b_val, d=float(data["d"]),
                    as_tension=float(data["as_tension"]),
                    as_compression=float(data.get("as_compression", 0)),
                    stirrup_spacing=float(data.get("stirrup_spacing", 150)),
                    db_stirrup=float(data.get("db_stirrup", 10)),
                    v_demand=float(data.get("v_demand", 0)),
                    code=backbone_code,
                )

        # Return all
        return jsonify({
            "status": "success",
            "result": {
                "mphi": mphi_result,
                "mtheta": mtheta_result,
                "backbone": backbone_result,
            },
        })
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/rebar-layout", methods=["POST"])
def rebar_layout():
    """Generate rebar layout coordinates."""
    data = request.get_json()
    try:
        result = section_generate_rect(
            width=float(data["width"]),
            height=float(data["height"]),
            cover=float(data["cover"]),
            ds=float(data["ds"]),
            db=float(data["db"]),
            nx=int(data["nx"]),
            ny=int(data["ny"]),
        )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400


@section_bp.route("/alternative-inertia", methods=["POST"])
def alternative_inertia():
    """Alternative moment of inertia calculation."""
    data = request.get_json()
    try:
        member_type = data.get("member_type", "column")
        if member_type in ("column", "wall"):
            result = alternative_inertia_column_wall(
                ig=float(data["ig"]),
                ast=float(data["ast"]),
                ag=float(data["ag"]),
                mu=float(data["mu"]),
                pu=float(data["pu"]),
                h=float(data["h"]),
                po=float(data["po"]),
            )
        else:
            result = alternative_inertia_beam_slab(
                ig=float(data["ig"]),
                rho=float(data["rho"]),
                b=float(data["b"]),
                d=float(data["d"]),
            )
        return jsonify({"status": "success", "result": result})
    except (KeyError, ValueError, TypeError) as e:
        return jsonify({"status": "error", "message": str(e)}), 400
