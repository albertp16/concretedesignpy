# MIT License
# Copyright (c) Albert Pamonag Engineering Consultancy

"""
SVG Diagram Generators
========================

Generates inline SVG diagrams for structural engineering visualizations.
Used by the web app to render beam cross-sections, hook geometry,
interaction diagrams, and torsion/shear diagrams.
"""

import math


def svg_hook_geometry(bend_diameter, lext, bar_size, angle, ldh=None):
    """
    Generate SVG for standard hook geometry.

    Returns
    -------
    str : SVG markup string
    """
    if ldh is None:
        ldh = lext + bend_diameter

    r = bend_diameter / 2.0
    w = 420
    h = 340
    pad = 40

    # Scale to fit
    max_dim = max(ldh + bar_size, lext + r + bar_size)
    scale = min((w - 2 * pad) / (ldh + bar_size + 40),
                (h - 2 * pad) / (lext + r + bar_size + 40))
    scale = min(scale, 2.0)

    ox = pad + 20
    oy = pad + 20

    svg = f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" style="max-width:420px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Horizontal bar (ldh)
    hx1 = ox
    hx2 = ox + ldh * scale
    hy = oy + bar_size * scale / 2

    svg += f'<rect x="{hx1}" y="{oy}" width="{ldh * scale}" height="{bar_size * scale}" '
    svg += 'fill="#d4d4d8" stroke="#333" stroke-width="1.5"/>\n'

    if angle == 90:
        # Arc at corner
        arc_cx = hx2 - r * scale
        arc_cy = hy + r * scale
        arc_r = r * scale

        # Vertical bar going down
        vx = hx2
        vy = hy
        svg += f'<rect x="{vx}" y="{vy}" width="{bar_size * scale}" height="{lext * scale}" '
        svg += 'fill="#d4d4d8" stroke="#333" stroke-width="1.5"/>\n'

        # Dashed bend circle
        svg += f'<circle cx="{arc_cx}" cy="{arc_cy}" r="{arc_r}" '
        svg += 'fill="none" stroke="#3b82f6" stroke-width="1" stroke-dasharray="4,3"/>\n'

        # Dimension: ldh
        dim_y = oy - 15
        svg += _svg_dim_line(hx1, dim_y, hx2, dim_y, f"ldh = {ldh:.0f} mm")
        # Dimension: lext
        dim_x = vx + bar_size * scale + 15
        svg += _svg_dim_line_v(dim_x, vy, dim_x, vy + lext * scale, f"lext = {lext:.0f} mm")
        # Dimension: bend
        svg += f'<text x="{arc_cx}" y="{arc_cy + arc_r + 14}" text-anchor="middle" '
        svg += f'font-size="11" fill="#3b82f6">D = {bend_diameter:.0f} mm</text>\n'

    elif angle == 180:
        # U-bend: arc goes back
        arc_cx = hx2
        arc_cy = hy
        arc_r = r * scale

        # Return bar going left
        ret_len = lext * scale
        svg += f'<rect x="{hx2 - ret_len}" y="{hy + 2 * arc_r - bar_size * scale}" '
        svg += f'width="{ret_len}" height="{bar_size * scale}" '
        svg += 'fill="#d4d4d8" stroke="#333" stroke-width="1.5"/>\n'

        # Semicircle
        svg += f'<path d="M {hx2} {hy} A {arc_r} {arc_r} 0 0 1 {hx2} {hy + 2 * arc_r}" '
        svg += 'fill="none" stroke="#3b82f6" stroke-width="1.5" stroke-dasharray="4,3"/>\n'

        # Dimensions
        dim_y = oy - 15
        svg += _svg_dim_line(hx1, dim_y, hx2, dim_y, f"ldh = {ldh:.0f} mm")
        dim_x = hx2 + arc_r + 15
        svg += _svg_dim_line_v(dim_x, hy, dim_x, hy + 2 * arc_r,
                               f"D = {bend_diameter:.0f} mm")
        # lext
        lext_x1 = hx2 - ret_len
        lext_y = hy + 2 * arc_r + bar_size * scale + 15
        svg += _svg_dim_line(lext_x1, lext_y, hx2, lext_y, f"lext = {lext:.0f} mm")

    else:
        # 135 degree - simplified as 90 with angled tail
        arc_cx = hx2 - r * scale
        arc_cy = hy + r * scale
        arc_r = r * scale
        svg += f'<circle cx="{arc_cx}" cy="{arc_cy}" r="{arc_r}" '
        svg += 'fill="none" stroke="#3b82f6" stroke-width="1" stroke-dasharray="4,3"/>\n'
        # Angled bar
        end_x = hx2 + lext * scale * 0.707
        end_y = hy + lext * scale * 0.707
        svg += f'<rect x="{hx2}" y="{hy}" width="{bar_size * scale}" height="{lext * scale}" '
        svg += f'fill="#d4d4d8" stroke="#333" stroke-width="1.5" '
        svg += f'transform="rotate(45 {hx2} {hy})"/>\n'
        dim_y = oy - 15
        svg += _svg_dim_line(hx1, dim_y, hx2, dim_y, f"ldh = {ldh:.0f} mm")

    # Title
    svg += f'<text x="{w / 2}" y="{h - 8}" text-anchor="middle" '
    svg += f'font-size="12" fill="#666">{angle}&deg; Hook &mdash; &empty;{bar_size:.0f} mm</text>\n'

    svg += '</svg>'
    return svg


def svg_beam_cross_section(b, h, rebar_forces, c, a, scale_factor=None):
    """
    Generate SVG of beam cross-section with neutral axis and rebar.

    Parameters
    ----------
    b, h : float
        Beam width and height (mm).
    rebar_forces : list of dict
        Each has 'd', 'area', 'force' keys.
    c : float
        Neutral axis depth (mm).
    a : float
        Whitney block depth (mm).

    Returns
    -------
    str : SVG markup
    """
    w = 320
    ht = 420
    pad = 50

    if scale_factor is None:
        scale_factor = min((w - 2 * pad) / b, (ht - 2 * pad) / h)

    bw = b * scale_factor
    bh = h * scale_factor
    ox = (w - bw) / 2
    oy = pad

    svg = f'<svg viewBox="0 0 {w} {ht}" xmlns="http://www.w3.org/2000/svg" style="max-width:320px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Concrete section
    svg += f'<rect x="{ox}" y="{oy}" width="{bw}" height="{bh}" '
    svg += 'fill="#e5e7eb" stroke="#333" stroke-width="2"/>\n'

    # Whitney stress block
    ah = a * scale_factor
    svg += f'<rect x="{ox}" y="{oy}" width="{bw}" height="{ah}" '
    svg += 'fill="rgba(59,130,246,0.2)" stroke="#3b82f6" stroke-width="1" stroke-dasharray="4,2"/>\n'
    svg += f'<text x="{ox + bw + 5}" y="{oy + ah / 2 + 4}" font-size="10" fill="#3b82f6">a={a:.1f}</text>\n'

    # Neutral axis
    cy = oy + c * scale_factor
    svg += f'<line x1="{ox - 10}" y1="{cy}" x2="{ox + bw + 10}" y2="{cy}" '
    svg += 'stroke="#dc2626" stroke-width="1.5" stroke-dasharray="6,3"/>\n'
    svg += f'<text x="{ox + bw + 5}" y="{cy - 4}" font-size="10" fill="#dc2626">c={c:.1f}</text>\n'

    # Rebar
    for rb in rebar_forces:
        ry = oy + rb["d"] * scale_factor
        # Distribute bars across width
        n_bars = max(1, round(rb["area"] / (math.pi * 64)))  # approx count assuming ~16mm
        bar_r = max(4, min(8, 6))
        if n_bars == 1:
            positions = [ox + bw / 2]
        else:
            margin = 25
            spacing = (bw - 2 * margin) / (n_bars - 1) if n_bars > 1 else 0
            positions = [ox + margin + i * spacing for i in range(n_bars)]

        color = "#16a34a" if rb.get("force", 0) > 0 else "#dc2626"
        for px in positions:
            svg += f'<circle cx="{px}" cy="{ry}" r="{bar_r}" fill="{color}" stroke="#333" stroke-width="1"/>\n'

        svg += f'<text x="{ox - 5}" y="{ry + 4}" text-anchor="end" font-size="9" fill="#666">d={rb["d"]:.0f}</text>\n'

    # Dimension labels
    # Width
    dy = oy + bh + 20
    svg += _svg_dim_line(ox, dy, ox + bw, dy, f"b = {b:.0f} mm")
    # Height
    dx = ox - 20
    svg += _svg_dim_line_v(dx, oy, dx, oy + bh, f"h = {h:.0f} mm")

    svg += '</svg>'
    return svg


def svg_interaction_diagram(points, title="P-M Interaction Diagram", demand=None):
    """
    Generate SVG of P-M interaction diagram.

    Parameters
    ----------
    points : list of dict
        Each has 'mu' and 'pu' keys (factored values).
    title : str
    demand : tuple or None
        (Pu_kN, Mu_kNm) demand point to plot.

    Returns
    -------
    str : SVG markup
    """
    w = 560
    h = 440
    pad = 60

    if not points:
        return '<svg viewBox="0 0 560 440"><text x="280" y="220" text-anchor="middle">No data</text></svg>'

    mu_vals = [abs(p["mu"]) for p in points]
    pu_vals = [p["pu"] for p in points]

    max_m = max(mu_vals) * 1.15 if max(mu_vals) > 0 else 1
    max_p = max(pu_vals) * 1.1 if max(pu_vals) > 0 else 1
    min_p = min(min(pu_vals) * 1.1, 0)

    # Include demand point in range
    if demand is not None:
        pu_d, mu_d = demand
        max_m = max(max_m, abs(mu_d) * 1.15)
        max_p = max(max_p, pu_d * 1.1)
        min_p = min(min_p, pu_d * 1.1)

    range_p = max_p - min_p if max_p != min_p else 1

    pw = w - 2 * pad
    ph = h - 2 * pad

    def tx(m):
        return pad + (abs(m) / max_m) * pw

    def ty(p):
        return pad + ph - ((p - min_p) / range_p) * ph

    svg = f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" style="max-width:560px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Grid
    for i in range(6):
        gy = pad + i * ph / 5
        svg += f'<line x1="{pad}" y1="{gy}" x2="{pad + pw}" y2="{gy}" stroke="#e5e7eb" stroke-width="0.5"/>\n'
        pval = max_p - i * range_p / 5 + min_p
        svg += f'<text x="{pad - 5}" y="{gy + 4}" text-anchor="end" font-size="9" fill="#999">{pval:.0f}</text>\n'
    for i in range(6):
        gx = pad + i * pw / 5
        svg += f'<line x1="{gx}" y1="{pad}" x2="{gx}" y2="{pad + ph}" stroke="#e5e7eb" stroke-width="0.5"/>\n'
        mval = i * max_m / 5
        svg += f'<text x="{gx}" y="{pad + ph + 15}" text-anchor="middle" font-size="9" fill="#999">{mval:.0f}</text>\n'

    # Axes
    svg += f'<line x1="{pad}" y1="{pad}" x2="{pad}" y2="{pad + ph}" stroke="#333" stroke-width="1.5"/>\n'
    svg += f'<line x1="{pad}" y1="{pad + ph}" x2="{pad + pw}" y2="{pad + ph}" stroke="#333" stroke-width="1.5"/>\n'

    # Zero axis if visible
    if min_p < 0:
        zy = ty(0)
        svg += f'<line x1="{pad}" y1="{zy}" x2="{pad + pw}" y2="{zy}" stroke="#999" stroke-width="0.5" stroke-dasharray="4,2"/>\n'

    # Nominal curve (lighter)
    path_n = "M"
    for pt in points:
        x = tx(pt.get("mn", pt.get("mu", 0)))
        y = ty(pt.get("pn", pt.get("pu", 0)))
        path_n += f" {x:.1f},{y:.1f}"
    svg += f'<path d="{path_n}" fill="none" stroke="#93c5fd" stroke-width="1.5" stroke-dasharray="4,2"/>\n'

    # Factored curve (φPn, φMn)
    path_f = "M"
    for pt in points:
        x = tx(pt["mu"])
        y = ty(pt["pu"])
        path_f += f" {x:.1f},{y:.1f}"
    svg += f'<path d="{path_f}" fill="none" stroke="#2563eb" stroke-width="2.5"/>\n'

    # Points on factored curve
    for pt in points:
        x = tx(pt["mu"])
        y = ty(pt["pu"])
        svg += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="2.5" fill="#2563eb"/>\n'

    # Demand point
    if demand is not None:
        pu_d, mu_d = demand
        dx = tx(mu_d)
        dy = ty(pu_d)
        svg += f'<circle cx="{dx:.1f}" cy="{dy:.1f}" r="6" fill="none" stroke="#dc2626" stroke-width="2"/>\n'
        svg += f'<circle cx="{dx:.1f}" cy="{dy:.1f}" r="2" fill="#dc2626"/>\n'
        svg += f'<text x="{dx + 10}" y="{dy - 4}" font-size="10" fill="#dc2626" font-weight="bold">'
        svg += f'Pu={pu_d:.0f}, Mu={mu_d:.0f}</text>\n'

    # Legend
    ly = pad + 10
    svg += f'<line x1="{pad + pw - 120}" y1="{ly}" x2="{pad + pw - 100}" y2="{ly}" stroke="#2563eb" stroke-width="2.5"/>\n'
    svg += f'<text x="{pad + pw - 95}" y="{ly + 4}" font-size="9" fill="#333">\u03c6Pn, \u03c6Mn</text>\n'
    svg += f'<line x1="{pad + pw - 120}" y1="{ly + 15}" x2="{pad + pw - 100}" y2="{ly + 15}" stroke="#93c5fd" stroke-width="1.5" stroke-dasharray="4,2"/>\n'
    svg += f'<text x="{pad + pw - 95}" y="{ly + 19}" font-size="9" fill="#333">Pn, Mn</text>\n'

    # Labels
    svg += f'<text x="{pad + pw / 2}" y="{pad + ph + 35}" text-anchor="middle" font-size="12" fill="#333">'
    svg += f'Moment, \u03c6Mn (kN-m)</text>\n'
    svg += f'<text x="14" y="{pad + ph / 2}" text-anchor="middle" font-size="12" fill="#333" '
    svg += f'transform="rotate(-90 14 {pad + ph / 2})">Axial, \u03c6Pn (kN)</text>\n'
    svg += f'<text x="{w / 2}" y="20" text-anchor="middle" font-size="13" font-weight="bold" fill="#1e3a8a">{title}</text>\n'

    svg += '</svg>'
    return svg


def matplotlib_interaction_diagram(points, title="P-M Interaction Diagram", demand=None):
    """
    Generate P-M interaction diagram using matplotlib, returned as base64 PNG.

    Parameters
    ----------
    points : list of dict
        Each has 'mu', 'pu', 'mn', 'pn' keys (factored and nominal values).
    title : str
    demand : tuple or None
        (Pu_kN, Mu_kNm) demand point to plot.

    Returns
    -------
    str : base64-encoded PNG image data URI
    """
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt
    import io
    import base64

    fig, ax = plt.subplots(figsize=(7, 5.5))

    # Extract data
    mu_vals = [abs(p["mu"]) for p in points]
    pu_vals = [p["pu"] for p in points]
    mn_vals = [abs(p.get("mn", p["mu"])) for p in points]
    pn_vals = [p.get("pn", p["pu"]) for p in points]

    # Plot nominal curve (dashed, lighter)
    ax.plot(mn_vals, pn_vals, color="#93c5fd", linewidth=1.5, linestyle="--",
            label="Pn, Mn (nominal)", zorder=2)

    # Plot factored curve (solid, bold)
    ax.plot(mu_vals, pu_vals, color="#1e40af", linewidth=2.5,
            label="\u03c6Pn, \u03c6Mn (factored)", zorder=3)
    ax.scatter(mu_vals, pu_vals, color="#1e40af", s=12, zorder=4)

    # Demand point
    if demand is not None:
        pu_d, mu_d = demand
        ax.plot(abs(mu_d), pu_d, "o", color="#dc2626", markersize=10,
                markerfacecolor="none", markeredgewidth=2.5, zorder=5)
        ax.plot(abs(mu_d), pu_d, "o", color="#dc2626", markersize=3, zorder=6)
        ax.annotate(f"Pu={pu_d:.0f}, Mu={abs(mu_d):.0f}",
                    xy=(abs(mu_d), pu_d), xytext=(12, 8),
                    textcoords="offset points", fontsize=9, color="#dc2626",
                    fontweight="bold")

    # Zero line
    ax.axhline(y=0, color="#999", linewidth=0.5, linestyle="--")

    # Labels and title
    ax.set_xlabel("Moment, \u03c6Mn (kN-m)", fontsize=11)
    ax.set_ylabel("Axial, \u03c6Pn (kN)", fontsize=11)
    ax.set_title(title, fontsize=13, fontweight="bold", color="#1e3a8a")
    ax.legend(loc="upper right", fontsize=9)
    ax.grid(True, alpha=0.3)
    ax.set_xlim(left=0)

    fig.tight_layout()

    # Convert to base64
    buf = io.BytesIO()
    fig.savefig(buf, format="png", dpi=150, bbox_inches="tight",
                facecolor="white", edgecolor="none")
    plt.close(fig)
    buf.seek(0)
    b64 = base64.b64encode(buf.read()).decode("ascii")
    return f"data:image/png;base64,{b64}"


def svg_rebar_section(section_data, width, height):
    """
    Generate SVG of rebar layout in rectangular section.

    Parameters
    ----------
    section_data : dict
        Output from section_generate_rect.
    width, height : float
        Section dimensions (mm).

    Returns
    -------
    str : SVG markup
    """
    w = 300
    ht = 400
    pad = 40

    scale = min((w - 2 * pad) / width, (ht - 2 * pad) / height)
    bw = width * scale
    bh = height * scale
    ox = (w - bw) / 2
    oy = (ht - bh) / 2

    def tx(x):
        return ox + x * scale

    def ty(y):
        return oy + (height - y) * scale

    svg = f'<svg viewBox="0 0 {w} {ht}" xmlns="http://www.w3.org/2000/svg" style="max-width:300px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Outer section
    svg += f'<rect x="{ox}" y="{oy}" width="{bw}" height="{bh}" fill="#f3f4f6" stroke="#333" stroke-width="2"/>\n'

    # Stirrup
    sc = section_data["stirrup_corners"]
    path = "M"
    for i, pt in enumerate(sc):
        x, y = tx(pt[0]), ty(pt[1])
        path += f" {x:.1f},{y:.1f}"
    path += " Z"
    svg += f'<path d="{path}" fill="none" stroke="#3b82f6" stroke-width="1.5" stroke-dasharray="5,3"/>\n'

    # Rebar
    for i, (rx, ry) in enumerate(section_data["rebar_coords"]):
        x, y = tx(rx), ty(ry)
        svg += f'<circle cx="{x:.1f}" cy="{y:.1f}" r="6" fill="#dc2626" stroke="#333" stroke-width="1"/>\n'
        svg += f'<text x="{x + 8}" y="{y - 5}" font-size="8" fill="#666">R{i + 1}</text>\n'

    # Dimensions
    svg += _svg_dim_line(ox, oy + bh + 15, ox + bw, oy + bh + 15, f"{width:.0f} mm")
    svg += _svg_dim_line_v(ox - 15, oy, ox - 15, oy + bh, f"{height:.0f} mm")

    # Count
    svg += f'<text x="{w / 2}" y="{ht - 5}" text-anchor="middle" font-size="11" fill="#333">'
    svg += f'Total: {section_data["total_bars"]} bars</text>\n'

    svg += '</svg>'
    return svg


def svg_torsion_section(width, height, x1, y1, aoh, ph, cover, db):
    """Generate SVG for torsion core geometry."""
    w = 300
    ht = 360
    pad = 40

    scale = min((w - 2 * pad) / width, (ht - 2 * pad) / height)
    bw = width * scale
    bh = height * scale
    ox = (w - bw) / 2
    oy = (ht - bh) / 2

    svg = f'<svg viewBox="0 0 {w} {ht}" xmlns="http://www.w3.org/2000/svg" style="max-width:300px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Outer
    svg += f'<rect x="{ox}" y="{oy}" width="{bw}" height="{bh}" fill="#f3f4f6" stroke="#333" stroke-width="2"/>\n'

    # Torsional core (x1, y1)
    core_offset = (cover + db / 2.0) * scale
    svg += f'<rect x="{ox + core_offset}" y="{oy + core_offset}" '
    svg += f'width="{x1 * scale}" height="{y1 * scale}" '
    svg += 'fill="rgba(59,130,246,0.1)" stroke="#3b82f6" stroke-width="1.5" stroke-dasharray="5,3"/>\n'

    # Labels
    mid_x = ox + bw / 2
    mid_y = oy + bh / 2
    svg += f'<text x="{mid_x}" y="{mid_y - 8}" text-anchor="middle" font-size="10" fill="#3b82f6">'
    svg += f'x1={x1:.0f}, y1={y1:.0f}</text>\n'
    svg += f'<text x="{mid_x}" y="{mid_y + 8}" text-anchor="middle" font-size="10" fill="#3b82f6">'
    svg += f'Aoh={aoh:.0f} mm&sup2;</text>\n'
    svg += f'<text x="{mid_x}" y="{mid_y + 22}" text-anchor="middle" font-size="10" fill="#3b82f6">'
    svg += f'ph={ph:.0f} mm</text>\n'

    svg += _svg_dim_line(ox, oy + bh + 15, ox + bw, oy + bh + 15, f"{width:.0f} mm")
    svg += _svg_dim_line_v(ox - 15, oy, ox - 15, oy + bh, f"{height:.0f} mm")

    svg += '</svg>'
    return svg


def svg_shear_diagram(vc_kn, vs_kn, vu_kn, phi):
    """Generate SVG bar chart for shear components."""
    w = 400
    h = 220
    pad = 50

    bars_data = [
        ("Vc", vc_kn, "#60a5fa"),
        ("Vs", vs_kn, "#34d399"),
        ("Vn", vc_kn + vs_kn, "#a78bfa"),
        (f"\u03c6Vn", phi * (vc_kn + vs_kn), "#f87171"),
    ]
    max_v = max(d[1] for d in bars_data) * 1.2 if any(d[1] for d in bars_data) else 1

    bw = 50
    gap = 20
    total_w = len(bars_data) * (bw + gap) - gap
    start_x = (w - total_w) / 2
    chart_h = h - 2 * pad

    svg = f'<svg viewBox="0 0 {w} {h}" xmlns="http://www.w3.org/2000/svg" style="max-width:400px;width:100%">\n'
    svg += '<rect width="100%" height="100%" fill="#fff"/>\n'

    # Axis
    svg += f'<line x1="{pad}" y1="{h - pad}" x2="{w - pad}" y2="{h - pad}" stroke="#333" stroke-width="1"/>\n'

    for i, (label, val, color) in enumerate(bars_data):
        x = start_x + i * (bw + gap)
        bar_h = (val / max_v) * chart_h if max_v > 0 else 0
        y = h - pad - bar_h

        svg += f'<rect x="{x}" y="{y}" width="{bw}" height="{bar_h}" fill="{color}" rx="3"/>\n'
        svg += f'<text x="{x + bw / 2}" y="{y - 5}" text-anchor="middle" font-size="10" fill="#333">{val:.1f}</text>\n'
        svg += f'<text x="{x + bw / 2}" y="{h - pad + 14}" text-anchor="middle" font-size="10" fill="#666">{label}</text>\n'

    svg += f'<text x="{w / 2}" y="16" text-anchor="middle" font-size="11" fill="#333">Shear Components (kN)</text>\n'
    svg += '</svg>'
    return svg


# --- Helpers ---

def _svg_dim_line(x1, y, x2, y2, label):
    """Horizontal dimension line with label."""
    mid = (x1 + x2) / 2
    s = f'<line x1="{x1}" y1="{y}" x2="{x2}" y2="{y2}" stroke="#333" stroke-width="0.8" marker-start="url(#arr)" marker-end="url(#arr)"/>\n'
    s += f'<defs><marker id="arr" viewBox="0 0 6 6" refX="3" refY="3" markerWidth="6" markerHeight="6" orient="auto"><path d="M0,0 L6,3 L0,6 Z" fill="#333"/></marker></defs>\n'
    s += f'<text x="{mid}" y="{y - 4}" text-anchor="middle" font-size="10" fill="#333">{label}</text>\n'
    return s


def _svg_dim_line_v(x, y1, x2, y2, label):
    """Vertical dimension line with label."""
    mid = (y1 + y2) / 2
    s = f'<line x1="{x}" y1="{y1}" x2="{x2}" y2="{y2}" stroke="#333" stroke-width="0.8"/>\n'
    s += f'<text x="{x - 4}" y="{mid}" text-anchor="end" font-size="10" fill="#333" '
    s += f'transform="rotate(-90 {x - 4} {mid})">{label}</text>\n'
    return s
