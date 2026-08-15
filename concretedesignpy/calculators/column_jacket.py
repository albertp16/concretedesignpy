"""
column_jacket.py
================
Strength assessment of an existing rectangular RC column strengthened by a
cast-in-place reinforced concrete jacket.

Albert Pamonag Engineering Consultancy (APEC)
Technical note: TN-RET-001, RC Column Concrete Jacketing
Created 2026-08-07.  Rev00.

UNITS -- strictly SI, internal: mm, MPa, N, N*mm.
        Reporting helpers convert to kN and kN*m.

SCOPE
    Rectangular existing column, four-sided jacket of uniform thickness.
    Purpose of the jacket is to INCREASE STRENGTH (axial, flexural, shear)
    and confinement of the existing column.

CODE BASIS
    ACI 318-19 sections 10.6, 16.4, 17.6, 21.2, 22.2, 22.4, 22.5, 22.9,
    25.4, 25.7.  NSCP 2015 Chapter 4 adopts ACI 318-14; clause numbers
    differ, values used here do not.  SI forms taken from ACI 318-19
    Appendix C and cross-checked against ACI 318-25M.
    ACI 562-16 sections 5.4, 7.3, 7.4 for the existing-structure and
    interface-bond framework.
    Confined concrete after Mander et al. (1988) as printed in
    fib Bulletin 14 (2001) Eq. (6-5), (6-34), (6-35), (6-36).

IMPORTANT
    This module implements code equations.  It does not exercise engineering
    judgement, does not know your structure, and does not replace the
    Engineer of Record.  Read TN-RET-001 sections "Limitations" and
    "Common Errors" before using any number it produces.
"""

from __future__ import annotations

import math
from dataclasses import dataclass, field
from typing import List, Literal, Optional

import numpy as np

__all__ = [
    "Concrete", "Steel", "Bar", "TieSet", "JacketedColumn",
    "rect_perimeter_bars", "ConfinedConcrete", "capacity_at_P", "utilisation",
]

# --------------------------------------------------------------------------
# 1.  MATERIALS
# --------------------------------------------------------------------------


@dataclass
class Concrete:
    """Normal-weight concrete.

    fc     : specified / lower-bound / expected cylinder strength, MPa.
             Which of the three you pass is YOUR decision -- see TN-RET-001
             section "Material properties".  For an EXISTING column,
             ASCE 41-17 Table 10-1 multiplies lower-bound f'c by 1.50 to get
             the expected value; use lower-bound for force-controlled actions.
    Ec     : modulus of elasticity, MPa.  If None, ACI 318-19 Eq. (19.2.2.1.b)
             Ec = 4700*sqrt(f'c) is used (normal-weight).
    lam    : lightweight modification factor lambda, ACI 318-19 19.2.4.
             1.0 for normal weight.
    """
    fc: float
    Ec: Optional[float] = None
    lam: float = 1.0
    label: str = ""

    def __post_init__(self):
        if self.Ec is None:
            self.Ec = 4700.0 * math.sqrt(self.fc)

    @property
    def beta1(self) -> float:
        """beta1, ACI 318-19 Table 22.2.2.4.3, SI form (318-19 App. C /
        318-25M):  0.85 for f'c <= 28 MPa; 0.85 - 0.05*(f'c-28)/7;
        not less than 0.65."""
        if self.fc <= 28.0:
            return 0.85
        return max(0.65, 0.85 - 0.05 * (self.fc - 28.0) / 7.0)

    @property
    def eps_c0(self) -> float:
        """Strain at peak unconfined stress. Popovics/Mander form used by
        fib B14 Eq. (6-3); 0.002 is the conventional value for normal
        strength concrete and is used here."""
        return 0.002

    def sigma_unconfined(self, eps: np.ndarray) -> np.ndarray:
        """Unconfined compressive stress-strain, Popovics curve
        (fib Bulletin 14 Eq. (6-3)/(6-4) with f_cc = f_co, eps_cc = eps_co).
        Compression positive.  Zero in tension (tensile strength neglected,
        ACI 318-19 22.2.2.2)."""
        eps = np.asarray(eps, dtype=float)
        out = np.zeros_like(eps)
        m = eps > 0.0
        if not np.any(m):
            return out
        x = eps[m] / self.eps_c0
        Esec = self.fc / self.eps_c0
        r = self.Ec / max(self.Ec - Esec, 1e-6)
        out[m] = self.fc * x * r / (r - 1.0 + x ** r)
        # crushed beyond eps_cu_unconfined -> no stress
        out[m] = np.where(eps[m] > 0.0038, 0.0, out[m])
        return out


@dataclass
class ConfinedConcrete:
    """Mander et al. (1988) confined concrete, as printed in
    fib Bulletin 14 (2001) Eq. (6-5), (6-34), (6-35), (6-36).

        f_l    = 0.5 * k_e * rho_st * f_yh                      Eq. (6-35)
        f_cc   = f_co * [2.254*sqrt(1 + 7.94*f_l/f_co)
                         - 2*f_l/f_co - 1.254]                  Eq. (6-5)
        eps_cc = eps_co * [1 + 5*(f_cc/f_co - 1)]               Eq. (6-3 note)
        eps_cu = 0.004 + 1.4*rho_st*f_yh*eps_su/f_cc            Eq. (6-36)

    k_e is the confinement effectiveness (arching) coefficient.  fib B14
    states "usually 0.8" for steel hoops; Priestley/Calvi/Kowalsky give
    C_e = 0.75-0.85 for rectangular sections.  Default 0.80.
    """
    base: Concrete
    rho_st: float           # volumetric ratio of confining (jacket) ties
    fyh: float              # tie yield strength, MPa
    ke: float = 0.80
    eps_su: float = 0.09    # tie steel strain at peak stress

    def __post_init__(self):
        fco = self.base.fc
        self.fl = 0.5 * self.ke * self.rho_st * self.fyh
        r = self.fl / fco
        self.fcc = fco * (2.254 * math.sqrt(1.0 + 7.94 * r) - 2.0 * r - 1.254)
        self.eps_cc = self.base.eps_c0 * (1.0 + 5.0 * (self.fcc / fco - 1.0))
        self.eps_cu = 0.004 + 1.4 * self.rho_st * self.fyh * self.eps_su / self.fcc
        # a confined section can never be worse than unconfined
        self.fcc = max(self.fcc, fco)
        self.eps_cu = max(self.eps_cu, 0.0038)

    @property
    def strength_gain(self) -> float:
        return self.fcc / self.base.fc

    def sigma(self, eps: np.ndarray) -> np.ndarray:
        """Popovics/Mander confined curve, fib B14 Eq. (6-3)/(6-4)."""
        eps = np.asarray(eps, dtype=float)
        out = np.zeros_like(eps)
        m = eps > 0.0
        if not np.any(m):
            return out
        x = eps[m] / self.eps_cc
        Esec = self.fcc / self.eps_cc
        r = self.base.Ec / max(self.base.Ec - Esec, 1e-6)
        out[m] = self.fcc * x * r / (r - 1.0 + x ** r)
        out[m] = np.where(eps[m] > self.eps_cu, 0.0, out[m])
        return out


@dataclass
class Steel:
    """Elastic-perfectly-plastic reinforcement.

    ACI 318-19 22.4.2.1 limits f_y to 550 MPa for axial strength
    calculations; 20.2.2.4 governs elsewhere.
    """
    fy: float
    Es: float = 200000.0
    eps_su: float = 0.09
    label: str = ""

    @property
    def eps_ty(self) -> float:
        """eps_ty, ACI 318-19 21.2.2.1 = f_y/E_s."""
        return self.fy / self.Es

    def sigma(self, eps):
        eps = np.asarray(eps, dtype=float)
        return np.clip(self.Es * eps, -self.fy, self.fy)


@dataclass
class Bar:
    """One longitudinal bar.

    y     : distance from the EXTREME COMPRESSION FIBRE of the jacketed
            section, mm (i.e. from the outer face of the jacket, top).
    z     : transverse coordinate, mm from the section centreline. Only used
            for plotting and for bending about the orthogonal axis.
    area  : bar area, mm^2
    steel : Steel
    zone  : 'existing' or 'jacket' -- controls whether the locked-in
            pre-strain applies to this bar.
    """
    y: float
    z: float
    area: float
    steel: Steel
    zone: Literal["existing", "jacket"] = "existing"


@dataclass
class TieSet:
    """Transverse reinforcement of the jacket.

    Av_x  : total area of tie legs crossing a plane normal to x, mm^2
            (i.e. legs that resist shear in the x/depth direction)
    Av_z  : total area of tie legs in the orthogonal direction, mm^2
    s     : centre-to-centre spacing along the member, mm
    fyt   : tie yield strength, MPa
    db    : tie bar diameter, mm
    """
    Av_x: float
    Av_z: float
    s: float
    fyt: float
    db: float = 10.0


def rect_perimeter_bars(b: float, h: float, cover_to_bar_centre: float,
                        n_b: int, n_h: int, dia: float, steel: Steel,
                        zone: str = "existing",
                        y_offset: float = 0.0,
                        z_offset: float = 0.0,
                        faces: tuple = (True, True, True, True)) -> List[Bar]:
    """Generate bars on the perimeter of a b (width, z) x h (depth, y)
    rectangle.  n_b bars along each face of width b, n_h bars along each
    face of depth h (corner bars counted once).

    y is measured downward from the top of THIS rectangle, then shifted by
    y_offset so it becomes a distance from the jacketed section's extreme
    compression fibre.
    """
    area = math.pi * dia ** 2 / 4.0
    c = cover_to_bar_centre
    ys = np.linspace(c, h - c, max(n_h, 2))
    zs = np.linspace(-(b / 2 - c), (b / 2 - c), max(n_b, 2))
    f_top, f_bot, f_left, f_right = faces
    coords = set()
    for z in zs:                       # top and bottom faces
        if f_top:
            coords.add((round(ys[0], 4), round(z, 4)))
        if f_bot:
            coords.add((round(ys[-1], 4), round(z, 4)))
    for y in ys:                       # left and right faces
        if f_left:
            coords.add((round(y, 4), round(zs[0], 4)))
        if f_right:
            coords.add((round(y, 4), round(zs[-1], 4)))
    return [Bar(y=y + y_offset, z=z + z_offset, area=area, steel=steel,
                zone=zone)
            for (y, z) in sorted(coords)]


# --------------------------------------------------------------------------
# 2.  THE JACKETED SECTION
# --------------------------------------------------------------------------


@dataclass
class JacketedColumn:
    """Existing rectangular column b_e (width) x h_e (depth) with a
    four-sided jacket of thickness t.

    Bending is taken about the axis normal to the depth direction, i.e.
    the neutral axis is horizontal and the section depth is H = h_e + 2t.

    Attributes
    ----------
    b_e, h_e : existing column width and depth, mm
    t        : jacket thickness (uniform, all four faces), mm
    conc_e   : existing concrete
    conc_j   : jacket concrete
    bars_e   : existing longitudinal bars (y measured from jacket outer face)
    bars_j   : new jacket longitudinal bars
    ties_j   : jacket transverse reinforcement
    cover_j  : clear cover from jacket face to jacket tie, mm
    spiral   : True if the jacket transverse steel is a spiral conforming to
               ACI 318-19 25.7.3 (raises P_n,max from 0.80P_o to 0.85P_o and
               phi from 0.65 to 0.75)
    """
    b_e: float
    h_e: float
    t: float
    conc_e: Concrete
    conc_j: Concrete
    bars_e: List[Bar]
    bars_j: List[Bar]
    ties_j: TieSet
    cover_j: float = 40.0
    spiral: bool = False
    ke_conf: Optional[float] = None      # None -> compute from geometry
    n_sup_b: int = 2                     # laterally supported jacket bars per
    n_sup_h: int = 2                     # face (corners only = 2; every bar
                                         # supported requires crossties)
    db_long: float = 20.0

    # ---- partial jacket geometry (CODE_IMPLEMENTATION_PARTIAL_JACKET s2.1) --
    # None -> t.  0.0 -> that face is NOT jacketed.  `t` stays field 3 and
    # every existing call site builds exactly the section it built before.
    t_top: Optional[float] = None
    t_bot: Optional[float] = None
    t_left: Optional[float] = None
    t_right: Optional[float] = None

    # ---- transverse steel continuity --------------------------------------
    hoop_closed: bool = True             # ties form a closed loop round the core
    s_closed: Optional[float] = None     # spacing of CLOSED/anchored ties; None -> ties_j.s
    ties_anchored: bool = True           # tie legs develop f_yt inside the core

    _bare: bool = False                  # set only by bare_column(); see s4.2(c)

    def __post_init__(self):
        for f in ("t_top", "t_bot", "t_left", "t_right"):
            if getattr(self, f) is None:
                setattr(self, f, self.t)
        if self.s_closed is None:
            self.s_closed = self.ties_j.s
        if self.n_faces == 0 and not self._bare:
            raise ValueError(
                "no face jacketed; use JacketedColumn.bare_column() for the "
                "un-jacketed section"
            )

    @classmethod
    def bare_column(cls, **kw) -> "JacketedColumn":
        """The un-jacketed section, for the existing_only path of s4.2(c).

        A dedicated constructor rather than a public flag: n_faces == 0 is a
        real error everywhere else, and a keyword that switches the guard off
        would eventually be passed by something that meant a partial jacket.
        """
        kw.update(t=0.0, t_top=0.0, t_bot=0.0, t_left=0.0, t_right=0.0)
        obj = cls(**kw, _bare=True)
        return obj


    # ---------------- geometry ------------------------------------------
    @property
    def faces(self) -> tuple:
        """(top, bottom, left, right) -- which faces carry jacket concrete."""
        return (self.t_top > 0, self.t_bot > 0,
                self.t_left > 0, self.t_right > 0)

    @property
    def n_faces(self) -> int:
        return sum(self.faces)

    @property
    def B(self) -> float:
        """Overall width of the jacketed section, mm."""
        return self.b_e + self.t_left + self.t_right

    @property
    def H(self) -> float:
        """Overall depth of the jacketed section, mm."""
        return self.h_e + self.t_top + self.t_bot

    @property
    def y_core(self) -> tuple:
        """(y at the top of the existing core, y at its soffit), mm."""
        return (self.t_top, self.t_top + self.h_e)

    @property
    def z_core(self) -> tuple:
        z0 = -self.B / 2 + self.t_left
        return (z0, z0 + self.b_e)

    @property
    def yc_e(self) -> float:
        """Centroid of the existing core, measured from the extreme fibre."""
        return self.t_top + self.h_e / 2

    @property
    def zc_e(self) -> float:
        """Centroid of the existing core about the COMPOSITE centreline.

        Not zero for a partial jacket -- that is why bars_e need a z_offset.
        """
        return self.z_core[0] + self.b_e / 2

    @property
    def yc_j(self) -> float:
        """Centroid of the jacket region (the composite less the core)."""
        return (self.Ag_tot * self.H / 2 - self.Ag_e * self.yc_e) / self.Ag_j

    @property
    def zc_j(self) -> float:
        return (0.0 - self.Ag_e * self.zc_e) / self.Ag_j

    @property
    def Ag_e(self) -> float:
        return self.b_e * self.h_e

    @property
    def Ag_tot(self) -> float:
        return self.B * self.H

    @property
    def Ag_j(self) -> float:
        """Gross area of the jacket annulus, mm^2."""
        return self.Ag_tot - self.Ag_e

    @property
    def Ast_e(self) -> float:
        return sum(bb.area for bb in self.bars_e)

    @property
    def Ast_j(self) -> float:
        return sum(bb.area for bb in self.bars_j)

    @property
    def Ast_tot(self) -> float:
        return self.Ast_e + self.Ast_j

    @property
    def bars(self) -> List[Bar]:
        return list(self.bars_e) + list(self.bars_j)

    @property
    def Ig_e(self) -> float:
        return self.b_e * self.h_e ** 3 / 12.0

    @property
    def Ig_tot(self) -> float:
        return self.B * self.H ** 3 / 12.0

    def widths(self, y: np.ndarray):
        """At depth y from the extreme compression fibre, return
        (width of jacket concrete, width of existing-core concrete), mm."""
        y = np.asarray(y, dtype=float)
        y0, y1 = self.y_core
        inside = (y >= y0) & (y <= y1)
        w_core = np.where(inside, self.b_e, 0.0)
        w_jack = np.where(inside, self.B - self.b_e, self.B)
        w_jack = np.where((y < 0) | (y > self.H), 0.0, w_jack)
        w_core = np.where((y < 0) | (y > self.H), 0.0, w_core)
        return w_jack, w_core

    # ---------------- confinement ---------------------------------------
    @property
    def core_dims(self):
        """Centre-to-centre dimensions of the jacket tie cage, mm."""
        self._require_four_sided("core_dims (the jacket tie cage)", "phase 5 (spec s4.5)")
        d = self.cover_j + self.ties_j.db / 2.0
        return self.B - 2 * d, self.H - 2 * d      # b_c, h_c

    @property
    def rho_s_ties(self) -> float:
        """Volumetric ratio of jacket ties, standard rectangular form
            rho_s = (A_sx*h_c + A_sz*b_c) / (s * b_c * h_c)
        """
        self._require_four_sided("rho_s_ties", "phase 5 (spec s4.5)")
        b_c, h_c = self.core_dims
        return (self.ties_j.Av_x * h_c + self.ties_j.Av_z * b_c) / \
               (self.ties_j.s * b_c * h_c)

    def confined(self, which: Literal["existing", "jacket"] = "existing") -> ConfinedConcrete:
        """Confined properties of the concrete inside the jacket tie cage.
        The jacket ties confine BOTH the thin inner ring of jacket concrete
        and the whole existing core.

        NOTE on which f'co to use.  f'cc scales with f'co, so taking the
        weaker (existing) concrete is conservative for STRENGTH -- but
        eps_cu = 0.004 + 1.4*rho_s*f_yh*eps_su/f'cc is INVERSELY proportional
        to f'cc, so the same choice INFLATES the ultimate strain.  Using the
        minimum is not uniformly conservative.  Check ductility-governed
        results both ways.
        """
        self._require_four_sided("confined", "phase 5 (spec s4.5)")
        base = self.conc_e if which == "existing" else self.conc_j
        return ConfinedConcrete(base=base, rho_st=self.rho_s_ties,
                                fyh=self.ties_j.fyt, ke=self.mander_ke())

    def mander_ke(self) -> float:
        """Confinement effectiveness coefficient k_e computed from the ACTUAL
        tie and bar geometry rather than assumed.  Mander et al. (1988),
        rectangular sections:

            k_e = [1 - sum(w_i^2)/(6*b_c*h_c)]
                  * [1 - s'/(2 b_c)] * [1 - s'/(2 h_c)] / (1 - rho_cc)

        w_i    = clear distance between adjacent LATERALLY SUPPORTED
                 longitudinal bars.  A bar is supported only where a tie or
                 crosstie corner engages it (ACI 318-19 25.7.2.3).
        s'     = clear vertical spacing between ties = s - db_tie
        rho_cc = longitudinal steel ratio referred to the core area

        This matters more than any other single input.  A plain perimeter
        hoop with no crossties supports only the four corner bars, so w_i is
        the full face width and k_e collapses to roughly 0.35 -- less than
        half the 0.8 that is routinely assumed.  Assuming 0.8 without
        crossties overstates the confined strength of the existing core and
        is unconservative.  Set n_sup_b / n_sup_h honestly, or pass ke_conf
        explicitly to override with a documented justification.
        """
        self._require_four_sided("mander_ke", "phase 5 (spec s4.5)")
        if self.ke_conf is not None:
            return self.ke_conf
        b_c, h_c = self.core_dims
        d_bar = self.cover_j + self.ties_j.db + self.db_long / 2.0
        span_b = self.B - 2 * d_bar          # c/c between corner bars
        span_h = self.H - 2 * d_bar
        nb, nh = max(self.n_sup_b, 2), max(self.n_sup_h, 2)
        w_b = span_b / (nb - 1) - self.db_long
        w_h = span_h / (nh - 1) - self.db_long
        sum_w2 = 2 * (nb - 1) * max(w_b, 0.0) ** 2 + \
                 2 * (nh - 1) * max(w_h, 0.0) ** 2
        s_clear = max(self.ties_j.s - self.ties_j.db, 0.0)
        rho_cc = self.Ast_tot / (b_c * h_c)
        ke = (1.0 - sum_w2 / (6.0 * b_c * h_c)) \
            * (1.0 - s_clear / (2.0 * b_c)) \
            * (1.0 - s_clear / (2.0 * h_c)) / (1.0 - rho_cc)
        return max(min(ke, 1.0), 0.0)

    # ---------------- axial ---------------------------------------------
    def Po(self) -> float:
        """Nominal axial strength at zero eccentricity, ACI 318-19
        Eq. (22.4.2.2), applied to each concrete and steel population:

            P_o = 0.85 f'c_e (A_g,e - A_st,e) + f_y,e A_st,e
                + 0.85 f'c_j (A_g,j - A_st,j) + f_y,j A_st,j

        f_y is capped at 550 MPa per 22.4.2.1.
        """
        def fy_cap(bars):
            return [(bb.area, min(bb.steel.fy, 550.0)) for bb in bars]

        Pe = 0.85 * self.conc_e.fc * (self.Ag_e - self.Ast_e) + \
            sum(a * f for a, f in fy_cap(self.bars_e))
        Pj = 0.85 * self.conc_j.fc * (self.Ag_j - self.Ast_j) + \
            sum(a * f for a, f in fy_cap(self.bars_j))
        return Pe + Pj

    def plastic_centroid(self) -> tuple:
        """Line of action of P_o: the point about which the interaction diagram
        of an UNSYMMETRIC section must be taken.

            y_pc = Sum(F_i y_i) / Sum(F_i)   over the same contributions as Po:
                0.85 f'c,e (Ag_e - Ast_e)  at  yc_e
                0.85 f'c,j (Ag_j - Ast_j)  at  yc_j
                f_y A_b                    at  each bar, f_y capped at 550 MPa
                                               (ACI 318-19 22.4.2.1)

        INTERPRETATION, no clause.  ACI 318 is understood in practice to require
        moments on a nonsymmetrical section to be referenced to the plastic
        centroid, but the clause is not verified in this workspace.
        See TN-RET-002 TODO B7(a) / C1 -- BLOCKING.

        For a four-sided symmetric jacket this returns exactly (H/2, 0.0),
        which is why the pre-existing H/2 convention is correct there.
        """
        F, Fy, Fz = 0.0, 0.0, 0.0

        Fe = 0.85 * self.conc_e.fc * (self.Ag_e - self.Ast_e)
        F += Fe
        Fy += Fe * self.yc_e
        Fz += Fe * self.zc_e

        if self.Ag_j > 0:
            Fj = 0.85 * self.conc_j.fc * (self.Ag_j - self.Ast_j)
            F += Fj
            Fy += Fj * self.yc_j
            Fz += Fj * self.zc_j

        for bb in self.bars:
            Fb = bb.area * min(bb.steel.fy, 550.0)
            F += Fb
            Fy += Fb * bb.y
            Fz += Fb * bb.z

        if F <= 0:
            return (self.H / 2.0, 0.0)
        return (Fy / F, Fz / F)

    def _require_four_sided(self, what: str, phase: str) -> None:
        """Refuse a method that still assumes a jacket closed on four faces.

        A partial jacket can now be BUILT (phases 1-3), but the confinement,
        shear, interface, stiffness and detailing routines have not been
        migrated yet. Left ungated they answer as though the jacket wrapped the
        column -- core_dims returns a tie cage on all four faces, mander_ke
        returns a confinement credit for a hoop that cannot close. That is the
        silent extrapolation this workspace exists to prevent, so these refuse
        until their phase lands.

        Not an open question being resolved: it is one being declined.
        """
        ts = (self.t_top, self.t_bot, self.t_left, self.t_right)
        uniform = max(ts) - min(ts) <= 1e-9
        if self.n_faces < 4 or not uniform:
            detail = (f"has {self.n_faces}" if self.n_faces < 4 else
                      f"has four faces of unequal thickness {ts}")
            raise NotImplementedError(
                f"{what} still assumes a jacket of UNIFORM thickness closed on "
                f"all four faces, and this section {detail}. Implementing it "
                f"for a "
                f"partial jacket is {phase} of "
                f"CODE_IMPLEMENTATION_PARTIAL_JACKET.md, which has not landed. "
                f"Answering now would silently apply four-sided behaviour to a "
                f"partial jacket."
            )

    @property
    def split_compression_fibre(self) -> bool:
        """True when BOTH concretes sit at the extreme compression fibre.

        The single-block rule takes one depth a = beta1*c with beta1 from the
        concrete at the extreme fibre.  That presumes ONE concrete is there.
        With side jackets and no top jacket the extreme fibre carries jacket
        strips and the core side by side, and no beta1 is 'the' one: picking
        the core's (larger) beta1 hands the STRONGER jacket concrete a deeper
        block, which is the inversion CLAUDE.md records as physically wrong and
        unconservative -- measurably so, the default returning up to 3.5% MORE
        than 'per_material', the rule kept only for comparison.

        Which beta1 governs a split fibre is not something this module may
        decide.  CODE_IMPLEMENTATION_PARTIAL_JACKET s4.2(b) prescribes
        `t_top > 0`, and that prescription is only valid where one concrete is
        at the fibre; for a split fibre the specification is silent.  Recorded
        as a specification defect, not resolved here.
        """
        return self.t_top <= 0 and (self.t_left > 0 or self.t_right > 0)

    @property
    def doubly_unsymmetric(self) -> bool:
        """True when the section is unsymmetric about BOTH principal axes.

        A uniaxial interaction diagram is only meaningful when the section is
        symmetric about the plane of bending. J2-L (one long face plus one
        short face) is unsymmetric both ways, so its capacity is a biaxial
        SURFACE, not a curve, and this module does not construct one.

        TN-RET-002 TODO C2 -- BLOCKING for J2-L. The tool refuses rather than
        returning the uniaxial curve, which would silently ignore M_z.
        """
        return (abs(self.t_top - self.t_bot) > 1e-9
                and abs(self.t_left - self.t_right) > 1e-9)

    def reversed(self) -> "JacketedColumn":
        """The same physical column with the bending sense reversed: t_top and
        t_bot swapped, and every bar's y mirrored about H/2.

        Required because an unsymmetric section has TWO interaction diagrams
        and seismic moments reverse. Returns a new object; never mutates self.
        H is unchanged by the swap, so the mirror is exact.

        For a section symmetric about the bending axis the result's
        interaction_code() output is identical to self's.
        """
        H = self.H

        def mirror(bars):
            return [Bar(H - bb.y, bb.z, bb.area, bb.steel, bb.zone)
                    for bb in bars]

        return JacketedColumn(
            b_e=self.b_e, h_e=self.h_e, t=self.t,
            conc_e=self.conc_e, conc_j=self.conc_j,
            bars_e=mirror(self.bars_e), bars_j=mirror(self.bars_j),
            ties_j=self.ties_j, cover_j=self.cover_j, spiral=self.spiral,
            ke_conf=self.ke_conf, n_sup_b=self.n_sup_b, n_sup_h=self.n_sup_h,
            db_long=self.db_long,
            t_top=self.t_bot, t_bot=self.t_top,
            t_left=self.t_left, t_right=self.t_right,
            hoop_closed=self.hoop_closed, s_closed=self.s_closed,
            ties_anchored=self.ties_anchored)

    def induced_eccentricity(self, P: float) -> dict:
        """Moment produced because the axial load arrives through the EXISTING
        column's centroid and is resisted about the composite plastic centroid.

            e_y = y_pc - yc_e        M_y = P * |e_y|
            e_z = z_pc - zc_e        M_z = P * |e_z|

        ADDITIVE in both bending senses -- the sign of the frame moment
        reverses under seismic demand and the eccentricity does not.
        Mechanics, no clause; TN-RET-002 s2, marked [derivation] there.

        e_geom_* are the GEOMETRIC offsets. They are reported because they are
        what an engineer estimates by eye, and they UNDERSTATE the real
        eccentricity, because the stronger new materials pull the plastic
        centroid further. Do not use them in the check.
        """
        y_pc, z_pc = self.plastic_centroid()
        e_y = y_pc - self.yc_e
        e_z = z_pc - self.zc_e
        return dict(e_y=e_y, e_z=e_z,
                    M_y=abs(P * e_y), M_z=abs(P * e_z),
                    e_geom_y=self.H / 2.0 - self.yc_e,
                    e_geom_z=0.0 - self.zc_e)

    def Po_existing_only(self) -> float:
        """P_o of the bare existing column, for the strength-gain ratio."""
        return 0.85 * self.conc_e.fc * (self.Ag_e - self.Ast_e) + \
            sum(bb.area * min(bb.steel.fy, 550.0) for bb in self.bars_e)

    @property
    def alpha_max(self) -> float:
        """0.80 for ties, 0.85 for spirals, ACI 318-19 Table 22.4.2.1."""
        return 0.85 if self.spiral else 0.80

    @property
    def phi_comp(self) -> float:
        """phi for compression-controlled sections, ACI 318-19
        Table 21.2.2 -- 0.75 spiral, 0.65 other."""
        return 0.75 if self.spiral else 0.65

    def phiPn_max(self, phi_override: Optional[float] = None) -> float:
        """phi * P_n,max, ACI 318-19 22.4.2.1 with Table 21.2.2."""
        phi = self.phi_comp if phi_override is None else phi_override
        return phi * self.alpha_max * self.Po()

    def Pnt_max(self) -> float:
        """Maximum axial tensile strength, ACI 318-19 Eq. (22.4.3.1),
        non-prestressed:  P_nt,max = f_y * A_st."""
        return sum(bb.area * bb.steel.fy for bb in self.bars)

    # ---------------- locked-in pre-strain ------------------------------
    def preload_strain(self, P0: float, creep_factor: float = 0.0) -> float:
        """Uniform compressive strain locked into the EXISTING column at the
        moment the jacket is cast, under sustained axial load P0 (N).

            eps_0 = P0 / (E_c,e * A_c,e + E_s * A_st,e)  * (1 + phi_creep)

        creep_factor is the creep coefficient phi_cc applied to the sustained
        part.  0.0 = short-term elastic only.  A value of 1.5-2.5 is typical
        for a long-standing column; it can double or triple eps_0 and is the
        realistic case.  This term is NOT a code equation -- it is mechanics,
        and it is why an unshored jacket does not deliver the monolithic
        capacity.  See TN-RET-001.
        """
        if P0 <= 0:
            return 0.0
        EA = self.conc_e.Ec * (self.Ag_e - self.Ast_e) + \
            sum(bb.steel.Es * bb.area for bb in self.bars_e)
        return P0 / EA * (1.0 + creep_factor)

    # ---------------- P-M interaction, ACI stress block -----------------
    def interaction_code(self, n_layers: int = 400, n_points: int = 60,
                         existing_only: bool = False,
                         stress_block: str = "single",
                         sense: str = "as_built"):
        """P-M interaction diagram using the ACI 318-19 equivalent
        rectangular stress block.

        Assumptions and interpretation:
          * eps_cu = 0.003 at the extreme compression fibre (22.2.2.1)
          * concrete stress 0.85 f'c over a depth a = beta1*c (22.2.2.4.1).
            ACI 318 does not address a section with two concrete strengths,
            so a rule has to be chosen.  A SINGLE block depth is used, with
            beta1 taken from the concrete AT THE EXTREME COMPRESSION FIBRE
            (the jacket), and each material contributing its own 0.85 f'c
            within that depth.
            The alternative -- giving each concrete its own beta1*c measured
            from y = 0 -- is physically wrong and is deliberately NOT used:
            beta1 and the 0.85 factor are jointly calibrated to the strain
            profile measured from the extreme fibre, and the existing core
            does not reach y = 0.  With a 45 MPa jacket on a 21 MPa core it
            would give the weaker interior core a DEEPER block (0.85c) than
            the stronger jacket at the extreme fibre (0.729c) -- inverted,
            and unconservative.  Set stress_block='per_material' to reproduce
            that behaviour for comparison only.
            Set the two f'c equal to recover the monolithic result exactly.
          * concrete tensile strength neglected (22.2.2.2)
          * bar area displaced by the stress block is deducted
          * moments taken about the PLASTIC CENTROID, which for a symmetric
            four-sided jacket is exactly H/2 but for a partial jacket is not.
            Taking them about the geometric mid-depth is TN-RET-002 Common
            Error 2; this docstring described that error for a while after the
            code stopped making it.
          * phi from Table 21.2.2 using eps_t at the extreme tension bar.
          * PRELOAD IS IGNORED -- this is the monolithic-equivalent code
            check.  Use interaction_fiber() to see the preload penalty.

        Returns dict with arrays P, M, phiP, phiM, eps_t, phi (N and N*mm).
        """
        # s4.2(d): an unsymmetric section has TWO diagrams. "reversed" mirrors
        # the section and re-enters as_built, so there is no second code path
        # to keep in step and no recursion beyond one level.
        if sense not in ("as_built", "reversed"):
            raise ValueError(f"sense must be as_built or reversed, got {sense!r}")
        if sense == "reversed":
            return self.reversed().interaction_code(
                n_layers, n_points, existing_only, stress_block,
                sense="as_built")

        if self.split_compression_fibre and not existing_only:
            raise NotImplementedError(
                "The extreme compression fibre of this section carries BOTH "
                f"concretes side by side (t_top={self.t_top:g} with "
                f"t_left={self.t_left:g}, t_right={self.t_right:g}), so the "
                "single stress block has no single beta1. Taking the core's "
                "beta1 would apply a deeper block to the stronger jacket "
                "concrete -- the inversion CLAUDE.md records as unconservative, "
                "and it makes the default exceed 'per_material'. "
                "CODE_IMPLEMENTATION_PARTIAL_JACKET s4.2(b) prescribes "
                "t_top > 0, which only decides the case where one concrete is "
                "at the fibre; for a split fibre the specification is silent. "
                "Refused rather than guessed -- report it as a spec defect."
            )

        if self.doubly_unsymmetric and not existing_only:
            raise NotImplementedError(
                "This section is unsymmetric about both principal axes "
                f"(t_top={self.t_top:g}, t_bot={self.t_bot:g}, "
                f"t_left={self.t_left:g}, t_right={self.t_right:g}), so its "
                "capacity is a biaxial SURFACE, not a uniaxial curve. How that "
                "surface is constructed is unresolved -- TN-RET-002 TODO C2, "
                "BLOCKING. Returning the uniaxial curve would silently ignore "
                "M_z, so this refuses instead of approximating."
            )

        # An unsymmetric section has two diagrams and the weaker governs. If
        # the OTHER sense cannot be produced, returning this one alone is the
        # single-sense error the sense= parameter exists to prevent -- it just
        # arrives wearing a refusal. Refuse both.
        if (sense == "as_built" and not existing_only
                and abs(self.t_top - self.t_bot) > 1e-9
                and self.reversed().split_compression_fibre):
            raise NotImplementedError(
                "This section is unsymmetric top-to-bottom, so BOTH bending "
                "senses are required and the weaker governs. The reversed "
                "sense cannot be produced: mirroring puts both concretes at the "
                "extreme compression fibre, where the single stress block has "
                "no single beta1 (see split_compression_fibre). Returning this "
                "sense alone would report the section as checked when only one "
                "of its two diagrams exists, so both are refused."
            )

        if existing_only:
            # s4.2(c): the bare section has no jacketed face, so it must be
            # built through bare_column() -- the n_faces == 0 guard is a real
            # error everywhere else.
            sub = JacketedColumn.bare_column(
                b_e=self.b_e, h_e=self.h_e,
                conc_e=self.conc_e, conc_j=self.conc_e,
                bars_e=self.bars_e, bars_j=[], ties_j=self.ties_j,
                cover_j=self.cover_j, spiral=self.spiral)
            # shift bars so y is measured from the bare column's own face.
            # s4.2(c): this was self.t, which is wrong the moment the top and
            # bottom thicknesses differ.
            shift = self.t_top
            sub.bars_e = [Bar(bb.y - shift, bb.z, bb.area, bb.steel, bb.zone)
                          for bb in self.bars_e]
            return sub.interaction_code(n_layers, n_points,
                                        stress_block=stress_block)

        H, B = self.H, self.B
        y_edges = np.linspace(0.0, H, n_layers + 1)
        y_mid = 0.5 * (y_edges[:-1] + y_edges[1:])
        dy = y_edges[1] - y_edges[0]
        w_j, w_e = self.widths(y_mid)
        # s4.2(a): moments on an unsymmetric section are referenced to the
        # plastic centroid, not the geometric mid-depth. Identical for a
        # four-sided symmetric jacket (s9.4), so the filed values cannot move.
        yc, _ = self.plastic_centroid()

        bars = self.bars
        by = np.array([bb.y for bb in bars])
        ba = np.array([bb.area for bb in bars])
        bfy = np.array([bb.steel.fy for bb in bars])
        bEs = np.array([bb.steel.Es for bb in bars])
        # concrete displaced by each bar
        bfc = np.array([(self.conc_e.fc if bb.zone == "existing"
                         else self.conc_j.fc) for bb in bars])
        if stress_block == "per_material":
            b1_e, b1_j = self.conc_e.beta1, self.conc_j.beta1
        else:
            # single block, beta1 from the concrete at the extreme fibre
            # s4.2(b): which concrete sits at the extreme compression fibre
            # now depends on the configuration and the bending sense. For
            # J1-c reversed, or J3-U as built, it is the existing core.
            fc_top_is_jacket = self.t_top > 0
            b1_e = b1_j = (self.conc_j.beta1 if fc_top_is_jacket
                           else self.conc_e.beta1)
        eps_cu = 0.003

        c_list = np.concatenate([
            np.linspace(0.02 * H, 1.2 * H, n_points - 6),
            np.array([1.5 * H, 2.0 * H, 3.0 * H, 6.0 * H, 20.0 * H, 1e6])])

        P, M, EPST, PHI = [], [], [], []
        for c in c_list:
            a_j, a_e = b1_j * c, b1_e * c
            # concrete compression
            sj = np.where(y_mid <= a_j, 0.85 * self.conc_j.fc, 0.0)
            se = np.where(y_mid <= a_e, 0.85 * self.conc_e.fc, 0.0)
            Fc = np.sum(sj * w_j * dy) + np.sum(se * w_e * dy)
            Mc = np.sum(sj * w_j * dy * (yc - y_mid)) + \
                np.sum(se * w_e * dy * (yc - y_mid))
            # steel
            eps_s = eps_cu * (c - by) / c            # + = compression
            fs = np.clip(bEs * eps_s, -bfy, bfy)
            a_mat = np.where(np.array([bb.zone == "existing" for bb in bars]),
                             a_e, a_j)
            fs_net = np.where(by <= a_mat, fs - 0.85 * bfc, fs)
            Fs = np.sum(fs_net * ba)
            Ms = np.sum(fs_net * ba * (yc - by))

            Pn, Mn = Fc + Fs, Mc + Ms
            eps_t = float(np.min(eps_s))             # most tensile
            eps_ty = float(bfy[int(np.argmin(eps_s))] /
                           bEs[int(np.argmin(eps_s))])
            et = -eps_t                              # tension positive
            if et <= eps_ty:
                phi = self.phi_comp
            elif et >= eps_ty + 0.003:
                phi = 0.90
            else:
                lo = self.phi_comp
                phi = lo + (0.90 - lo) * (et - eps_ty) / 0.003
            P.append(Pn); M.append(Mn); EPST.append(et); PHI.append(phi)

        P = np.array(P); M = np.array(M)
        PHI = np.array(PHI); EPST = np.array(EPST)
        phiP = np.minimum(PHI * P, self.phiPn_max())
        phiM = PHI * M
        # pure tension terminal point
        Pt = -self.Pnt_max()
        P = np.append(P, Pt); M = np.append(M, 0.0)
        phiP = np.append(phiP, 0.90 * Pt); phiM = np.append(phiM, 0.0)
        PHI = np.append(PHI, 0.90); EPST = np.append(EPST, 0.05)
        return dict(P=P, M=M, phiP=phiP, phiM=phiM, phi=PHI, eps_t=EPST,
                    c=np.append(c_list, 0.0))

    # ---------------- P-M interaction, fibre model with preload ---------
    def interaction_fiber(self, P0: float = 0.0, creep_factor: float = 0.0,
                          n_layers: int = 400, n_points: int = 50,
                          confine_core: bool = True):
        """Preload-aware fibre-model failure surface.

        Physics
        -------
        When the jacket is cast on a loaded column, the existing core already
        carries a locked-in compressive strain eps_0.  The jacket starts at
        zero strain.  Any further load produces an INCREMENT of strain
        d_eps(y) that is shared, so at any instant

            existing core fibres :  eps = eps_0 + d_eps(y)
            jacket fibres        :  eps =         d_eps(y)

        The section is exhausted when EITHER
            (a) the jacket's extreme compression fibre reaches eps_cu,unconf
            (b) the existing core's most compressed fibre reaches its
                CONFINED eps_cu, having started eps_0 ahead
        For a heavily preloaded column (b) governs -- which is exactly why
        the jacket ties matter: they raise eps_cu of the core and buy back
        the strain the preload consumed.

        Set P0 = 0 (or shore the column) to recover the fully composite case.

        Returns dict with arrays P, M (N, N*mm), plus kappa and the governing
        criterion at each point.
        """
        self._require_four_sided("interaction_fiber", "phase 7 (spec s4.3)")
        eps0 = self.preload_strain(P0, creep_factor)
        H = self.H
        y_edges = np.linspace(0.0, H, n_layers + 1)
        y_mid = 0.5 * (y_edges[:-1] + y_edges[1:])
        dy = y_edges[1] - y_edges[0]
        w_j, w_e = self.widths(y_mid)
        yc = H / 2.0
        d_tie = self.cover_j + self.ties_j.db / 2.0

        conf_e = self.confined("existing")
        conf_j = self.confined("jacket")
        eps_cu_j_unconf = 0.0038
        eps_cu_core = conf_e.eps_cu if confine_core else 0.0038

        # which layers are inside the jacket tie cage (confined)
        confined_zone = (y_mid >= d_tie) & (y_mid <= H - d_tie)

        bars = self.bars
        by = np.array([bb.y for bb in bars])
        ba = np.array([bb.area for bb in bars])
        is_e = np.array([bb.zone == "existing" for bb in bars])
        steels = [bb.steel for bb in bars]

        y_core_top = self.t            # first existing-core fibre

        kappas = np.concatenate([[0.0],
                                 np.geomspace(1e-7, 2.5e-4, n_points - 1)])
        P, M, GOV = [], [], []
        for kap in kappas:
            # candidate top strains from the two crush criteria
            e_top_a = eps_cu_j_unconf
            e_top_b = eps_cu_core - eps0 + kap * y_core_top
            if e_top_b < e_top_a:
                e_top, gov = e_top_b, "existing core crushes"
            else:
                e_top, gov = e_top_a, "jacket face crushes"
            e_top = max(e_top, 1e-6)

            d_eps = e_top - kap * y_mid                 # incremental field
            eps_j = d_eps
            eps_e = d_eps + eps0

            # jacket concrete: confined inside the tie cage, else unconfined
            s_j = np.where(confined_zone,
                           conf_j.sigma(eps_j), self.conc_j.sigma_unconfined(eps_j))
            s_e = conf_e.sigma(eps_e) if confine_core else \
                self.conc_e.sigma_unconfined(eps_e)
            Fc = np.sum(s_j * w_j * dy) + np.sum(s_e * w_e * dy)
            Mc = np.sum(s_j * w_j * dy * (yc - y_mid)) + \
                np.sum(s_e * w_e * dy * (yc - y_mid))

            eps_b = np.where(is_e, e_top - kap * by + eps0, e_top - kap * by)
            fs = np.array([st.sigma(e) for st, e in zip(steels, eps_b)])
            fcd = np.where(is_e, conf_e.sigma(np.maximum(eps_b, 0)),
                           np.where(by >= d_tie,
                                    conf_j.sigma(np.maximum(eps_b, 0)),
                                    self.conc_j.sigma_unconfined(np.maximum(eps_b, 0))))
            fs_net = np.where(eps_b > 0, fs - fcd, fs)
            Fs = np.sum(fs_net * ba)
            Ms = np.sum(fs_net * ba * (yc - by))

            P.append(Fc + Fs); M.append(Mc + Ms); GOV.append(gov)

        return dict(P=np.array(P), M=np.array(M), kappa=kappas,
                    governing=GOV, eps0=eps0,
                    fcc_core=conf_e.fcc, eps_cu_core=conf_e.eps_cu,
                    fl=conf_e.fl, rho_s=self.rho_s_ties)

    # ---------------- shear ---------------------------------------------
    def shear_capacity(self, Nu: float, fc_mode: str = "min",
                       d: Optional[float] = None, phi: float = 0.75):
        """One-way shear strength of the JACKETED section,
        ACI 318-19 22.5, SI (App. C).

            V_n = V_c + V_s                                     (22.5.1.1)
            V_c = (0.17*lam*sqrt(f'c) + N_u/(6 A_g)) b_w d      Table 22.5.5.1(a)
                  with N_u/(6A_g) <= 0.05 f'c                   (22.5.5.1.2)
                  and V_c <= 0.42*lam*sqrt(f'c) b_w d           (22.5.5.1.1)
            V_s = A_v f_yt d / s                                (22.5.8.5.3)
            V_u <= phi (V_c + 0.66 sqrt(f'c) b_w d)             (22.5.1.2)

        N_u positive in compression.  Only the JACKET ties are counted in
        V_s: the existing ties are usually the reason the retrofit is needed,
        they are lapped with 90-degree hooks in the cover that is about to be
        roughened off, and they cannot be relied on.  Add them deliberately
        if the condition assessment justifies it.

        fc_mode: 'min' (default, conservative) uses min(f'c_e, f'c_j);
                 'jacket' uses the jacket concrete;
                 'weighted' area-weights the two.
        """
        self._require_four_sided("shear_capacity", "phase 5 (spec s4.4)")
        if fc_mode == "min":
            fc = min(self.conc_e.fc, self.conc_j.fc)
        elif fc_mode == "jacket":
            fc = self.conc_j.fc
        else:
            fc = (self.conc_e.fc * self.Ag_e + self.conc_j.fc * self.Ag_j) / self.Ag_tot
        lam = min(self.conc_e.lam, self.conc_j.lam)
        bw = self.B
        if d is None:
            # d from the ACTUAL extreme tension bar, so the shear check and
            # the interaction diagram agree about where the steel is.  A
            # guessed cover stack (cover + tie + half bar) can differ from the
            # bar layout by 10-15 mm and silently desynchronise the two.
            d = max((bb.y for bb in self.bars), default=self.H * 0.9)
        Ag = self.Ag_tot

        axial = min(max(Nu, 0.0) / (6.0 * Ag), 0.05 * fc)
        Vc = (0.17 * lam * math.sqrt(fc) + axial) * bw * d
        Vc = min(Vc, 0.42 * lam * math.sqrt(fc) * bw * d)
        Vc = max(Vc, 0.0)
        Vs = self.ties_j.Av_x * self.ties_j.fyt * d / self.ties_j.s
        Vn_cap = Vc + 0.66 * math.sqrt(fc) * bw * d      # 22.5.1.2 ceiling
        Vn = min(Vc + Vs, Vn_cap)
        return dict(Vc=Vc, Vs=Vs, Vn=Vn, phiVn=phi * Vn, d=d, bw=bw,
                    fc_used=fc, capped=(Vc + Vs) > Vn_cap,
                    Vn_ceiling=Vn_cap, phi=phi)

    def shear_demand_from_hinging(self, Mn_top: float, Mn_bot: float,
                                  L_clear: float) -> float:
        """Shear associated with flexural hinging at both ends,
        V_p = (M_n,top + M_n,bot)/L_clear.  The jacketed column must not
        become shear-critical: a jacket that raises M_n without raising V_n
        converts a ductile failure into a brittle one.  ASCE 41-17
        Table 10-10a bands columns on V_yE/V_ColOE for exactly this reason.
        """
        return (Mn_top + Mn_bot) / L_clear

    # ---------------- interface ------------------------------------------
    def service_stress_split(self, P0: float, P_total: float,
                             creep_factor: float = 0.0):
        """Where the axial load actually sits after jacketing.

        A monolithic analysis assumes P_total is shared by core and jacket in
        proportion to axial stiffness.  It is not.  The sustained load P0 was
        locked into the core before the jacket existed and stays there; only
        the INCREMENT dP = P_total - P0 is shared.  The core therefore runs at
        a higher stress than a monolithic calculation reports, and the jacket
        at a lower one.  This is the practical consequence of an unshored
        jacket, and it is what limits how much of the jacket you may count on
        at service level.

        Shore the column (relieve P0) before casting and the two converge.
        """
        EA_c = self.conc_e.Ec / (1.0 + creep_factor) * (self.Ag_e - self.Ast_e) \
            + sum(bb.steel.Es * bb.area for bb in self.bars_e)
        EA_j = self.conc_j.Ec * (self.Ag_j - self.Ast_j) \
            + sum(bb.steel.Es * bb.area for bb in self.bars_j)
        dP = max(P_total - P0, 0.0)
        P_core = P0 + dP * EA_c / (EA_c + EA_j)
        P_jack = dP * EA_j / (EA_c + EA_j)
        # monolithic (fully composite from zero) comparison
        P_core_mono = P_total * EA_c / (EA_c + EA_j)
        Ac_e = self.Ag_e - self.Ast_e
        return dict(
            EA_core=EA_c, EA_jacket=EA_j,
            P_core=P_core, P_jacket=P_jack,
            P_core_monolithic=P_core_mono,
            core_overstress=P_core / P_core_mono if P_core_mono > 0 else 1.0,
            sigma_core=P_core / Ac_e,
            sigma_core_monolithic=P_core_mono / Ac_e,
            jacket_share=P_jack / P_total if P_total > 0 else 0.0,
            jacket_share_monolithic=EA_j / (EA_c + EA_j))

    def interface_check(self, C_jacket: float, dP_jacket: float,
                        L_shear_span: float, L_transfer: float,
                        rho_v: float, fy_dowel: float,
                        Vu: float = 0.0,
                        continuity: Literal["continuous", "discontinuous"]
                        = "discontinuous",
                        surface: str = "roughened_6mm",
                        phi: float = 0.75):
        """Longitudinal (interface) shear transfer between the existing
        column and the jacket.

        THIS IS WHAT MAKES THE JACKET STRUCTURAL.  A jacket that does not
        transfer force across the interface is formwork.

        The demand depends entirely on ONE detailing decision -- whether the
        jacket is continuous through the floor slab and bears directly on the
        foundation.  This is the single most consequential choice in jacket
        design and the one most often got wrong.

        continuity='continuous'
            The jacket passes through cored holes in the slab and bears on
            the footing, so axial load enters it by DIRECT BEARING and the
            jacket's own longitudinal bars are developed at both ends.  The
            interface then only has to carry the elastic compatibility shear
            flow:
                v_u = V_u * Q_j / (I_tr * b_v)          (VQ/Ib)
            with Q_j the first moment about the centroid of the jacket
            material outboard of the interface plane, and b_v = b_e, the
            contact width of the face being peeled off.

        continuity='discontinuous'
            The jacket stops at the slab soffit.  Nothing bears.  Every
            newton the jacket is credited with must migrate sideways across
            the interface:
                v_u,axial = dP_jacket / (p_int * L_transfer)
            and the jacket's flexural compression must be built up from the
            point of contraflexure, ACI 318-19 16.4.5.1 concept:
                v_u,flex  = C_jacket / (p_int * L_shear_span)
            These act on the same plane in the same direction and are added.

        A discontinuous jacket typically demands 3-5x the interface capacity
        of a continuous one.  If the calculator says your dowels are
        impossible, the answer is almost always to make the jacket
        continuous, not to add more dowels.

        Capacity -- ACI 318-19 16.4.4.1 SELECTS between two routes by the
        magnitude of the demand; they are alternatives, not a pair to be
        minimised:

            v_u <= 3.5 MPa  ->  Table 16.4.4.2 (SI), roughened surface with
                                A_v >= A_v,min:
                                v_nh = lam*(1.8 + 0.6*rho_v*f_yt) <= 3.5 MPa
                                v_nh = 0.55 MPa if not roughened
            v_u >  3.5 MPa  ->  shear friction, 22.9.4.2:
                                V_n = mu*A_vf*f_y   ->  v_n = mu*rho_v*f_y
                                mu = 1.4*lam monolithic
                                     1.0*lam roughened ~6 mm (1/4 in.)
                                     0.6     not intentionally roughened
                                                       Table 22.9.4.2
                                v_n <= min(0.2f'c, 3.3+0.08f'c, 11) MPa
                                                       Table 22.9.4.4 (SI)

        BOTH routes are returned.  The 1.8 MPa term in Table 16.4.4.2 is a
        BOND/cohesion contribution -- it is what makes column jacketing
        practical, because dowels alone need uneconomic densities.  But it is
        mobilised only by a clean, roughened, well-cured interface, and it is
        lost if the interface cracks.  Relying on it is a decision the
        Engineer of Record makes explicitly, and it is why ACI 562-16 7.4.3
        requires quantitative pull-off testing (ASTM C1583, minimum 3 tests).

        ACI 562-16 Table 7.4.1.2 thresholds:
            v_u > 0.21 MPa (30 psi) -> quantitative bond testing required
            v_u > 0.41 MPa (60 psi) -> interface reinforcement REQUIRED
        A structural column jacket is always past both.

        ACI 562-16 7.4.1 also requires restrained volume change (jacket
        shrinkage against the restraining core) to be added to v_u.  That
        term is NOT computed here -- it depends on the mix, curing and age
        and must be assessed separately.
        """
        self._require_four_sided("interface_check", "phase 6 (spec s4.7)")
        p_int = 2.0 * (self.b_e + self.h_e)          # interface perimeter, mm
        fc = min(self.conc_e.fc, self.conc_j.fc)
        lam = min(self.conc_e.lam, self.conc_j.lam)
        b_v = self.b_e                                # peeled contact width

        # Two idealisations of the flexural interface demand -- an ELASTIC
        # compatibility check (VQ/Ib) and a PLASTIC force-transfer check
        # (build the jacket's compression force up over the shear span).
        # They are not the same model and neither dominates in general, so
        # both are computed and the larger governs.
        n = self.conc_e.Ec / self.conc_j.Ec
        I_tr = self.Ig_tot - self.Ig_e + n * self.Ig_e
        Q_j = self.B * self.t * (self.H / 2.0 - self.t / 2.0)
        v_elastic = Vu * Q_j / (I_tr * b_v) if I_tr > 0 else 0.0
        v_plastic = C_jacket / (p_int * L_shear_span) if L_shear_span > 0 else 0.0
        v_flex = max(v_elastic, v_plastic)

        if continuity == "continuous":
            v_axial = 0.0
            note = ("jacket continuous through the slab and bearing on the "
                    "foundation; axial load enters by direct bearing, so the "
                    "interface carries the flexural demand only")
        else:
            v_axial = dP_jacket / (p_int * L_transfer) if L_transfer > 0 else 0.0
            note = ("jacket discontinuous; EVERY newton of axial load credited "
                    "to the jacket must also cross the interface")
        v_u = v_flex + v_axial

        # ACI 318-19 Table 22.9.4.2: EVERY mu is a multiple of lambda.  The
        # not-roughened entry was bare, which overstated it by 1/lam -- +33% at
        # lam = 0.75, on the one surface with the least capacity to spare.
        mu = {"monolithic": 1.4 * lam,
              "roughened_6mm": 1.0 * lam,
              "not_roughened": 0.6 * lam,
              "steel": 0.7 * lam}[surface]
        roughened = surface in ("monolithic", "roughened_6mm")

        # route A -- shear friction, 22.9 (dowels only, no cohesion)
        #
        # ACI 318-19 Table 22.9.4.4 splits the cap by surface condition:
        # monolithic, or hardened concrete intentionally roughened to a full
        # 1/4 in. amplitude, take the least of 0.2f'c, (3.3 + 0.08f'c) and
        # 11.0 MPa.  ALL OTHER CASES take the lesser of 0.2f'c and 5.5 MPa.
        # The single unbranched set credited the roughened ceiling to a
        # not-roughened interface -- +40% at f'c = 55 MPa.
        if roughened:
            sf_cap = min(0.2 * fc, 3.3 + 0.08 * fc, 11.0)
        else:
            sf_cap = min(0.2 * fc, 5.5)
        v_n_sf = min(mu * rho_v * fy_dowel, sf_cap)
        rho_req_sf = v_u / (phi * mu * fy_dowel)

        # route B -- composite flexural member, Table 16.4.4.2 (with bond)
        v_nh = min(lam * (1.8 + 0.6 * rho_v * fy_dowel), 3.5) if roughened \
            else 0.55
        rho_req_bond = max((v_u / phi / lam - 1.8) / (0.6 * fy_dowel), 0.0) \
            if roughened else float("inf")

        return dict(
            continuity=continuity, note=note, p_int=p_int, b_v=b_v,
            v_u=v_u, v_flex=v_flex, v_axial=v_axial, mu=mu, rho_v=rho_v,
            v_elastic=v_elastic, v_plastic=v_plastic,
            # route A
            v_n_dowel=v_n_sf, phi_v_n_dowel=phi * v_n_sf,
            DCR_dowel=v_u / (phi * v_n_sf) if v_n_sf > 0 else float("inf"),
            rho_v_req_dowel=rho_req_sf, sf_cap=sf_cap,
            # route B
            v_n_bond=v_nh, phi_v_n_bond=phi * v_nh,
            DCR_bond=v_u / (phi * v_nh) if v_nh > 0 else float("inf"),
            rho_v_req_bond=rho_req_bond,
            ok_dowel_only=(phi * v_n_sf >= v_u),
            ok_with_bond=(phi * v_nh >= v_u),
            aci562_bond_test_required=(v_u > 0.21),
            aci562_reinf_required=(v_u > 0.41))

    def dowel_layout(self, rho_v_req: float, dia: float, n_per_row: int):
        """Convert a required interface reinforcement ratio into a dowel
        spacing.  rho_v = A_vf/(p_int * s)  ->  s = A_vf/(rho_v * p_int).

        Dowels are post-installed adhesive anchors; their tensile capacity
        must also be checked to ACI 318-19 17.6 (steel 17.6.1, breakout
        17.6.2, bond 17.6.5 with tau_cr from ACI 355.4 qualification,
        x0.8 for seismic per Table 17.6.5.2.5 note 2).  This helper only
        gives the geometry.
        """
        self._require_four_sided("dowel_layout", "phase 6 (spec s4.8)")
        A = math.pi * dia ** 2 / 4.0 * n_per_row
        p_int = 2.0 * (self.b_e + self.h_e)
        s = A / (rho_v_req * p_int)
        return dict(A_row=A, s_required=s, dia=dia, n_per_row=n_per_row)

    # ---------------- detailing checks -----------------------------------
    def detailing_checks(self, db_long_min: float, L_clear: float):
        """ACI 318-19 detailing limits relevant to the jacketed column."""
        self._require_four_sided("detailing_checks", "phase 9 (spec s4.9)")
        rho_tot = self.Ast_tot / self.Ag_tot
        s_max = min(16 * db_long_min, 48 * self.ties_j.db,
                    min(self.B, self.H))
        db_tie_min = 10.0 if db_long_min <= 32.0 else 13.0
        return [
            dict(item="Total longitudinal ratio",
                 clause="ACI 318-19 10.6.1.1",
                 value=f"rho = {rho_tot*100:.2f}%",
                 limit="0.01 Ag <= Ast <= 0.08 Ag",
                 ok=0.01 <= rho_tot <= 0.08),
            dict(item="Jacket tie spacing",
                 clause="ACI 318-19 25.7.2.1(b)",
                 value=f"s = {self.ties_j.s:.0f} mm",
                 limit=f"<= least of 16db_l, 48db_t, least dim = {s_max:.0f} mm",
                 ok=self.ties_j.s <= s_max),
            dict(item="Jacket tie diameter",
                 clause="ACI 318-19 25.7.2.2",
                 value=f"{self.ties_j.db:.0f} mm",
                 limit=f">= {db_tie_min:.0f} mm",
                 ok=self.ties_j.db >= db_tie_min),
            dict(item="Jacket thickness",
                 clause="No code minimum for COLUMN jackets; ASCE 41-23 "
                        "C7.7 gives 3 in. (75 mm) for WALL boundary jackets",
                 value=f"t = {self.t:.0f} mm",
                 limit=">= 75 mm (interpretation); >= 100 mm practical for "
                       "cast-in-place with bars and ties",
                 ok=self.t >= 75.0),
            dict(item="Slenderness of jacketed column",
                 clause="ACI 318-19 6.2.5",
                 value=f"k*lu/r ~ {L_clear/(0.3*min(self.B,self.H)):.1f} (k=1)",
                 limit="<= 22 to neglect slenderness (braced, no moments)",
                 ok=L_clear / (0.3 * min(self.B, self.H)) <= 22.0),
        ]

    # ---------------- monolithic coefficients ----------------------------
    #
    # Source tables, transcribed exactly from Lampropoulos, Tsioulou & Dritsos
    # (2012), J. Earthquake Engineering 16(7):1023-1042, Tables 2-5
    # (PDF p17-18, printed 1038-1039).  Tuples are (K_F, K_k, K_dy, K_du).
    #
    # four-sided:   Table 2  r  = 1.5, versus nu      Table 3  nu = 0.2, versus r
    # three-sided:  Table 4  r  = 1.0, versus nu      Table 5  nu = 0.2, versus r
    _MONO_T2 = {0.2: (0.85, 0.55, 1.50, 1.00), 0.4: (0.70, 0.50, 1.30, 1.10)}
    _MONO_T3 = {0.5: (0.90, 0.70, 1.30, 0.95), 1.5: (0.85, 0.55, 1.50, 1.00),
                4.0: (0.75, 0.75, 1.05, 2.85)}
    _MONO_T4 = {0.2: (0.80, 0.55, 1.40, 1.10), 0.4: (0.70, 0.35, 1.90, 1.70)}
    _MONO_T5 = {0.5: (0.85, 0.55, 1.70, 1.15), 1.5: (0.80, 0.55, 1.40, 1.10),
                4.0: (0.75, 0.25, 3.00, 1.00)}
    # Reported by the paper for comparison only -- NOT design values.
    _MONO_EC8_P3 = (0.90, 0.95, 1.20, 1.00)
    _MONO_GRECO = (0.90, 0.80, 1.25, 0.80)

    _MONO_NU_RANGE = (0.1, 0.4)
    _MONO_R_RANGE = (0.5, 4.0)

    @property
    def n_faces_jacketed(self) -> int:
        """Number of enlarged faces.  Alias of n_faces, kept for phase 8.

        This was a stub reading the single thickness while the per-face fields
        of CODE_IMPLEMENTATION_PARTIAL_JACKET SS2.1 were unimplemented.  Phase 1
        landed those fields; leaving the stub in place made
        monolithic_coefficients see 4 faces for EVERY jacket, which killed its
        refusal branch and handed a one-sided jacket the four-sided
        Lampropoulos tables -- worse than the borrowing TN-RET-002 TODO C10
        forbids, and reported as "n_faces": 4, misstating the geometry.
        """
        return self.n_faces

    @staticmethod
    def _mono_interp(table: dict, x: float) -> tuple:
        """Piecewise-linear interpolation across a coefficient table.

        Clamped at the end anchors.  The paper states that intermediate values
        follow a linear distribution but publishes no values outside the
        tabulated span, so extrapolating past an anchor would be inventing
        data; the caller reports the clamp through the out_of_range flag.
        """
        xs = sorted(table)
        if x <= xs[0]:
            return table[xs[0]]
        if x >= xs[-1]:
            return table[xs[-1]]
        for lo, hi in zip(xs, xs[1:]):
            if lo <= x <= hi:
                f = (x - lo) / (hi - lo)
                return tuple(a + f * (b - a)
                             for a, b in zip(table[lo], table[hi]))
        raise AssertionError("unreachable")  # pragma: no cover

    def monolithic_coefficients(self, N: float) -> dict:
        """Monolithic coefficients K_F, K_k, K_dy, K_du after Lampropoulos,
        Tsioulou & Dritsos (2012), J. Earthquake Engineering 16(7):1023-1042,
        Tables 2-5 (PDF p17-18, printed 1038-1039).

            K_F  = F_max,strengthened / F_max,monolithic
            K_k  = (Fy/dy)_str / (Fy/dy)_mon
            K_dy, K_du = deflection or rotation at yield and at failure

        Enter with the normalized axial load, their Eq. (2) (PDF p9):

            nu = N / (A_co * f_co + A_cj * f_cj)      CONCRETE areas only

        and the concrete area ratio r = A_cj / A_co.  N is in newtons,
        compression positive, consistent with the rest of this module.

        RETURNS None WITH A REASON for anything other than a three- or
        four-sided jacket.  The paper studied nothing else, and the three-sided
        values are an UPPER BOUND for a smaller jacket, not an estimate.
        See TN-RET-002 TODO C10 -- BLOCKING for a one- or two-sided jacket.

        Valid range nu 0.1-0.4, r 0.5-4.0.  Outside it the clamped value is
        returned together with out_of_range; this method never extrapolates.

        The coefficients already include jacket shrinkage (800 microstrain
        free, 400 with creep) and an interface at friction 1.0 with COHESION
        0.0 MPa.  Do not add a separate shrinkage penalty on top without
        resolving TN-RET-002 TODO C12.

        THE CALLER MUST NOT APPLY THESE SILENTLY.  Whether K_F multiplies
        phi*M_n or M_n, and whether phi and K_F may both be applied, is
        unsettled -- TN-RET-001 TODO C8 and C9, both BLOCKING.  Report the
        monolithic and the adjusted result side by side, the same way
        interface_check() returns both capacity routes and picks neither.

        Two-way interpolation.  The paper gives one-dimensional tables and
        says intermediate values follow a linear distribution without saying
        how to combine them.  Both orders are evaluated:

            order_a: K = K_r(r)   * K_nu(nu) / K_nu(nu_anchor)
            order_b: K = K_nu(nu) * K_r(r)   / K_r(r_anchor)

        For a four-sided jacket the two tables agree at their crossing point
        and the orders coincide.  For a three-sided jacket they do NOT: Table 5
        interpolated to r = 1.0 gives K_F 0.825 and K_dy 1.550 against Table 4's
        own 0.80 and 1.40.  That is a source inconsistency, not an arithmetic
        one, so both raw branches are returned and the conservative one is
        selected -- min for the strength coefficients, max for the deflection
        coefficients.  TN-RET-002 TODO C11.
        """
        keys = ("K_F", "K_k", "K_dy", "K_du")
        n_faces = self.n_faces_jacketed

        if n_faces not in (3, 4):
            return {
                "applicable": False,
                "n_faces": n_faces,
                "reason": (
                    f"No published coefficients cover a {n_faces}-sided jacket. "
                    "Lampropoulos et al. (2012) studied four- and three-sided "
                    "jackets only; the three-sided values are an UPPER BOUND "
                    "for a smaller jacket, not an estimate, and must not be "
                    "borrowed. TN-RET-002 TODO C10 -- BLOCKING."
                ),
                "coefficients": None,
            }

        # Eq. (2): CONCRETE areas only -- gross concrete, steel not deducted.
        A_co, A_cj = self.Ag_e, self.Ag_j
        denom = A_co * self.conc_e.fc + A_cj * self.conc_j.fc
        nu = N / denom if denom > 0 else float("nan")
        r = A_cj / A_co if A_co > 0 else float("nan")

        if n_faces == 4:
            T_nu, T_r, nu_anchor, r_anchor = self._MONO_T2, self._MONO_T3, 0.2, 1.5
        else:
            T_nu, T_r, nu_anchor, r_anchor = self._MONO_T4, self._MONO_T5, 0.2, 1.0

        K_nu = self._mono_interp(T_nu, nu)
        K_nu_anchor = self._mono_interp(T_nu, nu_anchor)
        K_r = self._mono_interp(T_r, r)
        K_r_anchor = self._mono_interp(T_r, r_anchor)

        order_a = tuple(kr * (kn / kna) for kr, kn, kna
                        in zip(K_r, K_nu, K_nu_anchor))
        order_b = tuple(kn * (kr / kra) for kn, kr, kra
                        in zip(K_nu, K_r, K_r_anchor))

        # Conservative branch: a lower strength coefficient and a larger
        # deflection coefficient are both the unfavourable reading.
        chosen = (min(order_a[0], order_b[0]), min(order_a[1], order_b[1]),
                  max(order_a[2], order_b[2]), max(order_a[3], order_b[3]))

        nu_lo, nu_hi = self._MONO_NU_RANGE
        r_lo, r_hi = self._MONO_R_RANGE
        out_of_range = not (nu_lo <= nu <= nu_hi and r_lo <= r <= r_hi)
        # The tables' own anchors are narrower than the declared valid range.
        clamped_nu = not (min(T_nu) <= nu <= max(T_nu))
        clamped_r = not (min(T_r) <= r <= max(T_r))

        return {
            "applicable": True,
            "n_faces": n_faces,
            "nu": nu,
            "r": r,
            "coefficients": dict(zip(keys, chosen)),
            "order_a": dict(zip(keys, order_a)),
            "order_b": dict(zip(keys, order_b)),
            "orders_agree": all(abs(a - b) <= 5e-3
                                for a, b in zip(order_a, order_b)),
            "out_of_range": out_of_range,
            "clamped_nu": clamped_nu,
            "clamped_r": clamped_r,
            "reference_only": {"EC8_Part3": dict(zip(keys, self._MONO_EC8_P3)),
                               "Greco": dict(zip(keys, self._MONO_GRECO))},
            "reason": (
                "Coefficients returned for reporting only. They are defined on "
                "a member's lateral load-deflection response, not on a "
                "section's phi*M_n. Whether K_F multiplies phi*M_n or M_n, and "
                "whether phi and K_F may both be applied, is unsettled -- "
                "TN-RET-001 TODO C8/C9, BLOCKING. Report both results."
            ),
            "clause": ("Lampropoulos, Tsioulou & Dritsos (2012) Tables 2-5; "
                       "TN-RET-001 Core Technical Requirements 4; TN-RET-002 s10"),
        }

    # ---------------- summary --------------------------------------------
    def stiffness_summary(self):
        """Gross flexural stiffness before and after.  A jacket is a
        stiffness intervention as much as a strength one: EI rises roughly
        with H^3*B, the column attracts more seismic force, and the demand
        model MUST be re-run.  Jacketing some columns in a storey and not
        others redistributes demand and can create an irregularity."""
        self._require_four_sided("stiffness_summary", "phase 4 (spec s4.6)")
        EI_e = self.conc_e.Ec * self.Ig_e
        # transformed to jacket concrete
        n = self.conc_e.Ec / self.conc_j.Ec
        I_tr = self.Ig_tot - self.Ig_e + n * self.Ig_e
        EI_j = self.conc_j.Ec * I_tr
        return dict(Ig_existing=self.Ig_e, Ig_jacketed=self.Ig_tot,
                    I_transformed=I_tr, EI_existing=EI_e, EI_jacketed=EI_j,
                    EI_ratio=EI_j / EI_e, Ag_ratio=self.Ag_tot / self.Ag_e,
                    period_factor=math.sqrt(EI_e / EI_j))


# --------------------------------------------------------------------------
# 3.  REPORTING HELPERS
# --------------------------------------------------------------------------

kN = 1e-3
kNm = 1e-6


def envelope(P, M):
    """Return the P-M points sorted for plotting a closed envelope."""
    idx = np.argsort(-np.asarray(P))
    return np.asarray(P)[idx], np.asarray(M)[idx]


def capacity_at_P(P_arr, M_arr, P_target):
    """Interpolate the moment capacity at a given axial load.

    Returns NaN if P_target is outside the axial range of the surface --
    i.e. the section cannot carry that axial load at any eccentricity.
    Silently clamping would report a capacity the section does not have.
    """
    P_arr = np.asarray(P_arr, float)
    M_arr = np.asarray(M_arr, float)
    if P_target > P_arr.max() or P_target < P_arr.min():
        return float("nan")
    order = np.argsort(P_arr)
    return float(np.interp(P_target, P_arr[order], M_arr[order]))


def utilisation(P_arr, M_arr, Pu, Mu):
    """Demand/capacity ratio against a P-M surface, measured radially from
    the origin through the demand point (Pu, Mu).  Returns >1 if outside."""
    Mcap = capacity_at_P(P_arr, M_arr, Pu)
    if math.isnan(Mcap):
        return float("inf")
    return abs(Mu) / Mcap if Mcap > 0 else float("inf")
