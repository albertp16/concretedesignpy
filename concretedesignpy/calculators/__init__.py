# MIT License
#
# Copyright (c) Albert Pamonag Engineering Consultancy
#
# Permission is hereby granted, free of charge, to any person obtaining a copy
# of this software and associated documentation files (the "Software"), to deal
# in the Software without restriction, including without limitation the rights
# to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
# copies of the Software, and to permit persons to whom the Software is
# furnished to do so, subject to the following conditions:
#
# The above copyright notice and this permission notice shall be included in
# all copies or substantial portions of the Software.

"""
concretedesignpy.calculators
============================

Backend calculation modules for reinforced concrete design.
Extracted and consolidated from development notebooks.

Modules:
    - beam_moment: Beam flexural capacity (neutral axis solver)
    - beam_shear: Concrete and steel shear strength
    - beam_torsion: Torsion design checks per ACI/NSCP
    - beam_deflection: Short-term and long-term deflection
    - column_interaction: P-M interaction diagram generation
    - column_flexural: Minimum flexural strength ratio check
    - joint_shear: Joint shear verification for special moment frames
    - mander: Mander's confined concrete model
    - development_length: Hook geometry per NSCP 2015 Section 425
    - moment_curvature: 6-point moment-curvature relationship
    - rebar_layout: Rebar coordinate generation for rectangular sections
    - alternative_inertia: Alternative moment of inertia per NSCP 2015

Standards:
    NSCP 2015, ACI 318-19, ASCE 41
"""

from concretedesignpy.calculators.beam_moment import calculate_beam_moment
from concretedesignpy.calculators.beam_shear import (
    compute_concrete_shear_strength,
    compute_steel_shear_strength,
    compute_shear_spacing,
)
from concretedesignpy.calculators.beam_torsion import torsion_design
from concretedesignpy.calculators.beam_deflection import deflection_computation
from concretedesignpy.calculators.column_interaction import (
    generate_interaction_diagram,
)
from concretedesignpy.calculators.column_flexural import (
    check_min_flexural_strength,
)
from concretedesignpy.calculators.joint_shear import joint_shear_check
from concretedesignpy.calculators.mander import (
    confined_strength_ratio,
    confined_stress_strain,
)
from concretedesignpy.calculators.development_length import (
    bend_and_hook_deformed_bars,
    bend_and_hook_stirrups,
)
from concretedesignpy.calculators.moment_curvature import moment_curvature_analysis
from concretedesignpy.calculators.rebar_layout import section_generate_rect
from concretedesignpy.calculators.alternative_inertia import (
    alternative_inertia_column_wall,
    alternative_inertia_beam_slab,
)
