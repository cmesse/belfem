# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""Geometry of pancake HTS coils: Lame-spiral windings with intrinsically defined leads."""

from .basecurve import Basecurve, RigidTransform, rotation_x, rotation_z
from .spiral import PlanarSpiral, LameSpiral
from .lead import Segment, Straight, Twist, Bend, Release, Lead, smoothstep
from .pancake import PancakeCurve
from .export import frame_table, tape_corners, write_vtk_polyline, write_stl_tape
from .gmsh_geometry import CrossSection, Geometry, Tape, TapeBlock
from .tapestack import TapeStack
from .fem import PancakeMesh, BoxDomain, WedgeDomain, leads_to_plane, leads_to_cylinder, FemError

__all__ = [
    "Basecurve", "RigidTransform", "rotation_x", "rotation_z",
    "PlanarSpiral", "LameSpiral",
    "Segment", "Straight", "Twist", "Bend", "Release", "Lead", "smoothstep",
    "PancakeCurve",
    "frame_table", "tape_corners", "write_vtk_polyline", "write_stl_tape",
    "CrossSection", "Geometry", "Tape", "TapeBlock",
    "TapeStack", "PancakeMesh", "BoxDomain", "WedgeDomain", "leads_to_plane", "leads_to_cylinder", "FemError",
]
