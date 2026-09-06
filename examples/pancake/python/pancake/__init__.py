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
