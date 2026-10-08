# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""BELFEM input-file topology section, generated from the actual physical
tags of the final mesh.

Curves are defined for BELFEM as sideset intersections (plane @ tape):
odd ids on the back plane, even ids on the front plane, per tape. The
periodic block names the vertex triples (type-15 elements 1..10): three
front quarter points (2, 3, 4) map to the matching back quarter points
(7, 8, 9); BELFEM derives the affine cap-to-cap map from them, here the
translation by the cable length. Node 1 (front cap center) is not part
of the triple: three points already fix that map.

The boundary conditions section drives current through the tapes, ONE
condition per galvanically connected conductor, defined on ONE crossing
of that conductor: the back-plane curve of its first tape as the input
curve and the front-plane curve of the same tape as the output curve.
Tapes that continue into each other across the periodic planes are one
conductor, so the current is the same on every crossing and one crossing
fixes it (the corc_periodic_bc deck follows the same rule). Listing all
k crossings of one conductor in a single bracket group is NOT the same
condition: the group chain then links the conductor k times, and the
terminal transform of the cohomology module imposes the amplitude k times
per crossing.
"""


def topology_section(volume_count, tape_physicals, boundary_physical,
                     front_physical, back_physical, conductor_groups=None,
                     tape_current=30.0):
    if conductor_groups is None:
        conductor_groups = [[k] for k in range(len(tape_physicals))]
    lines = []
    a = lines.append
    a("topology")
    a("{")
    a("\tthinshell : tape")
    a("\t{")
    if len(tape_physicals) > 1:
        a("\t\tsidesets : {:d}:{:d} ;".format(tape_physicals[0], tape_physicals[-1]))
    else:
        a("\t\tsidesets : {:d} ;".format(tape_physicals[0]))
    a("\t}")
    a("")
    a("\t// the core section (block 1) is air: the copper core is not modeled")
    a("\tair")
    a("\t{")
    if volume_count > 1:
        a("\t\tblocks : 1:{:d} ;".format(volume_count))
    else:
        a("\t\tblocks : 1 ;")
    a("\t}")
    a("")
    a("\t// curves: intersection of a cap plane with a tape surface;")
    a("\t// odd = back plane ({:d}), even = front plane ({:d})".format(
        back_physical, front_physical))
    a("\tcurves")
    a("\t{")
    input_curves = []
    output_curves = []
    cid = 0
    for t in tape_physicals:
        cid += 1
        a("\t\t{:d} : {:d} @ {:d} ;".format(cid, back_physical, t))
        input_curves.append(cid)
        cid += 1
        a("\t\t{:d} : {:d} @ {:d} ;".format(cid, front_physical, t))
        output_curves.append(cid)
    a("\t}")
    a("")
    a("\t// periodicity: three front quarter points (type-15 elements 2, 3, 4)")
    a("\t// map to the matching back quarter points (7, 8, 9) by the")
    a("\t// translation along z; the front cap center (node 1) is not part")
    a("\t// of the triple")
    a("\tperiodic")
    a("\t{")
    a("\t\tsource : 2, 3, 4 ;")
    a("\t\ttarget : 7, 8, 9 ;")
    a("\t}")
    a("}")
    a("")
    a("boundary conditions")
    a("{")
    crossings = len(conductor_groups[0])
    a("\t// one condition per conductor, defined on ONE crossing of it: tapes")
    a("\t// that continue into each other across the periodic planes are one")
    a("\t// conductor crossing the plane {:d} time(s), the current is the same".format(crossings))
    a("\t// on every crossing, and the amplitude is the tape current")
    a("\tcurrent")
    a("\t{")

    def bracket(curves):
        return ",".join("[{}]".format(curves[g[0]]) for g in conductor_groups)

    a("\t\tinput curves  : {} ;".format(bracket(input_curves)))
    a("\t\toutput curves : {} ;".format(bracket(output_curves)))
    a("\t\ttype          : ramp ;")
    a("\t\tamplitude     : {:g} A ;".format(tape_current))
    a("\t\tperiod        : 10 s ;")
    a("\t\toffset        : -0.01 s ;")
    a("\t}")
    a("}")
    return "\n".join(lines) + "\n"


def write_topology(path, **kw):
    with open(path, "w") as f:
        f.write(topology_section(**kw))
