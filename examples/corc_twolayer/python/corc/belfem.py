"""BELFEM input-file topology section, generated from the actual physical
tags of the final mesh (see tmp/input.conf for the hand-written sample).

Curves are defined for BELFEM as sideset intersections (plane @ tape):
odd ids on the back plane, even ids on the front plane, per tape — the
same order as the sample. The periodic block names the vertex triples
(type-15 elements 1..10): three front quarter points (2, 3, 4) map to the
matching back quarter points (7, 8, 9); BELFEM derives the affine cap-to-cap
map from them. Node 1 (front cap center) is not part of the triple: three
points already fix that map, and the cap center lies on the axis.

The boundary conditions section drives current through the tapes: the odd
(back-plane) curves are the input curves, the even (front-plane) curves are
the output curves.
"""


def topology_section(volume_count, tape_physicals, boundary_physical,
                     front_physical, back_physical, solder_corrections=None,
                     solder_front_physicals=None, solder_back_physicals=None):
    lines = []
    a = lines.append
    # solder mode: the interlayer volumes 2..volume_count-1 conduct, the tapes
    # and the solder are ONE conductor, and the current condition is a single
    # bulk terminal pair on the solder cap faces ( see the emit comment in
    # postprocess.py and tapestack3d.geo "The current terminal" )
    solder_mode = bool(solder_corrections) and solder_front_physicals and solder_back_physicals
    if solder_corrections:
        a("// If the interlayer volumes are meshed as solder instead of air,")
        a("// the solder material needs this 'density correction' (the real-")
        a("// solder fraction of the annulus; each layer's tape stack lies")
        a("// radially outward of its shell, so annulus l carries layer l):")
        for l, f in enumerate(solder_corrections):
            a("//     Volume_{:d} : {:.4f}".format(l + 2, f))
        a("")
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
    a("\tair")
    a("\t{")
    if solder_mode:
        a("\t\tblocks : 1, {:d} ;".format(volume_count))
    elif volume_count > 1:
        a("\t\tblocks : 1:{:d} ;".format(volume_count))
    else:
        a("\t\tblocks : 1 ;")
    a("\t}")
    if solder_mode:
        a("")
        a("\tconductor : solder")
        a("\t{")
        if volume_count > 3:
            a("\t\tblocks : 2:{:d} ;".format(volume_count - 1))
        else:
            a("\t\tblocks : 2 ;")
        a("\t\tmaterial : solder ;")
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
    a("\t// map to the matching back quarter points (7, 8, 9); three points")
    a("\t// already fix the affine cap-to-cap map, so the front cap center")
    a("\t// (node 1) is not part of the triple")
    a("\tperiodic")
    a("\t{")
    a("\t\tsource : 2, 3, 4 ;")
    a("\t\ttarget : 7, 8, 9 ;")
    a("\t}")
    a("}")
    a("")
    a("boundary conditions")
    a("{")
    a("\tcurrent")
    a("\t{")
    if solder_mode:
        # one condition, TOTAL current: the terminals are the solder end
        # faces; input = back plane ( odd curves ), output = front plane.
        # The curves stay in the topology block ( thin-shell terminal
        # curves ) but must NOT appear here: a terminal key takes
        # precedence and a curve form would sum the tape loops only
        a("\t\tinput terminals : {} ;".format(",".join(str(c) for c in solder_back_physicals)))
        a("\t\toutput terminals : {} ;".format(",".join(str(c) for c in solder_front_physicals)))
    else:
        # insulated tapes: one bracket per tape, per-tape current
        a("\t\tinput curves : {} ;".format(",".join("[{:d}]".format(c) for c in input_curves)))
        a("\t\toutput curves : {} ;".format(",".join("[{:d}]".format(c) for c in output_curves)))
    a("\t\ttype : ramp ;")
    a("\t\tamplitude : 60 A ;")
    a("\t\tperiod : 10 s ;")
    a("\t\toffset : -0.01 s ;")
    a("\t}")
    a("}")
    return "\n".join(lines) + "\n"


def write_topology(path, **kw):
    with open(path, "w") as f:
        f.write(topology_section(**kw))
