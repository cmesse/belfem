# BELFEM -- The Berkeley Lab Finite Element Framework
# Copyright (c) 2026, The Regents of the University of California,
# through Lawrence Berkeley National Laboratory (subject to receipt of any required
# approvals from the U.S. Dept. of Energy).  All rights reserved.
#
# Developers: Christian Messe, Gregory Giard
#
# See the top-level LICENSE file for the complete license and disclaimer.

"""BELFEM topology and boundary-condition sections for the pancake.

The tape is one thin shell made of half-turn sidesets. The solder
between the turns, the sliver inside the first turn and the two lead
bars are conductor blocks. The terminals are the bars' end faces on
the domain boundary; the tape's terminal curves are their intersections
with the first and the last tape piece, and BELFEM attaches those to
the thin shell by itself, as in the racetrack example."""


def topology_section(ids, n_solder, n_tape, inner_piece, outer_piece, tape_current,
                     solder_label="solder"):
    first_tape = ids["Tape_1"]
    last_tape = ids["Tape_{:d}".format(n_tape)]
    tape_in = ids["Tape_{:d}".format(inner_piece + 1)]
    tape_out = ids["Tape_{:d}".format(outer_piece + 1)]
    lines = []
    a = lines.append
    a("topology")
    a("{")
    a("\tthinshell : tape")
    a("\t{")
    a("\t\tsidesets : {:d}:{:d} ;".format(first_tape, last_tape))
    a("\t}")
    a("")
    a("\t// half-turn pieces of the solder between the turns, the sliver")
    a("\t// inside the first turn (blocks 1:2), and the two lead bars")
    a("\t// (blocks {:d}:{:d})".format(ids["LeadIn"], ids["LeadOut"]))
    a("\tconductor : {}".format(solder_label))
    a("\t{")
    a("\t\tblocks   : 1:{:d} ;".format(ids["LeadOut"]))
    a("\t\tmaterial : {} ;".format(solder_label))
    a("\t}")
    a("")
    a("\tair")
    a("\t{")
    a("\t\tblock : {:d} ;".format(ids["Air"]))
    a("\t}")
    a("")
    a("\t// sideset {:d} (SolderAir), the solder/air contact, is not listed:".format(
        ids["SolderAir"]))
    a("\t// BELFEM types it as a conductor-air interface from its two blocks")
    a("")
    a("\t// terminal curves: the tape's edges on the bar end faces")
    a("\tcurves")
    a("\t{")
    a("\t\t1 : {:d} @ {:d} ;   // inner lead piece".format(ids["TerminalIn"], tape_in))
    a("\t\t2 : {:d} @ {:d} ;   // last piece of the outer lead".format(ids["TerminalOut"], tape_out))
    a("\t}")
    a("}")
    a("")
    a("boundary conditions")
    a("{")
    a("\tcurrent")
    a("\t{")
    a("\t\tinput terminals  : {:d} ;".format(ids["TerminalIn"]))
    a("\t\toutput terminals : {:d} ;".format(ids["TerminalOut"]))
    a("\t\ttype             : ramp ;")
    a("\t\tamplitude        : {:g} A ;".format(tape_current))
    a("\t\tperiod           : 10 s ;")
    a("\t\toffset           : -0.01 s ;")
    a("\t}")
    a("}")
    return "\n".join(lines) + "\n"


def write_topology(path, **kw):
    with open(path, "w") as f:
        f.write(topology_section(**kw))
