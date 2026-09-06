#!/usr/bin/env python3
"""
Build the BELFEM thermo.inp and trans.inp from the NASA CEA tables.

The output carries only the species BELFEM needs, listed in species.txt, and each
record is complete: where a low temperature fit can be produced, it is merged into
the record as an additional interval rather than kept in a second file. There is
therefore one table per property, no overlay, and no load order contract.

Why extract rather than ship the vendor file
--------------------------------------------
The CEA thermo.inp in circulation is not purely a US government work. Later
versions carry propellant components contributed by third parties, and one block
taken from a Russian database, none of which come with a licence. Those additions
are all solid and propellant materials, so a table restricted to the species below
contains none of them. Extracting is therefore cleaner than filtering, and it also
removes the need to redistribute a large file most of which is unused.

Provenance is recorded per record: the vendor intervals keep their original
reference codes, and any interval this tool generated is marked in the comment.

Usage
-----
    python build_tables.py --vendor-thermo <thermo.inp> --vendor-trans <trans.inp> \\
                           --species species.txt --out-dir ../../share/fluid

Requires numpy and CoolProp, and imports the fitting machinery from nasa9_lowT.py.
"""

import argparse
import os
import sys

sys.path.insert( 0, os.path.dirname( os.path.abspath( __file__ ) ) )

import nasa9_lowT as lowT

try:
    import CoolProp
    import CoolProp.CoolProp as CP
except ImportError:                                     # pragma: no cover
    sys.exit( "CoolProp is required: pip install CoolProp" )

#-------------------------------------------------------------------------------
# readers that keep the whole record, not just its lowest interval
#-------------------------------------------------------------------------------

def fix_label( raw ):
    """
    Map a CEA label to the one BELFEM uses, matching fn_GT_fix_label.cpp. The
    output keeps the CEA spelling, since the reader applies this itself; the
    mapping is only needed so species.txt can be written in BELFEM labels.
    """
    if raw == "C2H3,vinyl":
        return "C2H3"
    return raw.replace( "AR", "Ar" ).replace( "AL", "Al" ).replace( "CL", "Cl" )

def read_thermo_record( lines, label ):
    """
    Return a species record as
    ( comment_line, composition_line, [ ( header, coeff_line_a, coeff_line_b ), ... ] ).

    The lines are kept verbatim so that vendor data is reproduced byte for byte,
    apart from the two fields this tool has to change.
    """
    for i, line in enumerate( lines ):
        if not line.strip() or line[ 0 ] in " !-":
            continue
        if fix_label( line.split()[ 0 ] ) != label:
            continue

        count = int( lines[ i + 1 ][ 0:2 ] )
        intervals = []
        for k in range( count ):
            base = i + 2 + 3 * k
            intervals.append( ( lines[ base ], lines[ base + 1 ], lines[ base + 2 ] ) )

        return lines[ i ], lines[ i + 1 ], intervals

    return None

def read_trans_record( lines, label ):
    """
    Return ( header_line, { "V": [ row, ... ], "C": [ row, ... ] } ) for a pure
    species, keeping the rows verbatim. Interaction pairs carry a second name and
    are not returned here.
    """
    i = 0
    while i < len( lines ):
        line = lines[ i ]

        if not ( len( line ) > 38 and line[ 34 ] == "V" and line[ 36 ] == "C"
                 and line[ 35 ].isdigit() and line[ 37 ].isdigit() ):
            i += 1
            continue

        n_v, n_c = int( line[ 35 ] ), int( line[ 37 ] )

        if fix_label( line[ 0:15 ].strip() ) == label and not line[ 15:30 ].strip():
            rows = { "V": [], "C": [] }
            for k in range( n_v + n_c ):
                row = lines[ i + 1 + k ]
                rows[ row[ 1:2 ] ].append( row )
            return line, rows

        i += 1 + n_v + n_c

    return None

#-------------------------------------------------------------------------------
# transport by isotope scaling
#-------------------------------------------------------------------------------

#! child label -> parent label, for species whose transport is derived
ISOTOPES = { "D": "H", "OD": "OH", "HD": "H2" }

def molar_mass( composition_line ):
    """molar mass in g/mol, the second to last field of the composition line"""
    return float( composition_line.split()[ -2 ] )

def scale_transport( rows, header, m_child, m_parent ):
    """
    Derive a transport record from its lighter isotopologue.

    In the dilute gas the viscosity follows Chapman-Enskog,

        eta = 5/16 * sqrt( pi m k T ) / ( pi sigma^2 Omega )

    and isotopologues share the electronic potential, so sigma and Omega are the
    same and only the mass prefactor changes. Thermal conductivity carries the
    additional 1/M of eta * R/M, so

        eta_child / eta_parent    = sqrt( m_child / m_parent )
        lambda_child / lambda_parent = sqrt( m_parent / m_child )

    NASA generated the D2 transport from H2 this way: their ratio is 1.4134 to
    1.4137 against sqrt(2) = 1.41421 over 300 to 4000 K.

    Because the CEA form is ln( y/scale ) = A lnT + B/T + C/T^2 + D, a constant
    factor moves only D. The other three coefficients are the parent's unchanged,
    which makes the derivation easy to check by eye.

    The conductivity relation is exact only for the translational part. For a
    molecule the internal modes do not scale with mass, which is why the measured
    D2/H2 conductivity ratio drifts from 0.707 to 0.738. For OD from OH the mass
    changes by 6 percent, so that drift is small.
    """
    import math
    shift = { "V": 0.5 * math.log( m_child / m_parent ),
              "C": 0.5 * math.log( m_parent / m_child ) }

    out = { "V": [], "C": [] }
    for kind in ( "V", "C" ):
        for row in rows[ kind ]:
            c = [ lowT._f( row[ 20+15*j : 35+15*j ] ) for j in range( 4 ) ]
            c[ 3 ] += shift[ kind ]
            out[ kind ].append( row[ :20 ] + "".join( lowT._e( v ) for v in c ) )

    return out

#-------------------------------------------------------------------------------
# merging
#-------------------------------------------------------------------------------

def set_interval_tmin( header, t_min ):
    """rewrite the lower bound of an interval header, leaving every other column"""
    return f"{t_min:11.3f}" + header[ 11: ]

def merge_thermo( record, fluid, junction, samples, margin ):
    """
    Prepend a generated low temperature interval and lift the lowest vendor
    interval to the junction, so the intervals meet exactly and do not overlap.

    Returns ( intervals, note ) or ( None, reason ).
    """
    comment, composition, intervals = record

    t_triple = CP.PropsSI( "Ttriple", fluid )
    t_low = t_triple + margin

    if t_low >= junction:
        return None, f"triple point {t_triple:.2f} K is above the junction", None

    # the vendor interval the fit has to meet
    import numpy as np
    head = intervals[ 0 ][ 0 ]
    a_ref = [ lowT._f( intervals[ 0 ][ 1 ][ k:k+16 ] ) for k in range( 0, 80, 16 ) ]
    a_ref += [ lowT._f( intervals[ 0 ][ 2 ][ k:k+16 ] ) for k in range( 0, 32, 16 ) ]
    a_ref = np.array( a_ref[ :7 ] )
    b1_ref = lowT._f( intervals[ 0 ][ 2 ][ 48:64 ] )
    b2_ref = lowT._f( intervals[ 0 ][ 2 ][ 64:80 ] )

    blocks, residual = lowT.fit_heat( fluid, t_low, junction,
                                      a_ref, b1_ref, b2_ref, samples )

    jump = lowT.check_thermo_blocks( blocks, a_ref, b1_ref, b2_ref )
    if max( jump.values() ) > 1.0e-9:
        return None, "continuity check failed at an interval edge", None

    # the new intervals reuse the tail of the vendor header, which holds H(298)-H(0)
    tail = head[ 65: ].rstrip()
    generated = []

    for lo, hi, a, b1, b2 in blocks:
        generated.append( (
                f"{lo:11.3f}{hi:11.3f}"
                "7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0"
                + tail.rjust( 2 + len( tail ) ),
                "".join( lowT._d( a[ k ] ) for k in range( 5 ) ),
                "".join( lowT._d( a[ k ] ) for k in range( 5, 7 ) )
                + " " * 16 + lowT._d( b1 ) + lowT._d( b2 ) ) )

    lifted = ( set_interval_tmin( intervals[ 0 ][ 0 ], junction ),
               intervals[ 0 ][ 1 ], intervals[ 0 ][ 2 ] )

    count = f"{len( blocks )} intervals, " if len( blocks ) > 1 else ""

    return ( generated + [ lifted ] + list( intervals[ 1: ] ),
             f"{t_low:.2f} to {junction:.0f} K, {count}{100*residual:.2f}%",
             ( t_low, junction, len( blocks ), residual ) )

def merge_trans( record, fluid, junction, samples, dilute_density, tolerance,
                 allow_nist ):
    """
    Prepend generated low temperature intervals per property and lift the lowest
    vendor interval to the junction. Returns ( rows, notes ) or ( None, reason ).
    """
    import numpy as np
    header, rows = record

    t_triple = CP.PropsSI( "Ttriple", fluid )
    t_low = t_triple + 0.5

    if t_low >= junction:
        return None, "triple point above the junction", None

    out, notes, detail = { "V": [], "C": [] }, {}, {}

    for kind in ( "V", "C" ):
        if not rows[ kind ]:
            continue

        # reference polynomial covering the junction
        c_ref, index = None, 0
        for k, row in enumerate( rows[ kind ] ):
            lo, hi = float( row[ 2:11 ] ), float( row[ 11:20 ] )
            c = np.array( [ lowT._f( row[ 20+15*j : 35+15*j ] ) for j in range( 4 ) ] )
            if lo <= junction <= hi:
                c_ref, index = c, k
                break
        if c_ref is None:
            c_ref = np.array( [ lowT._f( rows[ kind ][ 0 ][ 20+15*j : 35+15*j ] )
                                for j in range( 4 ) ] )

        blocks, residual = lowT.fit_transport( fluid, kind, t_low, junction, c_ref,
                                               samples, dilute_density, tolerance,
                                               allow_nist = allow_nist )
        if blocks is None:
            out[ kind ] = list( rows[ kind ] )
            continue

        for lo, hi, c in blocks:
            out[ kind ].append( f" {kind}{lo:9.1f}{hi:9.1f}"
                                + "".join( lowT._e( v ) for v in c ) )

        # append the vendor rows, lifting the one that spans the junction so its
        # lower bound becomes the junction and the generated intervals meet it
        for k, row in enumerate( rows[ kind ] ):
            out[ kind ].append(
                    ( f" {kind}{junction:9.1f}" + row[ 11: ] ) if k == index
                    else row )

        notes[ kind ] = f"{len(blocks)}+{len(rows[kind])}, {100*residual:.2f}%"
        detail[ kind ] = ( len( blocks ), residual )

    return out, notes, detail

#-------------------------------------------------------------------------------

def read_species( path ):
    out = []
    for line in open( path ):
        line = line.split( "#" )[ 0 ].strip()
        if line and line not in out:
            out.append( line )
    return out

HEADER_THERMO = """\
!
! thermo.inp -- BELFEM fluid property tables, thermodynamic data
!
! Extracted from the NASA CEA thermodynamic database, a work of the United
! States government. Only the species BELFEM needs are carried; see species.txt
! for the list and why each is there. Records are reproduced byte for byte from
! the source, except that a low temperature interval has been added where one
! could be generated, and the lowest original interval then begins at the
! junction so the two meet exactly.
!
! This file holds caloric properties only. Viscosity and thermal conductivity
! are in trans.inp.
!
! Species carrying a generated interval are marked @1 in their comment; see the
! note at the end of this header. See scripts/fluidprop for the method.
!
! The format is unchanged, so a record may be added or replaced by hand, and a
! larger table from another source may be substituted wholesale.
!"""

HEADER_TRANS = """\
transport property coefficients
!
! trans.inp -- BELFEM fluid property tables, transport data
!
! Extracted from the NASA CEA transport database, a work of the United States
! government, restricted to the species in species.txt. Many species carry no
! transport data upstream, in particular the ions and most radicals; that is
! expected, and the readers treat a missing entry as such rather than as an
! error.
!
! This file holds viscosity and thermal conductivity only. Caloric properties
! are in thermo.inp.
!
! Records are marked @1 or @2 where they are not purely the original data; see
! the notes at the end of this header.
!"""

#-------------------------------------------------------------------------------
# the @n notes, which live in the header so the record lines stay short
#-------------------------------------------------------------------------------

def with_marks( line, marks, width = 80 ):
    """
    Append the @n markers to a record line, keeping it inside the fixed column
    layout. The separator narrows to a single space where two would not fit, and
    as a last resort the free text citation is trimmed: the marker is what the
    reader of the file needs, and the citation it qualifies is also carried in
    the note the marker points at.
    """
    if not marks:
        return line.rstrip()

    tail = " ".join( marks )

    for gap in ( "  ", " " ):
        if len( line.rstrip() ) + len( gap ) + len( tail ) <= width:
            return line.rstrip() + gap + tail

    return line.rstrip()[ : width - len( tail ) - 1 ].rstrip() + " " + tail

def note_thermo_block( entries ):
    """@1 for thermo.inp: which species carry a generated caloric interval"""
    if not entries:
        return ""

    out = [ "", "! NOTES", "!",
            "!   @1  A low temperature interval was generated for this species",
            "!       and merged into its record; the lowest original interval",
            "!       then begins at the junction. The fit is to the ideal gas",
            "!       heat capacity of a reference equation of state, with cp,",
            "!       its slope, and the absolute enthalpy and entropy held",
            "!       continuous where the two meet. Fitted range and worst",
            "!       relative error in cp:", "!" ]

    for label, t_low, junction, blocks, residual in entries:
        extra = f"   in {blocks} intervals" if blocks > 1 else ""
        out.append( f"!         {label:<10}{t_low:6.1f} to {junction:.0f} K"
                    f"   {100*residual:5.2f} %{extra}" )

    return "\n".join( out ) + "\n!"

def note_trans_block( fitted, derived ):
    """@1 and @2 for trans.inp: generated intervals, and mass scaled records"""
    out = []

    if fitted:
        out += [ "", "! NOTES", "!",
                 "!   @1  Low temperature intervals were generated for this",
                 "!       species and precede the original data, which then",
                 "!       begins at the junction. They are fitted against the",
                 "!       dilute gas limit, which is what this table",
                 "!       represents. Generated interval count and worst",
                 "!       relative error, per property:", "!" ]

        for label, detail in fitted:
            parts = []
            for kind, name in ( ( "V", "viscosity" ), ( "C", "conductivity" ) ):
                if kind in detail:
                    n, residual = detail[ kind ]
                    parts.append( f"{name} {n} {100*residual:5.2f} %" )
            out.append( f"!         {label:<8}" + "   ".join( parts ) )

    if derived:
        out += ( [ "!" ] if fitted else [ "", "! NOTES", "!" ] ) + [
                 "!   @2  No database carries transport for this species, so it",
                 "!       is derived from its lighter isotopologue by mass",
                 "!       scaling: in the dilute gas isotopologues share the",
                 "!       collision integral, so viscosity scales as sqrt(m)",
                 "!       and conductivity as 1/sqrt(m). Only the constant",
                 "!       term D changes; the other three coefficients are the",
                 "!       parent's unchanged. The citation shown on the record",
                 "!       is therefore the parent's measurement, not one of",
                 "!       this species.", "!" ]
        out.append( "!         "
                    + "   ".join( f"{c} from {p}" for c, p in derived ) )

    return ( "\n".join( out ) + "\n!" ) if out else ""

def main():
    here = os.path.dirname( os.path.abspath( __file__ ) )

    ap = argparse.ArgumentParser( description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter )
    ap.add_argument( "--vendor-thermo", required = True )
    ap.add_argument( "--vendor-trans", required = True )
    ap.add_argument( "--species", default = os.path.join( here, "species.txt" ) )
    ap.add_argument( "--out-dir", required = True )
    ap.add_argument( "--junction", type = float, default = 250.0 )
    ap.add_argument( "--samples", type = int, default = 200 )
    ap.add_argument( "--margin", type = float, default = 0.5 )
    ap.add_argument( "--tolerance", type = float, default = 0.02 )
    ap.add_argument( "--dilute-density", type = float, default = 1.0e-4 )
    ap.add_argument( "--nist-transport", action = "store_true" )
    args = ap.parse_args()

    species = read_species( args.species )
    vt = open( args.vendor_thermo, errors = "replace" ).read().split( "\n" )
    vr = open( args.vendor_trans, errors = "replace" ).read().split( "\n" )

    thermo_out, trans_out, missing, fitted = [], [], [], []

    # what the @1 and @2 notes in the two headers will say
    note_thermo, note_trans, note_derived = [], [], []

    print( f"{'species':<16}{'thermo':<8}{'low T interval':<38}transport" )
    print( "-" * 88 )

    for label in species:
        record = read_thermo_record( vt, label )
        if record is None:
            missing.append( label )
            print( f"{label:<16}{'MISSING':<8}" )
            continue

        comment, composition, intervals = record
        note = ""
        fluid = lowT.FLUIDS.get( label )

        if fluid is not None:
            merged, why, detail = merge_thermo( record, fluid, args.junction,
                                                args.samples, args.margin )
            if merged is not None:
                intervals = merged
                note = why
                fitted.append( label )
                note_thermo.append( ( label, ) + detail )

                # a marker, not the text: the note itself is a comment in the
                # header, so this line stays inside the eighty column layout
                comment = with_marks( comment, [ "@1" ] )
            else:
                note = why

        thermo_out.append( comment.rstrip() )
        thermo_out.append( f"{len(intervals):2d}" + composition[ 2: ].rstrip() )
        for head, a, b in intervals:
            thermo_out += [ head.rstrip(), a.rstrip(), b.rstrip() ]

        # transport, either from the vendor file or derived from a lighter isotope
        trec = read_trans_record( vr, label )
        tnote = "--"
        marks = []

        if trec is None and label in ISOTOPES:
            parent = ISOTOPES[ label ]
            prec = read_trans_record( vr, parent )
            precord = read_thermo_record( vt, parent )
            if prec is not None and precord is not None:
                pheader, prows = prec
                trec = ( pheader, scale_transport( prows, pheader,
                                                   molar_mass( composition ),
                                                   molar_mass( precord[ 1 ] ) ) )
                tnote = f"scaled from {parent}"
                marks.append( "@2" )
                note_derived.append( ( label, parent ) )

        if trec is not None:
            theader, trows = trec
            if fluid is not None:
                merged, notes, detail = merge_trans( trec, fluid, args.junction,
                                             args.samples, args.dilute_density,
                                             args.tolerance, args.nist_transport )
                if merged is not None:
                    trows = merged
                    tnote = ",".join( f"{k}:{v}" for k, v in notes.items() ) or "copied"
                    if detail:
                        marks.insert( 0, "@1" )
                        note_trans.append( ( label, detail ) )
            n_v, n_c = len( trows[ "V" ] ), len( trows[ "C" ] )
            trans_out.append( with_marks( f"{label:<34}V{n_v}C{n_c}"
                                          + theader[ 38: ].rstrip(), marks ) )
            trans_out += [ r.rstrip() for r in trows[ "V" ] + trows[ "C" ] ]
            if tnote == "--":
                tnote = f"V{n_v}C{n_c} copied"

        print( f"{label:<16}{'ok':<8}{note:<38}{tnote}" )

    os.makedirs( args.out_dir, exist_ok = True )

    with open( os.path.join( args.out_dir, "thermo.inp" ), "w" ) as f:
        f.write( HEADER_THERMO + note_thermo_block( note_thermo ) + "\nthermo\n"
                 + "\n".join( thermo_out ) + "\nEND PRODUCTS\nEND REACTANTS\n" )

    with open( os.path.join( args.out_dir, "trans.inp" ), "w" ) as f:
        f.write( HEADER_TRANS + note_trans_block( note_trans, note_derived )
                 + "\n" + "\n".join( trans_out ) + "\nend\n" )

    for name in ( "thermo.inp", "trans.inp" ):
        path = os.path.join( args.out_dir, name )
        long = [ n for n, l in enumerate( open( path ), 1 ) if len( l ) > 81 ]
        if long:
            print( f"WARNING: {name} has lines past column 80: {long}" )

    print( f"\n{len(species)-len(missing)} species written, "
           f"{len(fitted)} with a generated low temperature interval" )
    if missing:
        print( f"NOT FOUND in the vendor file: {missing}" )

if __name__ == "__main__":
    main()
