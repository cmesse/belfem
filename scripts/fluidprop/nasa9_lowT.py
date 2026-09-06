#!/usr/bin/env python3
"""
NASA-9 and CEA transport correlations for the low temperature range.

The CEA thermodynamic tables start at 200 K or 300 K, because the underlying fits
were not validated below that. This tool produces an additional interval that runs
from close to the triple point up to the start of the CEA data, so that a species
can be evaluated over the cryogenic range.

Method
------
Ideal gas heat capacity is sampled from CoolProp, which implements the same
reference equations of state that the usual property tables are generated from.
The NASA-9 polynomial is linear in its coefficients, and continuity of value and
slope at the junction is a pair of linear equality constraints, so the fit is an
equality constrained linear least squares problem with a closed form solution. No
iterative optimiser is involved and the result is reproducible.

The two integration constants are then solved exactly, so that enthalpy and
entropy are continuous with the CEA data at the junction. This matters: those
constants carry the integral below the start of the fit, and entropy on the
third law scale is what the equilibrium constants are built from.

Transport coefficients are fitted the same way, to the CEA form
ln( y / scale ) = A ln T + B/T + C/T^2 + D, against the dilute gas limit.

Conventions follow the BELFEM reader, cl_GT_HeatPoly and cl_GT_TransportPoly:

    cp( T ) = Rm * sum a_i * [ T^-2, T^-1, 1, T, T^2, T^3, T^4 ]
    h ( T ) = Rm * ( b1 + sum a_i * [ -T^-1, ln T, T, T^2/2, T^3/3, T^4/4, T^5/5 ] )
    s ( T ) = Rm * ( b2 + sum a_i * [ -T^-2/2, -T^-1, ln T, T, T^2/2, T^3/3, T^4/4 ] )

Usage
-----
    python nasa9_lowT.py --thermo <thermo.inp> --trans <trans.inp> --species N2 O2 ...

Requires numpy and CoolProp.
"""

import argparse
import math
import sys

import numpy as np

try:
    import CoolProp
    import CoolProp.CoolProp as CP
except ImportError:                                     # pragma: no cover
    sys.exit( "CoolProp is required: pip install CoolProp" )

#-------------------------------------------------------------------------------
# molar gas constant, matching src/core/constants.hpp
#-------------------------------------------------------------------------------

RM = 1.380649e-23 * 6.02214076e23

#! CEA transport polynomials carry these scale factors, see cl_GT_TransportPoly
TRANSPORT_SCALE = { "V": 1.0e-7, "C": 1.0e-4 }

#! species label in the CEA tables -> fluid name in CoolProp
FLUIDS = {
    "Ar":  "Argon",   "CH4": "Methane", "H2": "Hydrogen", "He": "Helium",
    "N2":  "Nitrogen","O2":  "Oxygen",  "Ne": "Neon",     "Xe": "Xenon",
    "CO":  "CarbonMonoxide", "CO2": "CarbonDioxide", "D2": "Deuterium",
    "Kr":  "Krypton", "F2":  "Fluorine", "D2O": "HeavyWater",
}

#-------------------------------------------------------------------------------
# NASA-9 basis functions
#-------------------------------------------------------------------------------

def cp_basis( T ):
    """cp / Rm as a linear form in the seven coefficients"""
    return np.array( [ T**-2, T**-1, 1.0, T, T**2, T**3, T**4 ] )

def dcpdT_basis( T ):
    """d( cp / Rm )/dT"""
    return np.array( [ -2.0*T**-3, -T**-2, 0.0, 1.0, 2.0*T, 3.0*T**2, 4.0*T**3 ] )

def h_basis( T ):
    """h / Rm without the enthalpy constant"""
    return np.array( [ -T**-1, math.log( T ), T, T**2/2.0, T**3/3.0,
                       T**4/4.0, T**5/5.0 ] )

def s_basis( T ):
    """s / Rm without the entropy constant"""
    return np.array( [ -0.5*T**-2, -T**-1, math.log( T ), T, T**2/2.0,
                       T**3/3.0, T**4/4.0 ] )

#-------------------------------------------------------------------------------
# CEA transport basis
#-------------------------------------------------------------------------------

def tr_basis( T ):
    return np.array( [ math.log( T ), 1.0/T, T**-2, 1.0 ] )

def dtr_basis( T ):
    return np.array( [ 1.0/T, -T**-2, -2.0*T**-3, 0.0 ] )

#-------------------------------------------------------------------------------

def constrained_lstsq( A, y, C, d ):
    """
    Minimise || A x - y ||^2 subject to C x = d.

    Solved through the Karush-Kuhn-Tucker system, which is exact and has no
    starting guess, so the result does not depend on an optimiser.
    """
    n = A.shape[ 1 ]
    m = C.shape[ 0 ]

    K = np.zeros( ( n + m, n + m ) )
    K[ :n, :n ] = 2.0 * A.T @ A
    K[ :n, n: ] = C.T
    K[ n:, :n ] = C

    r = np.concatenate( [ 2.0 * A.T @ y, d ] )

    return np.linalg.solve( K, r )[ :n ]

#-------------------------------------------------------------------------------
# readers for the CEA input files
#-------------------------------------------------------------------------------

def _f( text ):
    """
    CEA writes exponents with D in the thermo file, and in the transport file it
    writes a positive exponent with a space where the plus sign would be, as in
    "0.61205763E 00". Both forms are normalised here.
    """
    text = text.strip().replace( "D", "E" ).replace( "d", "e" )
    if not text:
        return 0.0
    if "E " in text:
        text = text.replace( "E ", "E+" )
    return float( text )

def read_thermo( path, label ):
    """
    Return the lowest temperature interval of a species as
    ( Tmin, Tmax, coefficients, enthalpy constant, entropy constant ).
    """
    lines = open( path, errors = "replace" ).read().split( "\n" )

    for i, line in enumerate( lines ):
        if not line.strip() or line[ 0 ] in " !-":
            continue
        if line.split()[ 0 ] != label:
            continue

        n_intervals = int( lines[ i + 1 ][ 0:2 ] )
        if n_intervals < 1:
            raise ValueError( f"{label} carries no temperature interval" )

        head = lines[ i + 2 ]
        t_min = float( head[ 0:11 ] )
        t_max = float( head[ 11:22 ] )

        a = [ _f( lines[ i + 3 ][ k:k+16 ] ) for k in range( 0, 80, 16 ) ]
        a += [ _f( lines[ i + 4 ][ k:k+16 ] ) for k in range( 0, 32, 16 ) ]

        b1 = _f( lines[ i + 4 ][ 48:64 ] )
        b2 = _f( lines[ i + 4 ][ 64:80 ] )

        # The composition line carries the elements, the molar mass and the
        # formation enthalpy, and the interval header carries H(298)-H(0).
        # The reader needs all of them, so they are handed back verbatim and
        # reused rather than reconstructed.
        return ( t_min, t_max, np.array( a[ :7 ] ), b1, b2,
                 lines[ i + 1 ], head )

    raise KeyError( f"{label} not found in {path}" )

def read_trans( path, label ):
    """
    Return { "V": [ ( Tmin, Tmax, coefficients ), ... ], "C": [ ... ] } for a
    pure species. Interaction pairs, which carry a second name, are skipped.
    """
    out = { "V": [], "C": [] }
    lines = open( path, errors = "replace" ).read().split( "\n" )

    def is_header( line ):
        """a record header carries V<n>C<m> at column 34"""
        return ( len( line ) > 38
                 and line[ 34 ] == "V" and line[ 36 ] == "C"
                 and line[ 35 ].isdigit() and line[ 37 ].isdigit() )

    i = 0
    while i < len( lines ):
        line = lines[ i ]

        if not is_header( line ):
            i += 1
            continue

        n_v = int( line[ 35 ] )
        n_c = int( line[ 37 ] )

        name = line[ 0:15 ].strip()
        partner = line[ 15:30 ].strip()

        if name == label and not partner:
            for k in range( n_v + n_c ):
                row = lines[ i + 1 + k ]
                kind = row[ 1:2 ]
                out[ kind ].append( (
                    float( row[ 2:11 ] ), float( row[ 11:20 ] ),
                    np.array( [ _f( row[ 20+15*j : 35+15*j ] ) for j in range( 4 ) ] ) ) )
            return out

        i += 1 + n_v + n_c

    return out

#-------------------------------------------------------------------------------
# evaluation of an existing CEA interval, used for the junction conditions
#-------------------------------------------------------------------------------

def cea_cp( a, T ):        return RM * float( cp_basis( T ) @ a )
def cea_dcpdT( a, T ):     return RM * float( dcpdT_basis( T ) @ a )
def cea_h( a, b1, T ):     return RM * ( b1 + float( h_basis( T ) @ a ) )
def cea_s( a, b2, T ):     return RM * ( b2 + float( s_basis( T ) @ a ) )

def cea_transport( c, T, kind ):
    return math.exp( float( tr_basis( T ) @ c ) ) * TRANSPORT_SCALE[ kind ]

#-------------------------------------------------------------------------------
# the fits
#-------------------------------------------------------------------------------

def ideal_cp( fluid, T ):
    """
    Ideal gas heat capacity in J/mol/K.

    Cp0 comes from the ideal gas part of the Helmholtz energy and is a function of
    temperature alone; evaluated across four decades of pressure it does not move
    at all. CoolProp nevertheless refuses the call when the state given by
    temperature and pressure lies below the melting line, which happens near the
    triple point of deuterium and hydrogen. Since the pressure is immaterial to
    the answer, a lower one is tried until the guard is satisfied.
    """
    last = None
    for pressure in ( 1.0e5, 1.0e4, 1.0e3, 1.0e2, 1.0e1 ):
        try:
            return CP.PropsSI( "Cp0molar", "T", float( T ), "P", pressure, fluid )
        except Exception as error:
            last = error
    raise RuntimeError( f"ideal gas cp unavailable for {fluid} at {T} K: {last}" )

def _fit_heat_global( fluid, edges, a_ref, b1_ref, b2_ref, n_samples ):
    """
    Fit every interval at once, for the same reason the transport fit does: the
    error at a shared edge would otherwise propagate downward. All coefficients
    are solved together, with continuity of cp and its slope at each internal
    edge, and the match to the reference polynomial at the junction, entering as
    linear equality constraints.
    """
    n = len( edges ) - 1
    rows, rhs = [], []

    for k in range( n ):
        T = np.linspace( edges[ k ], edges[ k + 1 ], max( 24, n_samples // n ) )

        # ideal gas heat capacity. Cp0 is the quantity the NASA-9 form
        # represents; the real gas cp at a finite pressure would double count
        # the departure.
        y = np.array( [ ideal_cp( fluid, t ) for t in T ] ) / RM

        for t, v in zip( T, y ):
            row = np.zeros( 7 * n )
            row[ 7*k : 7*k+7 ] = cp_basis( float( t ) )
            rows.append( row )
            rhs.append( v )

    A = np.array( rows )
    y = np.array( rhs )

    cons, d = [], []

    # match the reference polynomial at the junction, in value and slope
    for basis in ( cp_basis, dcpdT_basis ):
        row = np.zeros( 7 * n )
        row[ 7*(n-1) : 7*n ] = basis( edges[ -1 ] )
        cons.append( row )
        d.append( float( basis( edges[ -1 ] ) @ a_ref ) )

    # continuity between neighbouring intervals
    for k in range( 1, n ):
        for basis in ( cp_basis, dcpdT_basis ):
            row = np.zeros( 7 * n )
            row[ 7*(k-1) : 7*k   ] =  basis( edges[ k ] )
            row[ 7*k     : 7*k+7 ] = -basis( edges[ k ] )
            cons.append( row )
            d.append( 0.0 )

    x = constrained_lstsq( A, y, np.array( cons ), np.array( d ) )

    residual = float( np.max( np.abs( A @ x - y ) / np.abs( y ) ) )

    # The integration constants are not part of the least squares problem: once
    # cp is fixed, exact continuity of h and s determines them. They are chained
    # downward from the junction, where the reference constants are known.
    blocks = []
    a_above, b1, b2 = a_ref, b1_ref, b2_ref

    for k in range( n - 1, -1, -1 ):
        a    = x[ 7*k : 7*k+7 ]
        edge = edges[ k + 1 ]

        b1 = b1 + float( h_basis( edge ) @ a_above ) - float( h_basis( edge ) @ a )
        b2 = b2 + float( s_basis( edge ) @ a_above ) - float( s_basis( edge ) @ a )

        blocks.insert( 0, ( float( edges[ k ] ), float( edges[ k + 1 ] ),
                            a, b1, b2 ) )
        a_above = a

    return blocks, residual

def fit_heat( fluid, t_low, t_junction, a_ref, b1_ref, b2_ref, n_samples,
              tolerance = 0.01, max_intervals = 3 ):
    """
    Fit the NASA-9 form over [ t_low, t_junction ], adding intervals until the
    tolerance is met, and return [ ( t_lo, t_hi, a, b1, b2 ), ... ] with the
    worst relative error in cp.

    One interval is enough for most species. It is not enough where cp is not
    monotonic over the range: normal deuterium overshoots its classical limit
    near 100 K and comes back down, a shape the seven term form cannot hold
    across the whole span at once.

    Returns the best result found, so a species that cannot reach the tolerance
    still produces its closest fit together with an honest residual.
    """
    best = ( None, None )

    for n in range( 1, max_intervals + 1 ):
        edges = np.exp( np.linspace( math.log( t_low ),
                                     math.log( t_junction ), n + 1 ) )

        # exp( log( x ) ) is not exactly x. The ends are the triple point and
        # the junction, and the junction is where the constraint against the
        # vendor polynomial is imposed, so they are restored exactly.
        edges[ 0 ], edges[ -1 ] = t_low, t_junction

        blocks, residual = _fit_heat_global( fluid, edges, a_ref, b1_ref,
                                             b2_ref, n_samples )
        if best[ 0 ] is None or residual < best[ 1 ]:
            best = ( blocks, residual )
        if residual <= tolerance:
            break

    return best

#! per species dilute gas tables pulled from NIST, filled on demand
_NIST_TRANSPORT = {}

def _nist_transport( fluid, kind, t_a, t_b, dilute_pressure ):
    """
    Dilute gas transport from the NIST WebBook, for species CoolProp does not
    model. Fetched once per fluid and reused.
    """
    key = ( fluid, kind )
    if key not in _NIST_TRANSPORT:
        prop = "viscosity" if kind == "V" else "conductivity"
        table = fetch_nist( fluid, t_a, t_b, n_points = 60,
                            pressure = dilute_pressure )
        _NIST_TRANSPORT[ key ] = ( table or {} ).get( prop, [] )

    return [ ( t, v ) for t, v in _NIST_TRANSPORT[ key ]
             if t_a - 1.0e-6 <= t <= t_b + 1.0e-6 and v > 0.0 ]

def _sample_transport( fluid, kind, t_a, t_b, n_samples, dilute_density,
                       allow_nist = False, dilute_pressure = 0.001 ):
    """
    Dilute gas samples of ln( y / scale ) over [ t_a, t_b ].

    CoolProp is preferred, since it is reproducible and needs no network. Where
    it carries no transport model, and only then, the samples come from NIST at
    a low pressure instead.
    """
    prop = "viscosity" if kind == "V" else "conductivity"
    scale = TRANSPORT_SCALE[ kind ]

    T, y = [], []
    for t in np.linspace( t_a, t_b, n_samples ):
        try:
            v = CP.PropsSI( prop, "T", float( t ), "D", dilute_density, fluid )
        except Exception:
            continue
        if v > 0.0:
            T.append( float( t ) )
            y.append( math.log( v / scale ) )

    if not T and allow_nist:
        for t, v in _nist_transport( fluid, kind, t_a, t_b, dilute_pressure ):
            T.append( t )
            y.append( math.log( v / scale ) )

    return np.array( T ), np.array( y )

def _fit_transport_global( fluid, kind, edges, c_ref, n_samples, dilute_density,
                           allow_nist = False ):
    """
    Fit every interval at once.

    Fitting the intervals one after another, each constrained to the one above,
    lets the error at a shared edge propagate downward and can be worse than a
    single interval. Here all coefficients are solved together: the data of every
    interval enters one least squares problem, and continuity of value and slope
    at each internal edge, plus the match to the reference polynomial at the top,
    enter as linear equality constraints.
    """
    n = len( edges ) - 1
    rows, rhs = [], []

    for k in range( n ):
        T, y = _sample_transport( fluid, kind, edges[ k ], edges[ k + 1 ],
                                  max( 24, n_samples // n ), dilute_density,
                                  allow_nist )
        if len( T ) < 6:
            return None, None
        for t, v in zip( T, y ):
            row = np.zeros( 4 * n )
            row[ 4*k : 4*k+4 ] = tr_basis( t )
            rows.append( row )
            rhs.append( v )

    A = np.array( rows )
    y = np.array( rhs )

    cons, d = [], []

    # match the reference polynomial at the junction, in value and slope
    for basis in ( tr_basis, dtr_basis ):
        row = np.zeros( 4 * n )
        row[ 4*(n-1) : 4*n ] = basis( edges[ -1 ] )
        cons.append( row )
        d.append( float( basis( edges[ -1 ] ) @ c_ref ) )

    # continuity between neighbouring intervals
    for k in range( 1, n ):
        for basis in ( tr_basis, dtr_basis ):
            row = np.zeros( 4 * n )
            row[ 4*(k-1) : 4*k ]   =  basis( edges[ k ] )
            row[ 4*k     : 4*k+4 ] = -basis( edges[ k ] )
            cons.append( row )
            d.append( 0.0 )

    x = constrained_lstsq( A, y, np.array( cons ), np.array( d ) )

    residual = float( np.max( np.abs( np.exp( A @ x ) - np.exp( y ) )
                              / np.exp( y ) ) )

    blocks = [ ( float( edges[ k ] ), float( edges[ k + 1 ] ),
                 x[ 4*k : 4*k+4 ] ) for k in range( n ) ]
    return blocks, residual

def fit_transport( fluid, kind, t_low, t_junction, c_ref, n_samples,
                   dilute_density, tolerance = 0.02, max_intervals = 3,
                   allow_nist = False ):
    """
    Fit the CEA transport form over [ t_low, t_junction ], adding intervals until
    the tolerance is met. The four parameter form cannot span a very wide range,
    which is why the CEA files carry several intervals per species themselves.

    Returns the best result found, so a species that cannot reach the tolerance
    still produces its closest fit together with an honest residual.
    """
    best = ( None, None )

    for n in range( 1, max_intervals + 1 ):
        edges = np.exp( np.linspace( math.log( t_low ),
                                     math.log( t_junction ), n + 1 ) )

        # exp( log( x ) ) is not exactly x. The ends are the triple point and
        # the junction, and the junction is where the constraint against the
        # vendor polynomial is imposed, so they are restored exactly.
        edges[ 0 ], edges[ -1 ] = t_low, t_junction

        blocks, residual = _fit_transport_global( fluid, kind, edges, c_ref,
                                                  n_samples, dilute_density,
                                                  allow_nist )
        if blocks is None:
            break
        if best[ 0 ] is None or residual < best[ 1 ]:
            best = ( blocks, residual )
        if residual <= tolerance:
            break

    return best

def check_thermo_continuity( a_new, b1, b2, a_ref, b1_ref, b2_ref, t_junction ):
    """
    Relative jump in cp, its slope, h and s at the junction. The fit is only
    useful if these are at round off level: a step in h or s would corrupt the
    absolute scale that the equilibrium constants are built on.
    """
    # Each jump is scaled by a physical magnitude rather than by the reference
    # value itself. For a monatomic gas cp is constant, so dcp/dT is numerically
    # zero and scaling by it would turn a negligible absolute error into a large
    # relative one. cp/T is the natural scale for that derivative.
    cp_ref = cea_cp( a_ref, t_junction )

    def rel( new, ref, scale ):
        return abs( new - ref ) / abs( scale )

    return {
        "cp":     rel( cea_cp( a_new, t_junction ),
                       cp_ref, cp_ref ),
        "dcp/dT": rel( cea_dcpdT( a_new, t_junction ),
                       cea_dcpdT( a_ref, t_junction ), cp_ref / t_junction ),
        "h":      rel( cea_h( a_new, b1, t_junction ),
                       cea_h( a_ref, b1_ref, t_junction ), cp_ref * t_junction ),
        "s":      rel( cea_s( a_new, b2, t_junction ),
                       cea_s( a_ref, b2_ref, t_junction ), cp_ref ),
    }

def check_transport_continuity( blocks, c_ref, t_junction ):
    """Relative jump at the junction and at every internal edge."""
    worst = 0.0

    for basis in ( tr_basis, dtr_basis ):
        top = blocks[ -1 ][ 2 ]
        lhs = float( basis( t_junction ) @ top )
        rhs = float( basis( t_junction ) @ c_ref )
        worst = max( worst, abs( lhs - rhs ) / max( abs( rhs ), 1.0e-12 ) )

    for k in range( 1, len( blocks ) ):
        t_edge = blocks[ k ][ 0 ]
        for basis in ( tr_basis, dtr_basis ):
            lhs = float( basis( t_edge ) @ blocks[ k - 1 ][ 2 ] )
            rhs = float( basis( t_edge ) @ blocks[ k ][ 2 ] )
            worst = max( worst, abs( lhs - rhs ) / max( abs( rhs ), 1.0e-12 ) )

    return worst

def fetch_nist( fluid, t_low, t_high, n_points = 40, timeout = 60,
                pressure = 1.0 ):
    """
    Fetch an isobaric table from the NIST Chemistry WebBook and return
    it as { property: [ ( T, value ), ... ] }, vapour rows only.

    Used both to validate CoolProp and to supply transport data for species that
    CoolProp does not model, such as neon and xenon.

    For transport, pass a low pressure. The CEA transport form describes the
    dilute gas, and at 1 bar the density contribution is not negligible near the
    boiling point: it reaches 6 % in the conductivity of hydrogen at 25 K and
    about 2 % for neon at 30 K.
    """
    import urllib.parse
    import urllib.request

    cas = CP.get_fluid_param_string( fluid, "CAS" ).replace( "-", "" )
    step = max( ( t_high - t_low ) / max( n_points - 1, 1 ), 0.05 )

    query = {
        "Action": "Data", "Wide": "on", "ID": f"C{cas}", "Type": "IsoBar",
        "Digits": "8", "P": f"{pressure:g}", "TLow": f"{t_low:.4f}", "THigh": f"{t_high:.4f}",
        "TInc": f"{step:.4f}", "RefState": "DEF", "TUnit": "K", "PUnit": "bar",
        "DUnit": "kg/m3", "HUnit": "kJ/kg", "WUnit": "m/s",
        "VisUnit": "uPa*s", "STUnit": "N/m",
    }
    url = "https://webbook.nist.gov/cgi/fluid.cgi?" + urllib.parse.urlencode( query )

    try:
        request = urllib.request.Request(
            url, headers = { "User-Agent": "belfem-fluidprop-validation/1.0" } )
        with urllib.request.urlopen( request, timeout = timeout ) as response:
            text = response.read().decode( "utf-8", errors = "replace" )
    except Exception:
        return None

    lines = [ l for l in text.strip().split( "\n" ) if l.strip() ]
    if len( lines ) < 2 or not lines[ 0 ].lower().startswith( "temperature" ):
        return None

    COLUMNS = { "h": 5, "s": 6, "cp": 8, "viscosity": 11, "conductivity": 12 }
    UNITS   = { "h": 1.0e3, "s": 1.0e3, "cp": 1.0e3,
                "viscosity": 1.0e-6, "conductivity": 1.0 }
    PHASE   = 13

    out = { key: [] for key in COLUMNS }
    for line in lines[ 1: ]:
        field = line.split( "\t" )
        if len( field ) <= PHASE:
            continue
        # a state given by temperature and pressure is not unique on the
        # saturation line, so only vapour rows can be used
        if field[ PHASE ].strip().lower() != "vapor":
            continue
        try:
            T = float( field[ 0 ] )
        except ValueError:
            continue
        for key, column in COLUMNS.items():
            try:
                out[ key ].append( ( T, float( field[ column ] ) * UNITS[ key ] ) )
            except ValueError:
                pass

    return { k: v for k, v in out.items() if v } or None

def validate_against_nist( fluid, t_low, t_high, n_points = 9, timeout = 45 ):
    """
    Compare CoolProp with the NIST WebBook and return the largest relative
    difference per property.

    CoolProp is the sampling engine because it is scriptable, reproducible and
    can be pinned to a version. NIST is the authority where the two differ. For
    most fluids they implement the same reference correlations and agree to the
    printed digits; this makes that agreement visible per species rather than
    assumed.
    """
    table = fetch_nist( fluid, t_low, t_high, n_points, timeout )
    if table is None:
        return None

    NAMES = { "h": "Hmass", "s": "Smass", "cp": "Cpmass",
              "viscosity": "viscosity", "conductivity": "conductivity" }

    worst = {}
    for key, series in table.items():
        for T, reference in series:
            try:
                ours = CP.PropsSI( NAMES[ key ], "T", T, "P", 1.0e5, fluid )
            except Exception:
                continue
            if abs( reference ) > 0.0:
                worst[ key ] = max( worst.get( key, 0.0 ),
                                    abs( ours - reference ) / abs( reference ) )
    return worst or None

#-------------------------------------------------------------------------------
# emitters, in the fixed column layout the readers expect
#-------------------------------------------------------------------------------

def _d( value ):
    """sixteen character field in the D exponent notation CEA uses"""
    text = f"{value: .9E}"
    mantissa, exponent = text.split( "E" )
    return f"{mantissa}D{exponent}"

def emit_thermo( label, blocks, composition_line, reference_head, comment = "" ):
    """
    Emit one record in the fixed column CEA layout, from the intervals returned
    by fit_heat.

    The composition line is taken from the source record with only the interval
    count rewritten, so the elements, molar mass and formation enthalpy are
    carried over unchanged. The tail of the reference interval header, which
    holds H(298)-H(0), is carried over the same way.
    """
    out = [ f"{label:<18}{comment}",
            f"{len( blocks ):2d}" + composition_line[ 2: ].rstrip() ]

    tail = reference_head[ 65: ].rstrip()

    for t_low, t_high, a, b1, b2 in blocks:
        out.append( f"{t_low:11.3f}{t_high:11.3f}"
                    "7 -2.0 -1.0  0.0  1.0  2.0  3.0  4.0  0.0"
                    + tail.rjust( 65 - 63 + len( tail ) ) )
        out.append( "".join( _d( a[ k ] ) for k in range( 5 ) ) )
        out.append( "".join( _d( a[ k ] ) for k in range( 5, 7 ) )
                    + " " * 16 + _d( b1 ) + _d( b2 ) )

    return "\n".join( out )

def check_thermo_blocks( blocks, a_ref, b1_ref, b2_ref ):
    """
    Worst relative jump over every edge of a generated set of intervals: between
    neighbouring generated intervals, and between the topmost one and the
    reference polynomial it hands over to.
    """
    worst = {}
    a_above, b1_above, b2_above = a_ref, b1_ref, b2_ref

    for lo, hi, a, b1, b2 in reversed( blocks ):
        for key, value in check_thermo_continuity( a, b1, b2, a_above,
                                                   b1_above, b2_above,
                                                   hi ).items():
            worst[ key ] = max( worst.get( key, 0.0 ), value )
        a_above, b1_above, b2_above = a, b1, b2

    return worst

def _e( value ):
    """
    Fifteen character field in the notation the CEA transport file uses, where a
    positive exponent is written with a space instead of a plus sign.
    """
    text = f"{value: .8E}"
    mantissa, exponent = text.split( "E" )
    sign = " " if exponent[ 0 ] == "+" else "-"
    return f"{mantissa}E{sign}{exponent[ 1: ]}"

def emit_transport( label, blocks, source = "" ):
    n_v = len( blocks.get( "V", [] ) )
    n_c = len( blocks.get( "C", [] ) )

    # the trailing source text also keeps the header longer than the fixed
    # columns, which is what distinguishes a header from a polynomial row
    out = [ f"{label:<34}V{n_v}C{n_c}  {source}".rstrip() ]

    for kind in ( "V", "C" ):
        for t_low, t_high, c in blocks.get( kind, [] ):
            out.append( f" {kind}{t_low:9.1f}{t_high:9.1f}"
                        + "".join( _e( v ) for v in c ) )
    return "\n".join( out )

#-------------------------------------------------------------------------------

def main():
    ap = argparse.ArgumentParser( description = __doc__,
             formatter_class = argparse.RawDescriptionHelpFormatter )
    ap.add_argument( "--thermo", required = True, help = "path to thermo.inp" )
    ap.add_argument( "--trans", help = "path to trans.inp" )
    ap.add_argument( "--species", nargs = "+", required = True )
    ap.add_argument( "--samples", type = int, default = 200 )
    ap.add_argument( "--junction", type = float, default = 250.0,
                     help = "temperature where the fit hands over to the CEA data. "
                            "The default keeps 273.15 K and above on pure CEA. Pass 0 "
                            "to hand over at the start of the CEA interval instead." )
    ap.add_argument( "--nist-transport", action = "store_true",
                     help = "for species CoolProp cannot model, take the dilute gas "
                            "transport from the NIST WebBook instead" )
    ap.add_argument( "--validate-nist", action = "store_true",
                     help = "check CoolProp against the NIST WebBook per species" )
    ap.add_argument( "--margin", type = float, default = 0.5,
                     help = "start this many K above the triple point" )
    ap.add_argument( "--tolerance", type = float, default = 0.02,
                     help = "relative transport residual before splitting the range" )
    ap.add_argument( "--dilute-density", type = float, default = 1.0e-4,
                     help = "density in kg/m3 used for the dilute gas limit" )
    ap.add_argument( "--out-thermo", help = "write the thermo records here" )
    ap.add_argument( "--out-trans", help = "write the transport records here" )
    args = ap.parse_args()

    thermo_records, trans_records = [], []

    print( f"{'species':<8}{'fluid':<16}{'range [K]':>18}"
           f"{'cp resid':>11}{'visc':>10}{'cond':>10}{'jump':>11}" )
    print( "-" * 84 )

    for label in args.species:
        fluid = FLUIDS.get( label )
        if fluid is None:
            print( f"{label:<8}no CoolProp fluid is mapped to this label" )
            continue

        try:
            ( t_ref_min, _, a_ref, b1_ref, b2_ref,
              comp_line, ref_head ) = read_thermo( args.thermo, label )
        except KeyError as exc:
            print( f"{label:<8}{exc}" )
            continue

        t_triple = CP.PropsSI( "Ttriple", fluid )
        t_low = t_triple + args.margin

        # Hand over where the CEA data is still reliable. Measured against the
        # reference equations of state, the CEA polynomials agree to better than
        # 0.1 % down to about 200 K and degrade quickly below that, so a junction
        # at 250 K is comfortably inside their useful range even when the file
        # declares a higher lower bound.
        t_junction = t_ref_min if args.junction <= 0.0 else args.junction
        if t_junction < t_ref_min:
            print( f"{label:<8}note: handing over at {t_junction:.1f} K, below the "
                   f"declared start of the CEA interval at {t_ref_min:.1f} K" )

        t_ref_min = t_junction

        if t_low >= t_ref_min:
            print( f"{label:<8}{fluid:<16}"
                   f"triple point {t_triple:.2f} K is above the junction, skipped" )
            continue

        blocks, resid = fit_heat( fluid, t_low, t_ref_min,
                                  a_ref, b1_ref, b2_ref, args.samples )

        jump = check_thermo_blocks( blocks, a_ref, b1_ref, b2_ref )
        worst_jump = max( jump.values() )
        if worst_jump > 1.0e-9:
            print( f"{label:<8}CONTINUITY FAILED at an interval edge: "
                   + ", ".join( f"{k} {v:.2e}" for k, v in jump.items() ) )

        thermo_records.append(
            emit_thermo( label, blocks, comp_line, ref_head,
                         f"low temperature fit, CoolProp {CoolProp.__version__}" ) )

        v_txt = c_txt = "--"
        if args.trans:
            ref = read_trans( args.trans, label )
            blocks = {}
            for kind in ( "V", "C" ):
                if not ref.get( kind ):
                    continue
                # hand over at the same temperature as the heat capacity, using
                # whichever reference polynomial covers it
                t_ref = t_junction
                c_ref = ref[ kind ][ 0 ][ 2 ]
                for lo, hi, c in ref[ kind ]:
                    if lo <= t_ref <= hi:
                        c_ref = c
                        break
                if t_low >= t_ref:
                    continue
                got, r = fit_transport( fluid, kind, t_low, t_ref, c_ref,
                                        args.samples, args.dilute_density,
                                        args.tolerance,
                                        allow_nist = args.nist_transport )
                if got is None:
                    continue
                blocks[ kind ] = got
                tj = check_transport_continuity( got, c_ref, t_ref )
                if tj > 1.0e-9:
                    print( f"{label:<8}{kind} transport continuity failed: {tj:.2e}" )
                text = f"{100*r:.2f}%" + ( f"/{len(got)}" if len( got ) > 1 else "" )
                if kind == "V":
                    v_txt = text
                else:
                    c_txt = text
            if blocks:
                trans_records.append( emit_transport(
                    label, blocks,
                    f'COOLPROP {CoolProp.__version__} LOW T FIT' ) )

        print( f"{label:<8}{fluid:<16}{t_low:8.2f} to {t_ref_min:6.1f}"
               f"{100*resid:10.2f}%{v_txt:>10}{c_txt:>10}{worst_jump:>11.1e}" )

        if args.validate_nist:
            # stay in the vapour region, so the comparison is well posed
            try:
                t_boil = CP.PropsSI( "T", "P", 1.0e5, "Q", 1, fluid ) + 1.0
            except Exception:
                t_boil = t_low
            agreement = validate_against_nist( fluid, max( t_low, t_boil ), t_ref_min )
            if agreement is None:
                print( f"{'':<8}NIST check unavailable for this species" )
            else:
                print( f"{'':<8}NIST agreement: "
                       + ", ".join( f"{k} {100*v:.4f}%"
                                    for k, v in sorted( agreement.items() ) ) )

    if args.out_thermo and thermo_records:
        with open( args.out_thermo, "w" ) as f:
            f.write( "thermo\n" + "\n".join( thermo_records ) + "\nend\n" )
        print( f"\nwrote {args.out_thermo}" )

    if args.out_trans and trans_records:
        with open( args.out_trans, "w" ) as f:
            f.write( "transport property coefficients\n"
                     + "\n".join( trans_records ) + "\nend\n" )
        print( f"wrote {args.out_trans}" )

if __name__ == "__main__":
    main()
