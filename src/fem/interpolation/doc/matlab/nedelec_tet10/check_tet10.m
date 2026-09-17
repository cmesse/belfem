% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Asserting driver for the TET10 zero-tests. Runs tet10_function.m and
% tet10_derivatives.m and fails on any nonzero residual, so that
%   matlab -batch check_tet10
% exits 0 only if the generator reproduces every transcribed table.
tet10_function ;
r = [ expand( g - mG ); expand( h - mH ); ...
      expand( u( act ) - mU ); expand( v( act ) - mV ); expand( w( act ) - mW ) ];
assert( all( isAlways( r == 0 ) ), 'tet10_function: a value table differs from the generator' );
tet10_derivatives ;
r = [ expand( diff( g, xi ) - mGxi ); expand( diff( g, eta ) - mGeta ); expand( diff( g, zeta ) - mGzeta ); ...
      expand( diff( h, xi ) - mHxi ); expand( diff( h, eta ) - mHeta ); expand( diff( h, zeta ) - mHzeta ); ...
      expand( diff( u( act ), xi ) - mUxi ); expand( diff( u( act ), eta ) - mUeta ); expand( diff( u( act ), zeta ) - mUzeta ); ...
      expand( diff( v( act ), xi ) - mVxi ); expand( diff( v( act ), eta ) - mVeta ); expand( diff( v( act ), zeta ) - mVzeta ); ...
      expand( diff( w( act ), xi ) - mWxi ); expand( diff( w( act ), eta ) - mWeta ); expand( diff( w( act ), zeta ) - mWzeta ) ];
assert( all( isAlways( r == 0 ) ), 'tet10_derivatives: a derivative table differs from the generator' );
disp( 'check_tet10: 24 value rows and 15 derivative tables match the generator' );
