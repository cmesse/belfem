% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Asserting driver for penta18.m:
%   matlab -batch check_penta18
% exits 0 only if every entry of the transcribed d2NdXi2 equals the symbolic
% second derivative of the C++ shape functions.
penta18 ;
assert( all( isAlways( X( : ) == 0 ) ), 'penta18: d2NdXi2 differs from the symbolic second derivatives' );
disp( 'check_penta18: all 108 entries match' );
