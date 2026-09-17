% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Derivation note: the two-point-per-edge ( a, b ) ansatz for the six
% second-order triangle edge functions; fragment.m expands this form. It
% stops at the ansatz and prints nothing; nothing here is a C++ table.
% Needs the Symbolic Math Toolbox.
clear all ;

syms a b xi s eta real ;

%a = 0.5 * ( 1 - s );
%b = 0.5 * ( 1 + s );

syms dXi dEta real ;

zeta = 1 - xi - eta ;
dZeta = -dXi - dEta ;

N1 = ( a * xi + b ) * dEta + ( b * eta + a ) * dXi ;
N2 = ( b * xi + a ) * dEta + ( a * eta + b ) * dXi ;


N3 = ( a * eta + b ) * dZeta + ( b * zeta + a ) * dEta ;
N4 = ( b * eta + a ) * dZeta + ( a * zeta + b ) * dEta ;

N5 = ( a * zeta + b ) * dXi + ( b * xi + a ) * dZeta ;
N6 = ( b * zeta + a ) * dXi + ( a * xi + b ) * dZeta ;