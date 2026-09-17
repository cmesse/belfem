% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Shorthand used in the C++ tables: xi2 = 2*xi, xi4, xi8, xi16, and the
% corresponding forms for eta, zeta and tau. The two zero-tests run this
% after tet10_generate.m.
xi2 = xi+xi ;
xi4 = xi2+xi2 ;
xi8 = xi4+xi4 ;
xi16 = xi8+xi8 ;

eta2 = eta+eta ;
eta4 = eta2+eta2 ;
eta8 = eta4+eta4 ;
eta16 = eta8+eta8 ;

zeta2 = zeta+zeta ;
zeta4 = zeta2+zeta2 ;
zeta8 = zeta4+zeta4 ;
zeta16 = zeta8+zeta8 ;

tau = 1-xi-eta-zeta ;
tau2 = tau + tau ;
tau4 = tau2+tau2;
tau8 = tau4+tau4;