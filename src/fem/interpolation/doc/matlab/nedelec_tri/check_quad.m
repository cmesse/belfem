% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Derivation note: the Jacobians of the [ 0, 1 ]^2 ( alpha, beta ) and the
% [ -1, 1 ]^2 ( xi, eta ) parametrizations of a bilinear quad differ by a
% constant factor induced by alpha = 0.5*( 1 + xi ). This is the remap used
% by the facet integration ( fn_IF_initialize_integration_points_on_facet.cpp ).
% Needs the Symbolic Math Toolbox.
clear all;

syms x1 x2 x3 x4 real ;
syms y1 y2 y3 y4 real ;

x = [x1; x2; x3; x4 ];
y = [y1; y2; y3; y4 ];

syms alpha beta xi eta real ;

psi = [(1-alpha)*(1-beta) alpha*(1-beta) alpha*beta (1-alpha)*beta];
phi = [(1-xi)*(1-eta) (1+xi)*(1-eta) (1+xi)*(1+eta) (1-xi)*(1+eta)]*0.25;

Ja = vpa(zeros(2));
Jb = vpa(zeros(2));


Ja(1,1) = diff(psi,alpha)*x;
Ja(1,2) = diff(psi,alpha)*y;
Ja(2,1) = diff(psi,beta)*x;
Ja(2,2) = diff(psi,beta)*y;

alpha = 0.5*(1+xi);
beta = 0.5*(1+eta);

Ja = eval(Ja)

Jb(1,1) = diff(phi,xi)*x;
Jb(1,2) = diff(phi,xi)*y;
Jb(2,1) = diff(phi,eta)*x;
Jb(2,2) = diff(phi,eta)*y;

simplify(expand(Ja(1,1)/Jb(1,1)))
