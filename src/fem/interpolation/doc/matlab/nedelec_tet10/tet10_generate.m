% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Generator of the Nedelec edge and face tables in cl_EF_TET10.cpp
% precompute(): mG, mH for the twelve edge functions, mU, mV, mW for the
% face functions ( the Lagrange tables mNxi, mNeta, mNzeta there come from
% tet10_lagrange.m ). Run first; tet10_function.m and tet10_derivatives.m
% check the C++ tables against g, h, u, v, w. Needs the Symbolic Math
% Toolbox.
clear all ;
syms xi eta zeta real ;
tau = 1-xi-eta-zeta ;

% EXODUS node order: node 1 carries zeta, node 2 carries eta. The naive
% triangle-extension map [xi; eta; zeta; tau] is left-handed ( detJ < 0 )
% and yields tables that fail edge conformity. Must match defelement.m.
lambda = [xi; zeta; eta; tau ];

%% Edge functions for edges from i to j
% each edge has two functions:
% Ek: gk*nabla_j - hk*nabla_i

g = vpa(zeros(12,1));
h = vpa(zeros(12,1));

E = [1 2; 2 3; 3 1; 1 4; 2 4; 3 4 ];
c = 0 ;
for e =1:6
    i = E(e,1);
    j = E(e,2);

    c = c + 1 ;
    g(c) = simplify(4*lambda(i)*(2*lambda(i)-1));
    h(c) = simplify(2*lambda(j)*(4*lambda(i)-1));

    c = c + 1 ;
    g(c) = simplify(2*lambda(i)*(4*lambda(j)-1));
    h(c) = simplify(4*lambda(j)*(2*lambda(j)-1));
end

%% faces work as follows:
% each face has three functions, where F1+F2+F3=0:
% Fk : u * nabla_i + v * nabla_j + w*nabla_k
F = [0 1 3; 1 2 3; 0 3 2; 0 2 1]+1;
c = 0 ;

u = vpa(zeros(12,1));
v = vpa(zeros(12,1));
w = vpa(zeros(12,1));
for f = 1:4
    i = F(f,1);
    j = F(f,2);
    k = F(f,3);

    c = c + 1 ;
    u(c) = 16*lambda(j)*lambda(k);
    v(c) = -8*lambda(k)*lambda(i);
    w(c) = -8*lambda(i)*lambda(j);

    c = c + 1;
    u(c) = -8*lambda(j)*lambda(k);
    v(c) = 16*lambda(k)*lambda(i);
    w(c) = -8*lambda(i)*lambda(j);

    c = c + 1;
    u(c) = -8*lambda(j)*lambda(k);
    v(c) = -8*lambda(k)*lambda(i);
    w(c) = 16*lambda(i)*lambda(j);
end
