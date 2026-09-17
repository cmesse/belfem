% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Derivation note: the TET10 Lagrange shape functions of cl_IF_TET10.hpp by
% inverting a Vandermonde matrix over the Pascal basis, with the nodes in
% EXODUS order ( xi_hat, eta_hat, zeta_hat below: node 1 carries zeta, node
% 2 carries eta ). N is the C++ result. The plot afterward is a
% two-tet sketch and not part of the derivation.
clear all;
clc;

a = 1;
h = sqrt(3)*0.5*a;
t = a * sqrt(6)/3 ;
n = 4;
GRID = zeros(11, 3);

GRID(1,1) = -0.5*a;
GRID(1,2) = -2/3*h;

GRID(2,1) =  0.5*a;
GRID(2,2) = -2/3*h;

GRID(3,1) =      0;
GRID(3,2) =    h/3;

GRID(4,2) =  -1/3*h;
GRID(4,3) = t;

GRID(5,1) = 0;
GRID(5,2) = -2/3*h;
GRID(5,3) = 0 ;

GRID(6,1) = a/4 ;
GRID(6,2) = -h/6;
GRID(6,3) = 0 ;

GRID(7,1) = -a/4 ;
GRID(7,2) = -h/6;
GRID(7,3) = 0 ;

GRID(8,1) = -a/4 ;
GRID(8,2) = -h/2;
GRID(8,3) =  t/2 ;

GRID(9,1) = a/4 ;
GRID(9,2) = -h/2;
GRID(9,3) = t/2 ;

GRID(10,1) = 0 ;
GRID(10,2) = 0;
GRID(10,3) = t/2 ;

GRID(11,2) =  -1/3*h;
GRID(11,3) = -t;

GRID(12,1) = -0.25*a ;
GRID(12,2) = -h/2 ;
GRID(12,3) = -t/2 ;

GRID(13,1) = 0 ;
GRID(13,2) = 0 ;
GRID(13,3) = -t/2 ;

GRID(14,1) = 0.25*a ;
GRID(14,2) = -h/2 ;
GRID(14,3) = -t/2 ;

% midside parameter
b = 0.5;

% reference node parameters in EXODUS order ( row vectors )
xi_hat =    [1 0 0 0 b 0 b b 0 0];
eta_hat =   [0 0 1 0 0 b b 0 0 b];
zeta_hat =  [0 1 0 0 b b 0 0 b 0];

% polynomial basis from the Pascal tetrahedron ( row vector )
p = @(xi,eta, zeta)[1 xi eta zeta xi^2 eta^2 zeta^2 eta*zeta xi*zeta xi*eta];

%% empty if nothing is to be plotted
plotpattern = [];

%% end of user input
number_of_nodes = max(size(xi_hat));

%% build the coefficient matrix
V = zeros(number_of_nodes);
for i = 1:number_of_nodes
    V(i,:) = p(xi_hat(i),eta_hat(i),zeta_hat(i));
end


%% compute the shape functions ( row vector )
syms xi eta zeta real;


N = vpa(zeros(1, number_of_nodes));
for i = 1:number_of_nodes
    f = zeros(number_of_nodes, 1);
    f(i) = 1;
    coeff = V\f;
    N(i) = simplify(p(xi, eta, zeta)*coeff);
end

phi = @(xi,eta,zeta)[xi*(2*xi - 1), zeta*(2*zeta - 1), eta*(2*eta - 1), 2*eta^2 + 4*eta*xi + 4*eta*zeta - 3*eta + 2*xi^2 + 4*xi*zeta - 3*xi + 2*zeta^2 - 3*zeta + 1, 4*xi*zeta, 4*eta*zeta, 4*eta*xi, -4*xi*(eta + xi + zeta - 1), -4*zeta*(eta + xi + zeta - 1), -4*eta*(eta + xi + zeta - 1)];

%%
clf ;
plot3(GRID(:,1), GRID(:,2), GRID(:,3), 'ko', 'MarkerSize', 8  );
axis off ;
axis equal ;
hold on ;
eidx = [ 1 2; 2 3; 3 1 ; 1 4; 2 4; 3 4 ];
for e = 1:6
    plot3(GRID(eidx(e,:),1), GRID(eidx(e,:),2), GRID(eidx(e,:),3), '--k' ) ;
end

%%
N_tet4 = @(xi,eta,zeta)[ xi zeta eta 1-xi-eta-zeta];

idx = [1 3 2 11 ]
el2 = [1 3 2 11 7 6 5 12 13 14 ]

%% create the edges
Edges = [ 1 3 ; 3 2 ; 2 1 ; 1 11 ; 3 11 ; 2 11 ];
for e = 1:6
    plot3(GRID(Edges(e,:),1), GRID(Edges(e,:),2), GRID(Edges(e,:),3), '-k' ) ;
end

k = 6

el2(eidx(k,:))
%p = N_tet4(xi_hat(k), eta_hat(k), zeta_hat(k))*GRID(idx,:);
%plot3(p(1), p(2), p(3), 'rx', 'MarkerSize', 8 );
%plot3(GRID(k,1), GRID(k,2), GRID(k,3),'bo');


