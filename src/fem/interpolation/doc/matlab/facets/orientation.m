% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Interactive notebook, not a generator: the PENTA6 master faces in EXODUS
% order, every slave ordering of each face, and for each ( face, orientation )
% the parameter map eta = f( xi ) from the master face onto the slave
% element. These maps are what the slave intpoints_penta overload in
% fn_IF_initialize_integration_points_on_facet.cpp encodes. Pick a slave
% ordering with S = ..., keep only the matching eta block ( each later
% block overwrites the earlier one ), and the plot shows the master point p
% and the slave point q coinciding when the map is right.
clear all ;
clf ;
clc ;

a = (0.25 * sqrt(3))^(-1/3) ;
b = 0.5 * sqrt(3) * a;
h = a;

X = [ ...
   -0.5*a, -b/3, -0.5*h;
    0.5*a, -b/3, -0.5*h;
    0.0  ,  2*b/3, -0.5*h;
   -0.5*a, -b/3,  0.5*h;
    0.5*a, -b/3,  0.5*h;
    0.0  ,  2*b/3,  0.5*h ];

edges = [ ...
    1 2;
    2 3;
    3 1;
    4 5;
    5 6;
    6 4;
    1 4;
    2 5;
    3 6 ];

V = 0.5 * a * b * h ;



% faces on the element in EXODUS order
F1 = [ 1, 2, 5, 4 ];
F2 = [ 2, 3, 6, 5 ];
F3 = [ 1, 4, 6, 3 ];
F4 = [ 1, 3, 2 ];
F5 = [ 4, 5, 6 ];

% slave faces, pointing inward toward the element
S1a = [ 1, 4, 5, 2 ] ; % face 1, orientation 1
S1b = [ 2, 1, 4, 5 ] ; % face 1, orientation 2
S1c = [ 5, 2, 1, 4 ] ; % face 1, orientation 3
S1d = [ 4, 5, 2, 1 ] ; % face 1, orientation 4

S2a = [ 2, 5, 6, 3 ] ; % face 2, orientation 1
S2b = [ 3, 2, 5, 6 ] ; % face 2, orientation 2
S2c = [ 6, 3, 2, 5 ] ; % face 2, orientation 3
S2d = [ 5, 6, 3, 2 ] ; % face 2, orientation 4

S3a = [ 1, 3, 6, 4 ] ; % face 3, orientation 1
S3b = [ 4, 6, 3, 1 ] ; % face 3, orientation 2
S3c = [ 6, 3, 1, 4 ] ; % face 3, orientation 3
S3d = [ 3, 6, 4, 1 ] ; % face 3, orientation 4

S4a = [ 1, 2 , 3 ]; % face 4, orientation 1
S4b = [ 2, 3 , 1 ]; % face 4, orientation 2
S4c = [ 3, 1 , 2 ]; % face 4, orientation 3

S5a = [ 4, 6 , 5 ]; % face 5, orientation 1
S5b = [ 5, 4 , 6 ]; % face 5, orientation 2
S5c = [ 6, 5 , 4 ]; % face 5, orientation 3


% select the slave face here
S = S5c ;

% plot
clf;
plot3(X(:,1), X(:,2), X(:,3), 'ko', 'MarkerSize', 8, 'LineWidth', 1);
hold on;

for k = 1:size(edges,1)
    plot3(X(edges(k,:),1), X(edges(k,:),2), X(edges(k,:),3), '--k', 'LineWidth', 0.5);
end

N6 = @(xi,eta,zeta)[ -(xi*(zeta-1.0)), -(eta*(zeta-1.0)), ((zeta-1.0)*(eta+xi-1.0)), (xi*(zeta + 1.0)), (eta*(zeta + 1.0)), ((zeta+1.0)*(1.0-xi-eta))]*0.5 ;

N4 = @(xi,eta)[ (1-xi)*(1-eta), (1+xi)*(1-eta), (1+xi)*(1+eta), (1-xi)*(1+eta) ]/4 ;
N3 = @(xi,eta)[ xi, eta, (1-xi-eta) ] ;

x = X([ S, S(1) ], 1) ;
y = X([ S, S(1) ], 2) ;
z = X([ S, S(1) ], 3) ;

% reference point on the master surface
xi = [ 0.75, 0.05, 0 ];

eta = zeros( 3, 1 );

% face 1 orientation 1
eta(1) = 1 - 0.5 * (xi(2) + 1) ;
eta(2) = 0.5 * (xi(2) + 1) ;
eta(3) = xi(1);

% face 1 orientation 2
eta(1) = 0.5 * (xi(1) + 1) ;
eta(2) = 1-0.5 * (xi(1) + 1) ;
eta(3) = xi(2);

% face 1 orientation 3
eta(1) = 0.5 * (xi(2) + 1) ;
eta(2) = 1-0.5 * (xi(2) + 1) ;
eta(3) = -xi(1);

% face 1 orientation 4
eta(1) = 1-0.5 * (xi(1) + 1) ;
eta(2) = 0.5 * (xi(1) + 1) ;
eta(3) = -xi(2);

% face 2 orientation 1
eta(1) = 0 ;
eta(2) = 1-0.5 * (xi(2) + 1) ;
eta(3) = xi(1);

% face 2 orientation 2
eta(1) = 0 ;
eta(2) = 0.5 * (xi(1) + 1) ;
eta(3) = xi(2);

% face 2 orientation 3
eta(1) = 0 ;
eta(2) = 0.5 * (xi(2) + 1) ;
eta(3) = -xi(1);

% face 2 orientation 4
eta(1) = 0 ;
eta(2) = 1-0.5 * (xi(1) + 1) ;
eta(3) = -xi(2);

% face 3 orientation 1
eta(1) = 1-0.5 * (xi(1) + 1) ;
eta(2) = 0 ;
eta(3) = xi(2);

% face 3 orientation 2
eta(1) = 1 - 0.5 * (xi(1) + 1) ;
eta(2) = 0 ;
eta(3) = -xi(2);

%% face 3 orientation 3
eta(1) = 0.5 * (xi(2) + 1) ;
eta(2) = 0 ;
eta(3) = -xi(1);

%% face 3 orientation 4
eta(1) = 0.5 * (xi(2) + 1) ;
eta(2) = 0 ;
eta(3) = xi(1);
1
%% face 4 orientation 1
eta(1) = xi(1) ;
eta(2) = xi(2) ;
eta(3) = -1

%% face 4 orientation 2
eta(1) = 1-xi(1)-xi(2) ;
eta(2) = xi(1) ;
eta(3) = -1

%% face 4 orientation 3
eta(1) = xi(2) ;
eta(2) = 1-xi(1)-xi(2) ;
eta(3) = -1

%% face 5 orientation 1
eta(1) = xi(1) ;
eta(2) = 1-xi(1)-xi(2) ;
eta(3) = 1

%% face 5 orientation 2
eta(1) = xi(2) ;
eta(2) = xi(1) ;
eta(3) = 1

%% face 5 orientation 3
eta(1) = 1-xi(1)-xi(2) ;
eta(2) = xi(2);
eta(3) = 1

% map to the triangle

%p = N4(xi(1),xi(2)) * X(S,:) ;
p = N3(xi(1),xi(2)) * X(S,:) ;

q = N6(eta(1),eta(2),eta(3)) * X ;

fill3(x, y, z, 'y', 'FaceAlpha', 0.5);
plot3(p(1), p(2), p(3), 'rx', 'MarkerSize', 10, 'LineWidth', 2);
plot3(q(1), q(2), q(3), 'bo', 'MarkerSize', 10, 'LineWidth', 2);
axis equal;
axis off ;