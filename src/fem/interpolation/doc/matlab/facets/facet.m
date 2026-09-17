% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Interactive notebook, not a generator: for PENTA6 face f, it plots the
% outward normal at the face center, built from two rows of dN/dxi * X. The
% if-chain at the end records which rows, with which signs, give the
% outward normal on each of the five faces. Set f and run; the arrow must
% point out of the element.
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

Xi = [ 1, 0, -1 ; 
 0, 1, -1 ; 
 0, 0, -1 ; 
 1, 0, 1 ; 
 0, 1, 1 ; 
 0, 0, 1 ] ;

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
% faces on the element in EXODUS order
F1 = [ 1, 2, 5, 4 ];
F2 = [ 2, 3, 6, 5 ];
F3 = [ 1, 4, 6, 3 ];
F4 = [ 1, 3, 2 ];
F5 = [ 4, 5, 6 ];

N = @(xi,eta,zeta)[ -(xi*(zeta-1.0)), -(eta*(zeta-1.0)), ((zeta-1.0)*(eta+xi-1.0)), (xi*(zeta + 1.0)), (eta*(zeta + 1.0)), ((zeta+1.0)*(1.0-xi-eta))]*0.5 ;

N_xi = @(xi,eta,zeta)[1 / 2 - zeta / 2, 0, zeta / 2 - 1 / 2, zeta / 2 + 1 / 2, 0, -zeta / 2 - 1 / 2 ;
 0, 1 / 2 - zeta / 2, zeta / 2 - 1 / 2, 0, zeta / 2 + 1 / 2, -zeta / 2 - 1 / 2;
-xi / 2, -eta / 2, eta / 2 + xi / 2 - 1 / 2, xi / 2, eta / 2, -eta / 2 - xi / 2 + 1 / 2];

% identify points
Xi = zeros( 3, 6 );



% plot
clf;Xi = [ 1, 0, -1 ; 0, 1, -1 ; 0, 0, -1 ; 1, 0, 1 ; 0, 1, 1 ; 0, 0, 1 ] ;
%plot3(X(:,1), X(:,2), X(:,3), 'ko', 'MarkerSize', 8, 'LineWidth', 1);
hold on;
axis off ;
axis equal ;

for k = 1:size(edges,1)
 plot3(X(edges(k,:),1), X(edges(k,:),2), X(edges(k,:),3), '--k', 'LineWidth', 0.5);
end

f = 5 ;
if( f==1 )
 F = F1 ;
elseif( f==2 )
 F = F2 ;
elseif( f==3 )
 F = F3 ;
elseif( f==4 )
 F = F4 ;
elseif( f==5 )
 F = F5 ;
end

nk = length(F) ;

fill3(X(F,1), X(F,2), X(F,3), 'y', 'FaceAlpha', 0.5);

xi = [ sum(Xi(F,1))/nk , sum(Xi(F,2))/nk, sum(Xi(F,3))/nk ] ;

p = N(xi(1), xi(2), xi(3)) * X ;

plot3(p(1), p(2), p(3), 'rx', 'MarkerSize', 8, 'LineWidth', 1);

syms J00 J01 J02 J10 J11 J12 J20 J21 J22  real ;

dN = N_xi(xi(1), xi(2), xi(3));

J = dN * X ;

J = [ J00, J01, J02 ;
J10, J11, J12 ;
J20, J21, J22 ] ;

if( f==1 )
 u = J(2,:)-J(1,:) ;
 v = J(3,:) ;
elseif( f==2 )
 u = -J(2,:) ;
 v = J(3,:) ;
elseif( f==3 )
 u = J(1,:) ;
 v = J(3,:) ;
elseif( f==4 )
 u = J(1,:) ;
 v = -J(2,:) ;
elseif( f==5 )
 u = J(1,:) ;
 v = J(2,:) ;
end

n = cross(u,v);

quiver3(p(1), p(2), p(3),n(1), n(2), n(3), 'r', 'LineWidth', 1);

