% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Cross-check of the g, h, u, v, w construction against the published
% DefElement expressions phi0 .. phi19 of the degree-2 Nedelec ( first kind )
% tetrahedron. The URL below is DefElement's triangle page; the reference
% expressions typed in here are the tetrahedron's. Several pairs match with
% a sign flip and are marked so. Every expand(...) must print zeros.
% Standalone; needs the Symbolic Math Toolbox.
clear all ;
clc ;

syms x y z real ;

% https://defelement.org/elements/examples/triangle-nedelec2-lagrange-2.html

X = [1 0 0 0 ];
Y = [0 0 1 0 ];
Z = [0 1 0 0 ];

N_xi = [1 0 0 -1; 0 0 1 -1; 0 1 0 -1];
J = N_xi*[X' Y' Z'];
iJ = inv(J);
nabla_xi = iJ(:,1);
nabla_eta = iJ(:,2);
nabla_zeta = iJ(:,3);
nabla_tau = -nabla_xi -nabla_eta-nabla_zeta ;

nabla = [nabla_xi, nabla_eta, nabla_zeta, nabla_tau ];

xi  = x ;
zeta = y ;
eta = z ;
tau = 1-xi-eta-zeta ;

g = vpa(zeros(12,1));
h = vpa(zeros(12,1));

lambda = vpa(zeros(4,1));
lambda(1) = xi ;
lambda(2) = zeta ;
lambda(3) = eta ;
lambda(4) = tau ;
%% edge 0 : xi->zeta
i = 1 ;
j = 2 ;
g(1) = 4*lambda(i)*(2*lambda(i)-1);
h(1) = 2*lambda(j)*(4*lambda(i)-1);


g(2) = 2*lambda(i)*(4*lambda(j)-1);
h(2) = 4*lambda(j)*(2*lambda(j)-1);

psi1 = g(1)*nabla(:,j) - h(1)*nabla(:,i) ;
psi2 = g(2)*nabla(:,j) - h(2)*nabla(:,i) ;

phi4 = [2*y*(1-4*x); 4*x*(2*x-1); 0 ];
phi5 = [4*y*(1-2*y); 2*x*(4*y-1); 0 ];

expand(phi4-psi1)
expand(phi5-psi2)

%% edge 1 : zeta->eta

phi0 = [0; 2*z*(1-4*y); 4*y*(2*y-1)];
phi1 = [0; 4*z*(1-2*z); 2*y*(4*z-1)];


i = 2 ;
j = 3 ;

g(3) = 4*lambda(i)*(2*lambda(i)-1);
h(3) = 2*lambda(j)*(4*lambda(i)-1);


g(4) = 2*lambda(i)*(4*lambda(j)-1);
h(4) = 4*lambda(j)*(2*lambda(j)-1);

psi3 = g(3)*nabla(:,j) - h(3)*nabla(:,i) ;
psi4 = g(4)*nabla(:,j) - h(4)*nabla(:,i) ;

expand(phi0-psi3)
expand(phi1-psi4)


%% edge 3 eta -> xi
i = 3 ;
j = 1 ;

g(5) = 4*lambda(i)*(2*lambda(i)-1);
h(5) = 2*lambda(j)*(4*lambda(i)-1);


g(6) = 2*lambda(i)*(4*lambda(j)-1);
h(6) = 4*lambda(j)*(2*lambda(j)-1);

psi5 = g(5)*nabla(:,j) - h(5)*nabla(:,i) ;
psi6 = g(6)*nabla(:,j) - h(6)*nabla(:,i) ;

phi2 = [ 2*z*(1-4*x); 0 ; 4*x*(2*x-1)];
phi3 = [ 4*z*(1-2*z); 0 ; 2*x*(4*z-1)];

% careful: these are swapped
expand(psi5+phi3)
expand(psi6+phi2)


%% edge 4 xi -> tau

i = 1 ;
j = 4 ;

g(7) = 4*lambda(i)*(2*lambda(i)-1);
h(7) = 2*lambda(j)*(4*lambda(i)-1);


g(8) = 2*lambda(i)*(4*lambda(j)-1);
h(8) = 4*lambda(j)*(2*lambda(j)-1);

psi7 = g(7)*nabla(:,j) - h(7)*nabla(:,i) ;
psi8 = g(8)*nabla(:,j) - h(8)*nabla(:,i) ;

phi10 = [8*x*y+8*x*z-6*x+8*y^2+16*y*z-12*y+8*z^2-12*z+4 ;
         2*x*(-4*x-4*y-4*z+3);
         2*x*(-4*x-4*y-4*z+3)];

phi11 = [-8*x*y-8*x*z+6*x+2*y+2*z-2; 4*x*(2*x-1); 4*x*(2*x-1)];

% careful: these are swapped
expand(psi7+phi11)
expand(psi8+phi10)

%% edge 5 : 2->4
i = 2 ;
j = 4 ;

g(9) = 4*lambda(i)*(2*lambda(i)-1);
h(9) = 2*lambda(j)*(4*lambda(i)-1);


g(10) = 2*lambda(i)*(4*lambda(j)-1);
h(10) = 4*lambda(j)*(2*lambda(j)-1);

psi9 = g(9)*nabla(:,j) - h(9)*nabla(:,i) ;
psi10 = g(10)*nabla(:,j) - h(10)*nabla(:,i) ;

phi8 = [2*y*(-4*x-4*y-4*z+3);
        8*x^2+8*x*y+16*x*z-12*x+8*y*z-6*y+8*z^2-12*z+4;
        2*y*(-4*x-4*y-4*z+3)];

phi9 = [4*y*(2*y-1); -8*x*y+2*x-8*y*z+6*y+2*z-2; 4*y*(2*y-1)];

% careful: these are swapped
expand(psi9+phi9)
expand(psi10+phi8)

%% edge 6 : 3->4
i = 3 ;
j = 4 ;

g(11) = 4*lambda(i)*(2*lambda(i)-1);
h(11) = 2*lambda(j)*(4*lambda(i)-1);


g(12) = 2*lambda(i)*(4*lambda(j)-1);
h(12) = 4*lambda(j)*(2*lambda(j)-1);

psi11 = g(11)*nabla(:,j) - h(11)*nabla(:,i) ;
psi12 = g(12)*nabla(:,j) - h(12)*nabla(:,i) ;

phi6 = [2*z*(-4*x-4*y-4*z+3);
    2*z*(-4*x-4*y-4*z+3); 
    8*x^2+16*x*y+8*x*z-12*x+8*y^2+8*y*z-12*y-6*z+4];
phi7 = [4*z*(2*z-1); 4*z*(2*z-1); -8*x*z+2*x-8*y*z+2*y+6*z-2];

% careful: these are swapped
expand(psi11+phi7)
expand(psi12+phi6)

%% surface 1 : 1->2->4
U = vpa(zeros(12,1));
V = vpa(zeros(12,1));
W = vpa(zeros(12,1));

i = 1 ;
j = 2 ;
k = 4 ;

U(1) = 16*lambda(j)*lambda(k);
V(1) = -8*lambda(i)*lambda(k);
W(1) = -8*lambda(i)*lambda(j);


U(2) = -8*lambda(j)*lambda(k);
V(2) = 16*lambda(i)*lambda(k);
W(2) = -8*lambda(i)*lambda(j);

%U(3) = -8*lambda(j)*lambda(k);
%V(3) = -8*lambda(i)*lambda(k);
%W(3) = 16*lambda(i)*lambda(j);

phiA = U(1)*nabla(:,i)+V(1)*nabla(:,j)+W(1)*nabla(:,k);
phiB = U(2)*nabla(:,i)+V(2)*nabla(:,j)+W(2)*nabla(:,k);
phiC = -phiA-phiB ;

phi18 = [8*y*(-x-2*y-2*z+2); 8*x*(x+2*y+z-1); 8*x*y];
phi19 = [8*y*(2*x+y+z-1); 8*x*(-2*x-y-2*z+2); 8*x*y];

expand(phiA-phi18)
expand(phiB-phi19)

%% surface 2 : 2->3->4
i = 2 ;
j = 3 ;
k = 4 ;

U(4) = 16*lambda(j)*lambda(k);
V(4) = -8*lambda(i)*lambda(k);
W(4) = -8*lambda(i)*lambda(j);


U(5) = -8*lambda(j)*lambda(k);
V(5) = 16*lambda(i)*lambda(k);
W(5) = -8*lambda(i)*lambda(j);

phiA = U(4)*nabla(:,i)+V(4)*nabla(:,j)+W(4)*nabla(:,k);
phiB = U(5)*nabla(:,i)+V(5)*nabla(:,j)+W(5)*nabla(:,k);
phiC = -phiA-phiB ;

phi14 = [8*y*z; 8*z*(-2*x-y-2*z+2); 8*y*(x+y+2*z-1)];
phi15 = [8*y*z; 8*z*(x+2*y+z-1); 8*y*(-2*x-2*y-z+2)];
expand(phiA-phi14)
expand(phiB-phi15)

%% surface 3 : 1->4->3
i = 1 ;
j = 4 ;
k = 3 ;

U(7) = 16*lambda(j)*lambda(k);
V(7) = -8*lambda(i)*lambda(k);
W(7) = -8*lambda(i)*lambda(j);


U(8) = -8*lambda(j)*lambda(k);
V(8) = 16*lambda(i)*lambda(k);
W(8) = -8*lambda(i)*lambda(j);

phiA = U(7)*nabla(:,i)+V(7)*nabla(:,j)+W(7)*nabla(:,k);
phiB = U(8)*nabla(:,i)+V(8)*nabla(:,j)+W(8)*nabla(:,k);
phiC = -phiA-phiB ;

phi16 = [8*z*(-x-2*y-2*z+2); 8*x*z; 8*x*(x+y+2*z-1)];
phi17 = [8*z*(2*x+y+z-1); 8*x*z; 8*x*(-2*x-2*y-z+2)];
expand(phi16-phiA)
expand(phi17-phiC)

%% surface 4 : 1->3->2

i = 1;
j = 3 ;
k = 2 ;

U(10) = 16*lambda(j)*lambda(k);
V(10) = -8*lambda(i)*lambda(k);
W(10) = -8*lambda(i)*lambda(j);
%phi15 = [8*y*z; 8*x*(x+2*y+z-1); 8*y*(-2*x-2*y-z+2)] ;


U(11) = -8*lambda(j)*lambda(k);
V(11) = 16*lambda(i)*lambda(k);
W(11) = -8*lambda(i)*lambda(j);

%U(3) = -8*lambda(j)*lambda(k);
%V(3) = -8*lambda(i)*lambda(k);
%W(3) = 16*lambda(i)*lambda(j);

phiA = U(10)*nabla(:,i)+V(10)*nabla(:,j)+W(10)*nabla(:,k);
phiB = U(11)*nabla(:,i)+V(11)*nabla(:,j)+W(11)*nabla(:,k);
phiC = -phiA-phiB ;

phi12 = [-8*y*z; 16*x*z; -8*x*y ];
phi13 = [-8*y*z; -8*x*z; 16*x*y ];

expand(phi12+phiA+phiB)
expand(phi13-phiB)