% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Zero-test of cl_IF_PENTA18.hpp d2NdXi2. N is the C++ shape-function set,
% differentiated twice symbolically; ad2NdXi2 below is a transcription of
% the C++ table ( compare_tables.py checks the transcription against the
% header ). X must be all zeros; check_penta18.m asserts it. Needs the
% Symbolic Math Toolbox.
clear all ;
clc ;

syms xi eta zeta real ;
zeta2 = zeta*zeta ;


N = vpa(zeros(18,1));

N(  1 ) = 0.5*(xi*zeta*(xi*2.0-1.0)*(zeta-1.0));
N(  2 ) = 0.5*(eta*zeta*(eta*2.0-1.0)*(zeta-1.0));
N(  3 ) = 0.5*((1.0-xi-eta)*(2*(xi+eta)-1.0)*(1.0-zeta)*zeta);
N(  4 ) = 0.5*(xi*zeta*(xi*2.0-1.0)*(zeta+1.0));
N(  5 ) = 0.5*(eta*zeta*(eta*2.0-1.0)*(zeta+1.0));
N(  6 ) = 0.5*((1.0-xi-eta)*(1.0-2.0*(xi+eta))*zeta*(1.0+zeta));
N(  7 ) = xi*eta*zeta*(zeta-1.0)*2.0;
N(  8 ) = 2.0*eta*zeta*(1.0-zeta)*(eta+xi-1.0);
N(  9 ) = 2.0*xi*zeta*(1.0-zeta)*(eta+xi-1.0);
N( 10 ) = xi*(1.0-2.0*xi)*(zeta2-1.0);
N( 11 ) = eta*(1.0-2.0*eta)*(zeta2-1.0);
N( 12 ) = (1-xi-eta)*(2.0*(xi+eta)-1)*(zeta2-1.0);
N( 13 ) = xi*eta*zeta*(1.0+zeta)*2.0;
N( 14 ) = 2.0*eta*zeta*(1.0+zeta)*(1.0-xi-eta);
N( 15 ) = 2.0*xi*zeta*(1.0+zeta)*(1.0-xi-eta);
N( 16 ) = 4.0*eta*xi*(1.0-zeta2);
N( 17 ) = 4.0*eta*(1.0-zeta2)*(1.0-xi-eta);
N( 18 ) = 4.0*xi*(1.0-zeta2)*(1.0-xi-eta);

dNdxi = vpa(zeros(18,3));

dNdxi(:,1) = diff(N,xi);
dNdxi(:,2) = diff(N,eta);
dNdxi(:,3) = diff(N,zeta);

d2Ndxi2 = vpa(zeros(18,6));
d2Ndxi2(:,1) = diff(N,xi,xi);
d2Ndxi2(:,2) = diff(N,eta,eta);
d2Ndxi2(:,3) = diff(N,zeta,zeta);

d2Ndxi2(:,4) = diff(N,eta,zeta);
d2Ndxi2(:,5) = diff(N,zeta,xi);
d2Ndxi2(:,6) = diff(N,xi,eta);
 
ad2NdXi2 = vpa(zeros(6,18));

ad2NdXi2( 1,  1 ) = 2.0*(zeta-1.0)*zeta;
ad2NdXi2( 2,  1 ) = 0.0;
ad2NdXi2( 3,  1 ) = xi*(2.0*xi-1.0);
ad2NdXi2( 4,  1 ) = 0.0;
ad2NdXi2( 5,  1 ) = 0.5-zeta+2.0*xi*(2.0*zeta-1.0);
ad2NdXi2( 6,  1 ) = 0.0;

ad2NdXi2( 1,  2 ) = 0.0;
ad2NdXi2( 2,  2 ) = 2.0*(zeta-1.0)*zeta;
ad2NdXi2( 3,  2 ) = eta*(2.0*eta-1.0);
ad2NdXi2( 4,  2 ) = 0.5-zeta+2.0*eta*(2.0*zeta-1.0);
ad2NdXi2( 5,  2 ) = 0.0;
ad2NdXi2( 6,  2 ) = 0.0;

ad2NdXi2( 1,  3 ) = 2.0*zeta*(zeta-1.0);
ad2NdXi2( 2,  3 ) = 2.0*zeta*(zeta-1.0);
ad2NdXi2( 3,  3 ) = 2.0*(1.0-xi-eta)*(0.5-eta-xi);
ad2NdXi2( 4,  3 ) = 1.5-3.0*zeta+2.0*(xi+eta)*(2.0*zeta-1.0);
ad2NdXi2( 5,  3 ) = 1.5-3.0*zeta+2.0*(xi+eta)*(2.0*zeta-1.0);
ad2NdXi2( 6,  3 ) = 2.0*zeta*(zeta-1.0);

ad2NdXi2( 1,  4 ) = 2.0*zeta*(1.0+zeta);
ad2NdXi2( 2,  4 ) = 0.0;
ad2NdXi2( 3,  4 ) = xi*(2.0*xi-1.0);
ad2NdXi2( 4,  4 ) = 0.0;
ad2NdXi2( 5,  4 ) = 2.0*xi*(1.0+2.0*zeta)-zeta-0.5;
ad2NdXi2( 6,  4 ) = 0.0;

ad2NdXi2( 1,  5 ) = 0.0;
ad2NdXi2( 2,  5 ) = 2.0*zeta*(1.0+zeta);
ad2NdXi2( 3,  5 ) = eta*(2.0*eta-1.0);
ad2NdXi2( 4,  5 ) = 2.0*eta*(1.0+2.0*zeta)-0.5-zeta;
ad2NdXi2( 5,  5 ) = 0.0;
ad2NdXi2( 6,  5 ) = 0.0;

ad2NdXi2( 1,  6 ) = 2.0*zeta*(1.0+zeta);
ad2NdXi2( 2,  6 ) = 2.0*zeta*(1.0+zeta);
ad2NdXi2( 3,  6 ) = 2.0*(1.0-xi-eta)*(0.5-xi-eta);
ad2NdXi2( 4,  6 ) = 2.0*(xi+eta)*(1.0+2.0*zeta)-3.0*zeta-1.5;
ad2NdXi2( 5,  6 ) = 2.0*(xi+eta)*(1.0+2.0*zeta)-3.0*zeta-1.5;
ad2NdXi2( 6,  6 ) = 2.0*zeta*(1.0+zeta);

ad2NdXi2( 1,  7 ) = 0.0;
ad2NdXi2( 2,  7 ) = 0.0;
ad2NdXi2( 3,  7 ) = 4.0*xi*eta;
ad2NdXi2( 4,  7 ) = 2.0*xi*(2.0*zeta-1.0);
ad2NdXi2( 5,  7 ) = 2.0*eta*(2.0*zeta-1.0);
ad2NdXi2( 6,  7 ) = 2.0*zeta*(zeta-1.0);

ad2NdXi2( 1,  8 ) = 0.0;
ad2NdXi2( 2,  8 ) = 4.0*(1.0-zeta)*zeta;
ad2NdXi2( 3,  8 ) = 4.0*eta*(1.0-xi-eta);
ad2NdXi2( 4,  8 ) = 2.0*(xi+2.0*eta)*(1.0-2.0*zeta)+4.0*zeta-2.0;
ad2NdXi2( 5,  8 ) = eta*(2.0-4.0*zeta);
ad2NdXi2( 6,  8 ) = 2.0*(1.0-zeta)*zeta;

ad2NdXi2( 1,  9 ) = 4.0*(1.0-zeta)*zeta;
ad2NdXi2( 2,  9 ) = 0.0;
ad2NdXi2( 3,  9 ) = 4.0*xi*(1.0-xi-eta);
ad2NdXi2( 4,  9 ) = 2.0*xi*(1.0-2.0*zeta);
ad2NdXi2( 5,  9 ) = 4.0*zeta+2.0*(1.0-2.0*zeta)*(2.0*xi+eta)-2.0;
ad2NdXi2( 6,  9 ) = 2.0*(1.0-zeta)*zeta;

ad2NdXi2( 1, 10 ) = 4.0*(1.0-zeta2);
ad2NdXi2( 2, 10 ) = 0.0;
ad2NdXi2( 3, 10 ) = (2.0-4.0*xi)*xi;
ad2NdXi2( 4, 10 ) = 0.0;
ad2NdXi2( 5, 10 ) = (2.0-8.0*xi)*zeta;
ad2NdXi2( 6, 10 ) = 0.0;

ad2NdXi2( 1, 11 ) = 0.0;
ad2NdXi2( 2, 11 ) = 4.0*(1.0-zeta2);
ad2NdXi2( 3, 11 ) = (2.0-4.0*eta)*eta;
ad2NdXi2( 4, 11 ) = (2.0-8.0*eta)*zeta;
ad2NdXi2( 5, 11 ) = 0.0;
ad2NdXi2( 6, 11 ) = 0.0;

ad2NdXi2( 1, 12 ) = 4.0*(1.0-zeta2);
ad2NdXi2( 2, 12 ) = 4.0*(1.0-zeta2);
ad2NdXi2( 3, 12 ) = 4.0*(1.0-xi-eta)*(xi+eta-0.5);
ad2NdXi2( 4, 12 ) = (6.0-8.0*(xi+eta))*zeta;
ad2NdXi2( 5, 12 ) = (6.0-8.0*(xi+eta))*zeta;
ad2NdXi2( 6, 12 ) = 4.0*(1.0-zeta2);

ad2NdXi2( 1, 13 ) = 0.0;
ad2NdXi2( 2, 13 ) = 0.0;
ad2NdXi2( 3, 13 ) = 4.0*xi*eta;
ad2NdXi2( 4, 13 ) = xi*(2.0+4.0*zeta);
ad2NdXi2( 5, 13 ) = eta*(2.0+4.0*zeta);
ad2NdXi2( 6, 13 ) = 2.0*zeta*(1.0+zeta);

ad2NdXi2( 1, 14 ) = 0.0;
ad2NdXi2( 2, 14 ) = -4.0*zeta*(1.0+zeta);
ad2NdXi2( 3, 14 ) = 4.0*eta*(1.0-xi-eta);
ad2NdXi2( 4, 14 ) = 2.0*(1.0+2.0*zeta)*(1.0-xi-2.0*eta);
ad2NdXi2( 5, 14 ) = -2.0*eta*(1.0+2.0*zeta);
ad2NdXi2( 6, 14 ) = -2.0*zeta*(1.0+zeta);

ad2NdXi2( 1, 15 ) = -4.0*zeta*(1.0+zeta);
ad2NdXi2( 2, 15 ) = 0.0;
ad2NdXi2( 3, 15 ) = 4.0*xi*(1.0-xi-eta);
ad2NdXi2( 4, 15 ) = -xi*(2.0+4.0*zeta);
ad2NdXi2( 5, 15 ) = 2.0*(1.0+2.0*zeta)*(1.0-2.0*xi-eta);
ad2NdXi2( 6, 15 ) = -2.0*zeta*(1.0+zeta);

ad2NdXi2( 1, 16 ) = 0.0;
ad2NdXi2( 2, 16 ) = 0.0;
ad2NdXi2( 3, 16 ) = -8.0*xi*eta;
ad2NdXi2( 4, 16 ) = -8.0*xi*zeta;
ad2NdXi2( 5, 16 ) = -8.0*eta*zeta;
ad2NdXi2( 6, 16 ) = 4.0-4.0*zeta2;

ad2NdXi2( 1, 17 ) = 0.0;
ad2NdXi2( 2, 17 ) = 8.0*zeta2-8.0;
ad2NdXi2( 3, 17 ) = 8.0*eta*(xi+eta-1.0);
ad2NdXi2( 4, 17 ) = (8.0*xi+16.0*eta-8.0)*zeta;
ad2NdXi2( 5, 17 ) = 8.0*eta*zeta;
ad2NdXi2( 6, 17 ) = 4.0*zeta2-4.0;

ad2NdXi2( 1, 18 ) = 8.0*zeta2-8.0;
ad2NdXi2( 2, 18 ) = 0.0;
ad2NdXi2( 3, 18 ) = 8.0*xi*(xi+eta-1.0);
ad2NdXi2( 4, 18 ) = 8.0*xi*zeta;
ad2NdXi2( 5, 18 ) = (16.0*xi+8.0*eta-8.0)*zeta;
ad2NdXi2( 6, 18 ) = 4.0*zeta2-4.0;

X = expand(ad2NdXi2' - d2Ndxi2);
