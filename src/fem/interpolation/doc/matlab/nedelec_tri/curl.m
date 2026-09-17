% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Derivation note, not a generator: the curl of the three signed TRI3
% Whitney functions in parameter space, transformed with the inverse
% Jacobian. The Whitney sign convention here ( eta*xi_x - xi*eta_x ) is the
% opposite of cl_EF_TRI3.cpp, which uses G*nabla_eta - H*nabla_xi. The
% collapse to the constant operator C = 2/detJ [ s1 s2 s3 ] of
% nedelec_derivation.md, section 3.2, is left as the commented line at the
% end. Needs the Symbolic Math Toolbox.
clear all;

syms xi eta real ;
syms xi_x xi_y eta_x eta_y real ;
syms s1 s2 s3 real ;

iJ = [xi_x xi_y; eta_x eta_y]';

zeta = 1 - xi - eta ;
zeta_x = -xi_x -eta_x ;
zeta_y = -xi_y - eta_y ;

phi_1x = s1 * ( eta*xi_x - xi * eta_x );
phi_2x = s2 * ( zeta*eta_x - eta * zeta_x );
phi_3x = s3 * ( xi*zeta_x - zeta * xi_x );


phi_1y = s1 * ( eta*xi_y - xi * eta_y );
phi_2y = s2 * ( zeta*eta_y - eta * zeta_y );
phi_3y = s3 * ( xi*zeta_y - zeta * xi_y );

phi1_x_xi = diff( phi_1x, xi );
phi1_x_eta = diff( phi_1x, eta );

phi2_x_xi = diff( phi_2x, xi );
phi2_x_eta = diff( phi_2x, eta );

phi3_x_xi = diff( phi_3x, xi );
phi3_x_eta = diff( phi_3x, eta );

phi1_y_xi = diff( phi_1y, xi );
phi1_y_eta = diff( phi_1y, eta );

phi2_y_xi = diff( phi_2y, xi );
phi2_y_eta = diff( phi_2y, eta );

phi3_y_xi = diff( phi_3y, xi );
phi3_y_eta = diff( phi_3y, eta );

a = [phi1_x_xi phi2_x_xi phi3_x_xi;
     phi1_x_eta phi2_x_eta phi3_x_eta];
 
A = iJ*a

b = [phi1_y_xi phi2_y_xi phi3_y_xi;
     phi1_y_eta phi2_y_eta phi3_y_eta];

B = iJ*b
 



%simplify([ (B(2,:);-A(1,:) ])/(eta_x*xi_y-eta_y*xi_x))
 