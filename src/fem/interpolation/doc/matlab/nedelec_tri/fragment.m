% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Derivation note: the TRI6 Nedelec ansatz behind the mG, mH tables of
% cl_EF_TRI6.cpp. Edge function k is theta = a*nabla_j + b*nabla_i with
% a = 4*xi_i*( 2*xi_i - 1 ) and b = 2*xi_j*( 1 - 4*xi_i ); k = 7, 8 are the
% two face functions. The C++ carries the expanded polynomials at half these
% coefficients ( first pair checked: mG( 0, k ) = eta*( 1 - xi4 ) is b/2 ),
% and its normalization is the unit circulation that
% tests/fem/test_EdgeFunctions.cpp checks. Originally a loop body; the
% preamble and the loop around it make it runnable.
syms xi eta real ;
syms nabla_xi_1 nabla_xi_2 nabla_eta_1 nabla_eta_2 real ;
nabla_xi  = [ nabla_xi_1;  nabla_xi_2  ];
nabla_eta = [ nabla_eta_1; nabla_eta_2 ];
theta_all = sym( zeros( 2, 8 ) );
for k = 1:8
		if ( k == 1 || k==2 )
                xi_i = xi ;
                xi_j = eta ;
                nabla_i = nabla_xi ;
                nabla_j = nabla_eta ; 
        elseif ( k == 3 || k==4 )
                xi_i = eta ;
                xi_j = 1-xi-eta;
                nabla_i = nabla_eta ;
                nabla_j = -nabla_xi-nabla_eta ;  
        elseif ( k == 5 || k==6 )
                xi_i = 1-xi-eta;
                xi_j = xi ;
                nabla_i = -nabla_xi-nabla_eta ;
                nabla_j = nabla_xi ;  
        end
        
        if ( k == 1 || k==3 || k==5 )
            a = 4*xi_i*(2*xi_i - 1);
            b = 2*xi_j*(1-4*xi_i);
            theta = a*nabla_j + b*nabla_i ;
        elseif ( k==2 || k==4 || k==6 )
            c = 2*(4*xi_j-1)*xi_i;
            d = 4*xi_j*(1-2*xi_j);
            theta = c*nabla_j + d*nabla_i ;
        elseif k==7
            e =  8*eta*(eta + 2*xi - 1);
            f =  -8*xi*(eta + 2*xi - 2);
            theta = e*nabla_xi + f*nabla_eta;
         elseif k==8
            g =  -8*eta*(xi - eta + 1);
            h =  -8*xi*(eta - xi + 1);
             theta = g*nabla_xi + h*nabla_eta;
        end
        theta_all( :, k ) = simplify( theta );
end
theta_all
