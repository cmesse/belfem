% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Prints one derivative table ( mWzeta ) from the generator in the C++
% syntax of cl_EF_TET10.cpp, one line per entry. The pattern for
% regenerating any other table: replace w and zeta.
clear all ;

tet10_generate ;

%for i = 1:12
%    disp( sprintf('mG( %i, k ) = %s ;', i-1, char(diff(w(i),zeta))) );
%end

for i = 1:12
    disp( sprintf('mWzeta( %i, k ) = %s ;', i-1, char(diff(w(i),zeta) )));
end