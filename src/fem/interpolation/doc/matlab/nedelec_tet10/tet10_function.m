% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Zero-test of the mG, mH, mU, mV, mW tables of cl_EF_TET10.cpp against the
% generator. The tables below are transcriptions of the C++ under the pinned
% node map; compare_tables.py checks the transcription against the source
% file. Every expand(...) must print zeros; check_tet10.m asserts it. Runs
% tet10_generate.m and tet10_shortcuts.m.
clc ;
tet10_generate ;

%%
tet10_shortcuts ;

% eta and zeta are exchanged relative to the naive node map; see the
% pin in tet10_generate.m
mG = [xi4*(xi2-1)
xi2*(zeta4-1)
zeta4*(zeta2-1)
zeta2*(eta4-1)
eta4*(eta2-1)
eta2*(xi4-1)
xi4*(xi2-1)
xi2*(3-eta4-xi4-zeta4)
zeta4*(zeta2-1)
zeta2*(3-eta4-xi4-zeta4)
eta4*(eta2-1)
eta2*(3-eta4-xi4-zeta4)];

mH = [zeta2*(xi4-1)
zeta4*(zeta2-1)
eta2*(zeta4-1)
eta4*(eta2-1)
xi2*(eta4-1)
xi4*(xi2-1)
(xi4-1)*tau2
(1-eta2-xi2-zeta2)*tau4
(zeta4-1)*tau2
(1-eta2-xi2-zeta2)*tau4
(eta4-1)*tau2
(1-eta2-xi2-zeta2)*tau4];

mU = [16*zeta*tau
-zeta8*tau
16*eta*tau
-eta8*tau
16*eta*tau
-eta8*tau
16*eta*zeta
-eta8*zeta];

mV = [-xi8*tau
16*xi*tau
-zeta8*tau
16*zeta*tau
-xi8*eta
16*xi*eta
-xi8*zeta
16*xi*zeta];

mW = [-zeta8*xi
-zeta8*xi
-zeta8*eta
-zeta8*eta
-xi8*tau
-xi8*tau
-xi8*eta
-xi8*eta];

% faces: the generator emits 3 candidates per face; only the first two are active
act = [1 2 4 5 7 8 10 11];
expand(g-mG)'
expand(h-mH)'
expand(u(act)-mU)'
expand(v(act)-mV)'
expand(w(act)-mW)'