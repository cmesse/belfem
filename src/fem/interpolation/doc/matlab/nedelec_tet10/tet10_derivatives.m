% BELFEM -- The Berkeley Lab Finite Element Framework
% Copyright (c) 2026, The Regents of the University of California,
% through Lawrence Berkeley National Laboratory (subject to receipt of any required
% approvals from the U.S. Dept. of Energy).  All rights reserved.
%
% Developers: Christian Messe, Gregory Giard
%
% See the top-level LICENSE file for the complete license and disclaimer.

% Zero-test of the fifteen derivative tables of cl_EF_TET10.cpp ( mGxi ...
% mWzeta ) against the derivatives of the generator's g, h, u, v, w. The
% tables below are transcriptions of the C++; compare_tables.py checks the
% transcription against the source file. Every expand(...) must print zeros;
% check_tet10.m asserts it. Runs tet10_generate.m and tet10_shortcuts.m.
tet10_generate ;
tet10_shortcuts ;

% eta and zeta are exchanged relative to the naive node map ( see the pin
% in tet10_generate.m ); consequently their derivative tables have
% exchanged contents.

% faces: the generator emits 3 candidates per face; only the first two are active
act = [1 2 4 5 7 8 10 11];

%% G
mGxi = [xi16-4
zeta8-2
0
0
0
eta8
xi16-4
6-xi16-eta8-zeta8
0
-zeta8
0
-eta8];

expand(diff(g,xi)-mGxi)'

mGeta = [0
0
0
zeta8
eta16-4
xi8-2
0
-xi8
0
-zeta8
eta16-4
6-xi8-eta16-zeta8];
expand(diff(g,eta)-mGeta)'

mGzeta = [0
xi8
zeta16-4
eta8-2
0
0
0
-xi8
zeta16-4
6-xi8-eta8-zeta16
0
-eta8];

expand(diff(g,zeta)-mGzeta)'

%% h
mHxi = [zeta8
0
0
0
eta8-2
xi16-4
10-xi16-eta8-zeta8
eta16+xi16+zeta16-12
2-zeta8
eta16+xi16+zeta16-12
2-eta8
eta16+xi16+zeta16-12];
expand(mHxi-diff(h,xi))'

mHeta = [0
0
zeta8-2
eta16-4
xi8
0
2-xi8
eta16+xi16+zeta16-12
2-zeta8
eta16+xi16+zeta16-12
10-xi8-eta16-zeta8
eta16+xi16+zeta16-12];
expand(mHeta-diff(h,eta))'

mHzeta = [xi8-2
zeta16-4
eta8
0
0
0
2-xi8
eta16+xi16+zeta16-12
10-xi8-eta8-zeta16
eta16+xi16+zeta16-12
2-eta8
eta16+xi16+zeta16-12 ];

expand(mHzeta-diff(h,zeta))'

%% U

mUxi = [-zeta16
zeta8
-eta16
eta8
-eta16
eta8
0
0];

expand(mUxi-diff(u(act),xi))'

mUeta = [-zeta16
zeta8
16-xi16-eta16-eta16-zeta16
zeta8+xi8+eta16-8
16-xi16-eta16-eta16-zeta16
zeta8+xi8+eta16-8
zeta16
-zeta8];

expand(mUeta-diff(u(act),eta))'

mUzeta = [16-xi16-eta16-zeta16-zeta16
zeta16+xi8+eta8-8
-eta16
eta8
-eta16
eta8
eta16
-eta8];
expand(mUzeta-diff(u(act),zeta))'

%% V
mVxi = [zeta8+xi16+eta8-8
16-xi16-xi16-eta16-zeta16
zeta8
-zeta16
-eta8
eta16
-zeta8
zeta16];
expand(mVxi-diff(v(act),xi))'

mVeta = [xi8
-xi16
zeta8
-zeta16
-xi8
xi16
0
0];
expand(mVeta-diff(v(act),eta))'

mVzeta = [xi8
-xi16
zeta16+xi8+eta8-8
16-xi16-eta16-zeta16-zeta16
0
0
-xi8
xi16];
expand(mVzeta-diff(v(act),zeta))'

%% W
mWxi = [-zeta8
-zeta8
0
0
zeta8+xi16+eta8-8
zeta8+xi16+eta8-8
-eta8
-eta8];
expand(mWxi-diff(w(act),xi))'

mWeta = [0
0
-zeta8
-zeta8
xi8
xi8
-xi8
-xi8];
expand(mWeta-diff(w(act),eta))'

mWzeta = [-xi8
-xi8
-eta8
-eta8
xi8
xi8
0
0];

expand(mWzeta-diff(w(act),zeta))'
