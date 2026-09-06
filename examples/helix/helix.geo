// ============================================================
// Four helical conductors in cylindrical air domain
// 90° twist over domain length — Pure translation periodic BCs
// ============================================================

// ============================================
// Parameters
// ============================================
r_cond   = 0.5 ;    // Conductor radius
r_offset = 1.0 ;    // Distance of conductor center from z-axis
r_air    = 10.0 ;   // Air domain radius
L        = 5.0 ;    // Length along z
twist    = Pi/2 ;   // Total twist angle (90°)

h_cond0  = 0.1 ;   // Mesh size at conductor 0 (finest)
h_cond1  = 0.1  ;   // Mesh size at conductor 1
h_cond2  = 0.1  ;   // Mesh size at conductor 2
h_cond3  = 0.1  ;   // Mesh size at conductor 3 (coarsest)
h_air    = 1.0  ;   // Mesh size at air boundary

// ============================================
// Points at z = 0 (front face)
// ============================================
// Origin (center of air domain)
Point(1) = {0, 0, 0, h_air};

// Air boundary
Point(2) = { r_air,  0,     0, h_air};
Point(3) = { 0,      r_air, 0, h_air};
Point(4) = {-r_air,  0,     0, h_air};
Point(5) = { 0,     -r_air, 0, h_air};

// Conductor 0 — centered at (r_offset, 0)
Point(6)  = { r_offset,          0,       0, h_cond0};  // center
Point(7)  = { r_offset + r_cond, 0,       0, h_cond0};
Point(8)  = { r_offset,          r_cond,  0, h_cond0};
Point(9)  = { r_offset - r_cond, 0,       0, h_cond0};
Point(10) = { r_offset,         -r_cond,  0, h_cond0};

// Conductor 1 — centered at (0, r_offset)
Point(11) = { 0,      r_offset,          0, h_cond1};  // center
Point(12) = { r_cond, r_offset,          0, h_cond1};
Point(13) = { 0,      r_offset + r_cond, 0, h_cond1};
Point(14) = {-r_cond, r_offset,          0, h_cond1};
Point(15) = { 0,      r_offset - r_cond, 0, h_cond1};

// Conductor 2 — centered at (-r_offset, 0)
Point(16) = {-r_offset,          0,       0, h_cond2};  // center
Point(17) = {-r_offset + r_cond, 0,       0, h_cond2};
Point(18) = {-r_offset,          r_cond,  0, h_cond2};
Point(19) = {-r_offset - r_cond, 0,       0, h_cond2};
Point(20) = {-r_offset,         -r_cond,  0, h_cond2};

// Conductor 3 — centered at (0, -r_offset)
Point(21) = { 0,      -r_offset,          0, h_cond3};  // center
Point(22) = { r_cond, -r_offset,          0, h_cond3};
Point(23) = { 0,      -r_offset + r_cond, 0, h_cond3};
Point(24) = {-r_cond, -r_offset,          0, h_cond3};
Point(25) = { 0,      -r_offset - r_cond, 0, h_cond3};

// ============================================
// Arcs at z = 0
// ============================================
// Air boundary (quarter arcs, center = Point 1)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Conductor 0 (quarter arcs, center = Point 6)
Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

// Conductor 1 (quarter arcs, center = Point 11)
Circle(9)  = {12, 11, 13};
Circle(10) = {13, 11, 14};
Circle(11) = {14, 11, 15};
Circle(12) = {15, 11, 12};

// Conductor 2 (quarter arcs, center = Point 16)
Circle(13) = {17, 16, 18};
Circle(14) = {18, 16, 19};
Circle(15) = {19, 16, 20};
Circle(16) = {20, 16, 17};

// Conductor 3 (quarter arcs, center = Point 21)
Circle(17) = {22, 21, 23};
Circle(18) = {23, 21, 24};
Circle(19) = {24, 21, 25};
Circle(20) = {25, 21, 22};

// ============================================
// Curve loops and plane surfaces (front face)
// ============================================
// Conductor disks
Curve Loop(1) = {5, 6, 7, 8};
Plane Surface(1) = {1};

Curve Loop(2) = {9, 10, 11, 12};
Plane Surface(2) = {2};

Curve Loop(3) = {13, 14, 15, 16};
Plane Surface(3) = {3};

Curve Loop(4) = {17, 18, 19, 20};
Plane Surface(4) = {4};

// Air annulus (outer boundary with 4 conductor holes)
Curve Loop(5) = {1, 2, 3, 4};
Plane Surface(5) = {5, 1, 2, 3, 4};

// ============================================
// Twist extrude — 90° rotation about z-axis
// ============================================
out[] = Extrude { {0,0,L}, {0,0,1}, {0,0,0}, twist } {
    Surface{1, 2, 3, 4, 5};
};

// Output layout (per surface: top, volume, lateral...):
//   Conductor 0: out[0..5]   — top=out[0],  vol=out[1]
//   Conductor 1: out[6..11]  — top=out[6],  vol=out[7]
//   Conductor 2: out[12..17] — top=out[12], vol=out[13]
//   Conductor 3: out[18..23] — top=out[18], vol=out[19]
//   Air:         out[24..45] — top=out[24], vol=out[25]

// ============================================
// Periodic boundary conditions along z
// Pure translation: each front face maps to the back face
// of a DIFFERENT conductor (the one that twisted into its
// position due to 4-fold symmetry + 90° twist)
// ============================================
// Front cond 0 (0°)   → Back cond 3 (270°→360°=0°)
Periodic Surface {out[18]} = {1} Translate {0, 0, L};
// Front cond 1 (90°)  → Back cond 0 (0°→90°)
Periodic Surface {out[0]}  = {2} Translate {0, 0, L};
// Front cond 2 (180°) → Back cond 1 (90°→180°)
Periodic Surface {out[6]}  = {3} Translate {0, 0, L};
// Front cond 3 (270°) → Back cond 2 (180°→270°)
Periodic Surface {out[12]} = {4} Translate {0, 0, L};
// Air
Periodic Surface {out[24]} = {5} Translate {0, 0, L};

// ============================================
// Entity summary
// ============================================
Printf("--- Volumes ---");
Printf("  Conductor 0 : %g", out[1]);
Printf("  Conductor 1 : %g", out[7]);
Printf("  Conductor 2 : %g", out[13]);
Printf("  Conductor 3 : %g", out[19]);
Printf("  Air          : %g", out[25]);
Printf("--- Front faces (z = 0) ---");
Printf("  Conductor 0 : 1");
Printf("  Conductor 1 : 2");
Printf("  Conductor 2 : 3");
Printf("  Conductor 3 : 4");
Printf("  Air          : 5");
Printf("--- Back faces (z = L) ---");
Printf("  Conductor 0 : %g", out[0]);
Printf("  Conductor 1 : %g", out[6]);
Printf("  Conductor 2 : %g", out[12]);
Printf("  Conductor 3 : %g", out[18]);
Printf("  Air          : %g", out[24]);
