// ============================================================
// 1 layer CORC with translated periodic BC
// ============================================================

// ============================================
// Parameters
// ============================================

r_air    = 20.0 ;   // Air domain radius


pitch_angle = Pi/4;
num_pitch = 1/3;
tape_width = 4.0; // In mm
gap = 0.1; // In mm



h_cond0  = 0.05 ;   // Mesh size at conductor 0 (finest)
h_cond1  = 0.25 ;   // Mesh size at conductor 0 (coarse)
h_air    = 2  ;   // Mesh size at air boundary

// ============================================
// Points at z = 0 (front face)
// ============================================
// Origin (center of air domain)
Point(1) = {0, 0, 0, h_air/10};

// Air boundary
Point(2) = { r_air,  0,     0, h_air};
Point(3) = { r_air*Cos(2*Pi/3),      r_air*Sin(2*Pi/3), 0, h_air};
Point(4) = { r_air*Cos(4*Pi/3),      r_air*Sin(4*Pi/3), 0, h_air};
//Point(5) = { 0,     -r_air, 0, h_air};

// CORC points
r = 3*(tape_width + gap)/(2*Pi*Cos(pitch_angle)) ;
P = 2*Pi*r/Tan(pitch_angle) ;
Point(6)  = { r,          0,       0, h_cond0};
Point(7)  = { r*Cos(0.5*tape_width/(r*Cos(pitch_angle))),          r*Sin(0.5*tape_width/(r*Cos(pitch_angle))),       0, h_cond1};
Point(8)  = { r*Cos(tape_width/(r*Cos(pitch_angle))),          r*Sin(tape_width/(r*Cos(pitch_angle))),       0, h_cond0};
Point(9)  = { r*Cos((tape_width+gap)/(r*Cos(pitch_angle))),          r*Sin((tape_width+gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(10)  = { r*Cos((1.5*tape_width+gap)/(r*Cos(pitch_angle))),          r*Sin((1.5*tape_width+gap)/(r*Cos(pitch_angle))),       0, h_cond1};
Point(11)  = { r*Cos((2*tape_width+gap)/(r*Cos(pitch_angle))),          r*Sin((2*tape_width+gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(12)  = { r*Cos((2*tape_width+2*gap)/(r*Cos(pitch_angle))),          r*Sin((2*tape_width+2*gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(13)  = { r*Cos((2.5*tape_width+2*gap)/(r*Cos(pitch_angle))),          r*Sin((2.5*tape_width+2*gap)/(r*Cos(pitch_angle))),       0, h_cond1};
Point(14)  = { r*Cos((3*tape_width+2*gap)/(r*Cos(pitch_angle))),          r*Sin((3*tape_width+2*gap)/(r*Cos(pitch_angle))),       0, h_cond0};

// ============================================
// Arcs at z = 0
// ============================================
// Air boundary (quarter arcs, center = Point 1)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 2};
//Circle(4) = {5, 1, 2};

// Conductor 0 (quarter arcs, center = Point 6)
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 10};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 12};
Circle(11) = {12, 1, 13};
Circle(12) = {13, 1, 14};
Circle(13) = {14, 1, 6};


// ============================================
// Curve loops and plane surfaces (front face)
// ============================================
// Conductor disks
Curve Loop(1) = {5, 6, 7, 8, 9, 10,11,12,13};
Plane Surface(1) = {1};
Point{1} In Surface{1} ;


// Air annulus (outer boundary with 4 conductor holes)
Curve Loop(2) = {1, 2, 3,-5,-6,-7,-8,-9,-10,-11,-12,-13};
Plane Surface(2) = {2};

// ============================================
// Twist extrude — 90° rotation about z-axis
// ============================================
out[] = Extrude { {0,0,P*num_pitch}, {0,0,1}, {0,0,0}, 2*Pi*num_pitch } {
    Surface{1, 2} ; Curve{5,6,7,8,9,10,11,12,13};
};
Point{16} In Surface{60} ;

tx = 0; ty = 0; tz = P*num_pitch; // Translation vector


Periodic Surface {60} = {1} Translate {0, 0, tz};
Periodic Surface {122} = {2} Translate {0, 0, tz};
