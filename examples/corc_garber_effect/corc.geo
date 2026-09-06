// ============================================================
// 1 layer CORC with twisted periodic BC
// ============================================================

// ============================================
// Parameters
// ============================================

r_air    = 20.0 ;   // Air domain radius


pitch_angle = Pi/4;
num_pitch = 0.01;
tape_width = 4.0; // In mm
gap = 0.01; // In mm



h_cond0  = 0.01 ;   // Mesh size at conductor 0 (finest)
h_air    = 2  ;   // Mesh size at air boundary

// ============================================
// Points at z = 0 (front face)
// ============================================
// Origin (center of air domain)
Point(1) = {0, 0, 0, h_air/2};

// Air boundary
Point(2) = { r_air,  0,     0, h_air};
Point(3) = { 0,      r_air, 0, h_air};
Point(4) = {-r_air,  0,     0, h_air};
Point(5) = { 0,     -r_air, 0, h_air};

// CORC points
r = 3*(tape_width + gap)/(2*Pi*Cos(pitch_angle)) ;
P = 2*Pi*r/Tan(pitch_angle) ;
Point(6)  = { r,          0,       0, h_cond0};
Point(7)  = { r*Cos(tape_width/(r*Cos(pitch_angle))),          r*Sin(tape_width/(r*Cos(pitch_angle))),       0, h_cond0};
Point(8)  = { r*Cos((tape_width+gap)/(r*Cos(pitch_angle))),          r*Sin((tape_width+gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(9)  = { r*Cos((2*tape_width+gap)/(r*Cos(pitch_angle))),          r*Sin((2*tape_width+gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(10)  = { r*Cos((2*tape_width+2*gap)/(r*Cos(pitch_angle))),          r*Sin((2*tape_width+2*gap)/(r*Cos(pitch_angle))),       0, h_cond0};
Point(11)  = { r*Cos((3*tape_width+2*gap)/(r*Cos(pitch_angle))),          r*Sin((3*tape_width+2*gap)/(r*Cos(pitch_angle))),       0, h_cond0};

// ============================================
// Arcs at z = 0
// ============================================
// Air boundary (quarter arcs, center = Point 1)
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};

// Conductor 0 (quarter arcs, center = Point 6)
Circle(5) = {6, 1, 7};
Circle(6) = {7, 1, 8};
Circle(7) = {8, 1, 9};
Circle(8) = {9, 1, 10};
Circle(9) = {10, 1, 11};
Circle(10) = {11, 1, 6};


// ============================================
// Curve loops and plane surfaces (front face)
// ============================================
// Conductor disks
Curve Loop(1) = {5, 6, 7, 8, 9, 10};
Plane Surface(1) = {1};
Point{1} In Surface{1} ;


// Air annulus (outer boundary with 4 conductor holes)
Curve Loop(2) = {1, 2, 3, 4,-5,-6,-7,-8,-9,-10};
Plane Surface(2) = {2};

// ============================================
// Twist extrude — 90° rotation about z-axis
// ============================================
out[] = Extrude { {0,0,P*num_pitch}, {0,0,1}, {0,0,0}, 2*Pi*num_pitch } {
    Surface{1, 2} ; Curve{5,6,7,8,9,10};
};
Point{13} In Surface{42} ;
// Define Rotation and Translation values
angle = 2*Pi*num_pitch; // e.g., 90 degrees
ax = 0; ay = 0; az = 1; // Axis of rotation
px = 0; py = 0; pz = 0; // Center of rotation
tx = 0; ty = 0; tz = P*num_pitch; // Translation vector


Periodic Surface {42} = {1} Affine {
  Cos(angle) + ax*ax*(1-Cos(angle)), ax*ay*(1-Cos(angle)) - az*Sin(angle), ax*az*(1-Cos(angle)) + ay*Sin(angle), tx,
  ay*ax*(1-Cos(angle)) + az*Sin(angle), Cos(angle) + ay*ay*(1-Cos(angle)), ay*az*(1-Cos(angle)) - ax*Sin(angle), ty,
  az*ax*(1-Cos(angle)) - ay*Sin(angle), az*ay*(1-Cos(angle)) + ax*Sin(angle), Cos(angle) + az*az*(1-Cos(angle)), tz
};
Periodic Surface {94} = {2} Affine {
  Cos(angle) + ax*ax*(1-Cos(angle)), ax*ay*(1-Cos(angle)) - az*Sin(angle), ax*az*(1-Cos(angle)) + ay*Sin(angle), tx,
  ay*ax*(1-Cos(angle)) + az*Sin(angle), Cos(angle) + ay*ay*(1-Cos(angle)), ay*az*(1-Cos(angle)) - ax*Sin(angle), ty,
  az*ax*(1-Cos(angle)) - ay*Sin(angle), az*ay*(1-Cos(angle)) + ax*Sin(angle), Cos(angle) + az*az*(1-Cos(angle)), tz
};
