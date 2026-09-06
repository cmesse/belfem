//Solver.AutoMesh = 0;
Mesh.Algorithm = 5;
//Geometry.AutoCoherence = 1;
Mesh.Format = 1;
// fine mesh value
mm = 1e-3 ;
cm  =1e-2 ;
lc = 0.5*cm;
lc_air = 5*cm;
air_r=60*cm;
layers = 20 ;
coil_r1 = 4.5*cm ;
coil_r2 = 7*cm ;
coil_heigth = 10*cm ;

// Center point with fine mesh size 
Point(13) = {0, 0, 0, lc};

Point(1) = {0,0,-air_r,lc_air};
Point(2) = {air_r,0,0,lc_air};
Point(3) = {0,0,air_r,lc_air};

Line(1) = {13,1};
Circle(2) = {1,13,2} ;
Circle(3) = {2,13,3};
Line(4) = {3,13} ;

// Lines connecting the points of the magnet
//Line(1) = {1, 2};
//Line(2) = {2, 3};
//Line(3) = {3, 4};
//Line(4) = {4, 5};
//Line(5) = {5, 6};
//Line(6) = {6, 7};
//Line(7) = {7, 8};
//Line(8) = {8, 1};

// Coil points
//Point(4) = {0.5*coil_r1+0.5*coil_r2,0,0} ;
Point(9) = {coil_r1, 0, -coil_heigth/2, lc};
Point(10) = {coil_r2, 0, -coil_heigth/2, lc};
Point(11) = {coil_r2, 0, coil_heigth/2, lc};
Point(12) = {coil_r1, 0, coil_heigth/2, lc};

// Lines for the coil
Line(9) = {9, 10};
Line(10) = {10, 11};
Line(11) = {11, 12};
Line(12) = {12, 9};




// Line loops to generate surfaces
Line Loop(1) = {1, 2, 3, 4,-9,-10,-11,-12}; // Outer boundary
Line Loop(3) = {9, 10, 11, 12}; // Coil boundary

// Define the surface with a hole
Plane Surface(2) = {1};
Plane Surface(4) = {3};

//volIron=Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {Surface{2} ; Layers{layers}; };
volCoil=Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {Surface{2,4} ; };
volCoil2=Extrude {{0, 0, 1}, {0, 0, 0}, -Pi/2} {Surface{2,4}  ; };
volCoil3=Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {Surface{68,46}  ; };
volCoil4=Extrude {{0, 0, 1}, {0, 0, 0}, Pi/2} {Surface{146,180}  ; Recombine ; };
//volAir=Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {Surface{2}  ; };
//volAir2=Extrude {{0, 0, 1}, {0, 0, 0}, -Pi} {Surface{2}  ; };



// Air

//Sphere(1000) = {0, 0, 0, air_r, 0, 2*Pi, 2*Pi};

// Remove volumes from air
//volAirTmp = BooleanDifference{ Volume{1000}; Delete;}{ Volume{1};  };
// Remove overlapping surfaces
//volumes = BooleanFragments{ Volume{volAirTmp}; Delete;}{Volume{1};  Delete;}; 

//Coherence ;
