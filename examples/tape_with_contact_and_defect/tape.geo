// One REBCO tape in a periodic air cylinder, derived from tapestack3d.geo.
//
// The tape lies in the plane y = 0, spans x in [ -tapeWidth/2, tapeWidth/2 ]
// and runs along z from 0 to domainLength. The plane y = 0 splits the air
// into an upper ( volume 1 ) and a lower ( volume 2 ) half, so the tape is a
// shared face of both.
//
// Entity ids used by input.conf:
//   surfaces 1, 2     tape, left and right half ( thin shell )
//   surfaces 9, 10    air end faces at z = 0, upper and lower
//   surfaces 11, 12   air end faces at z = domainLength, upper and lower
//   points 2, 3, 5    periodic source ( z = 0 )
//   points 9, 10, 12  periodic target ( z = domainLength )
//
// Point ids double as node ids in the mesh: every point is meshed and gmsh
// numbers the point nodes first, contiguously. Keep the ids gap-free.

tapeWidth    = 4 ;
domainRadius = 25 ;
domainLength = 10 ;

tapeResolutionTip    = 0.1 ;
tapeResolutionCenter = 0.5 ;
domainResolution     = 5 ;

// refinement around the defect in src/matlib.cpp ( centered at z = domainLength/2 )
defectResolution = 0.15 ;
defectHalfLength = 2.0 ;

// --- points -------------------------------------------------------------
For k In {0:1}
	z = k*domainLength ;
	p = 7*k ;
	Point(p+1) = { 0,               0,             z, tapeResolutionCenter } ;
	Point(p+2) = { -0.5*tapeWidth,  0,             z, tapeResolutionTip } ;
	Point(p+3) = {  0.5*tapeWidth,  0,             z, tapeResolutionTip } ;
	Point(p+4) = {  domainRadius,   0,             z, domainResolution } ;
	Point(p+5) = { 0,               domainRadius,  z, domainResolution } ;
	Point(p+6) = { -domainRadius,   0,             z, domainResolution } ;
	Point(p+7) = { 0,              -domainRadius,  z, domainResolution } ;
EndFor

// --- curves -------------------------------------------------------------
Line(1)  = { 2, 1 } ;      // tape front edge, left
Line(2)  = { 1, 3 } ;      // tape front edge, right
Line(3)  = { 9, 8 } ;      // tape back edge, left
Line(4)  = { 8, 10 } ;     // tape back edge, right
Line(5)  = { 2, 9 } ;      // tape tip, x = -tapeWidth/2
Line(6)  = { 1, 8 } ;      // tape center line
Line(7)  = { 3, 10 } ;     // tape tip, x = +tapeWidth/2
Line(8)  = { 3, 4 } ;
Line(9)  = { 6, 2 } ;
Line(10) = { 10, 11 } ;
Line(11) = { 13, 9 } ;
Line(12) = { 4, 11 } ;
Line(13) = { 6, 13 } ;

Circle(14) = { 4, 1, 5 } ;
Circle(15) = { 5, 1, 6 } ;
Circle(16) = { 6, 1, 7 } ;
Circle(17) = { 7, 1, 4 } ;
Circle(18) = { 11, 8, 12 } ;
Circle(19) = { 12, 8, 13 } ;
Circle(20) = { 13, 8, 14 } ;
Circle(21) = { 14, 8, 11 } ;

Line(22) = { 5, 12 } ;
Line(23) = { 7, 14 } ;

// --- surfaces -----------------------------------------------------------
// both tape halves are traversed in the same sense, so their normals agree
Curve Loop(1) = { 1, 6, -3, -5 } ;
Plane Surface(1) = {1} ;
Curve Loop(2) = { 2, 7, -4, -6 } ;
Plane Surface(2) = {2} ;

// air-air parts of the plane y = 0
Curve Loop(3) = { 8, 12, -10, -7 } ;
Plane Surface(3) = {3} ;
Curve Loop(4) = { 9, 5, -11, -13 } ;
Plane Surface(4) = {4} ;

// outer cylinder, one quadrant each
Curve Loop(5) = { 14, 22, -18, -12 } ;
Surface(5) = {5} ;
Curve Loop(6) = { 15, 13, -19, -22 } ;
Surface(6) = {6} ;
Curve Loop(7) = { 16, 23, -20, -13 } ;
Surface(7) = {7} ;
Curve Loop(8) = { 17, 12, -21, -23 } ;
Surface(8) = {8} ;

// end faces
Curve Loop(9) = { 9, 1, 2, 8, 14, 15 } ;
Plane Surface(9) = {9} ;
Curve Loop(10) = { 9, 1, 2, 8, -17, -16 } ;
Plane Surface(10) = {10} ;
Curve Loop(11) = { 11, 3, 4, 10, 18, 19 } ;
Plane Surface(11) = {11} ;
Curve Loop(12) = { 11, 3, 4, 10, -21, -20 } ;
Plane Surface(12) = {12} ;

// --- volumes ------------------------------------------------------------
Surface Loop(1) = { 1, 2, 3, 4, 5, 6, 9, 11 } ;
Volume(1) = {1} ;
Surface Loop(2) = { 1, 2, 3, 4, 7, 8, 10, 12 } ;
Volume(2) = {2} ;

// --- mesh periodicity ---------------------------------------------------
// The solver pairs the z = 0 and z = domainLength faces facet by facet.
Periodic Surface { 11 } = { 9 }  Translate { 0, 0, domainLength } ;
Periodic Surface { 12 } = { 10 } Translate { 0, 0, domainLength } ;

// --- defect refinement --------------------------------------------------
Field[1] = Box ;
Field[1].VIn  = defectResolution ;
Field[1].VOut = domainResolution ;
Field[1].XMin = -0.5*tapeWidth - 0.5 ;
Field[1].XMax =  0.5*tapeWidth + 0.5 ;
Field[1].YMin = -0.5 ;
Field[1].YMax =  0.5 ;
Field[1].ZMin = 0.5*domainLength - defectHalfLength ;
Field[1].ZMax = 0.5*domainLength + defectHalfLength ;
Field[1].Thickness = 2 ;
Background Field = 1 ;
