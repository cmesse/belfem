// domain radius
R = 2000 ;

// resolution for domain edge
re = 200 ;

// resolution for symmetry line
rs = 5 ;

// resolution along yoke
ry0 = 25 ;
ry1 = 1 ;
ry2 = 1 ;

// radius of the rounded yoke corners. The sharp pole-face corners are a
// field-singularity source; rounding them bounds the peak. The binding
// constraint is the corner at point 5, whose segment to point 18 is only
// 1.5 mm long, so the tangent offset there must stay well under that.
r_fillet = 0.5 ;

// width of one tape
w = 4.6 ;

// resolution along coil
rc =  w/20; 
rcm = 2*rc ;

// number of tapes -1 ( 58 )
n = 57 ;

// length of one coil
lc = 23.6 ;

// distance between tape
d = lc/n ;


// domain points

yoke_x0 = 381.1302 ;
yoke_x1 = 59.7 ;
yoke_x2 = 314.9302 ;
yoke_x3 = 30 ;


coil_x1 = 31.5 ;
coil_x2 = 336.4302 ;
coil_y1 = 29 ;
coil_y2 = 34.5 ;
coil_y3 = 42.5 ;
coil_y4 = 48 ;
coil_y5 = 29 ;
coil_y6 = 34.5 ;
coil_y7 = 42.5 ;
coil_y8 = 48 ;

// offset of center point
offx =  yoke_x0 ;

Point(1) = { offx , 0, 0, ry2 } ;
Point(2) = { R+offx , 0, 0, re } ;
Point(3) = { offx , R, 0, re } ;
Point(4) = { offx-R, 0, 0, re } ;

// domain boundary
Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
l = 2 ;
b = 0 ;

// yoke points
Point( 5 ) = { 30, 52, 0, ry2 };
Point( 6 ) = { 30, 228.936375, 0, ry0 };
Point( 7 ) = { 289.9618304, 371.70178, 0, ry0 };
Point( 8 ) = { 481.0889124, 371.70178, 0, ry0 };
Point( 9 ) = { 686.6666667, 143.1956775, 0, ry0 };
Point( 10 ) = { 686.6666667, 0, 0, ry0 };
Point( 11 ) = { yoke_x0, 52, 0, ry2 };
Point( 12 ) = { 203.318856, 52, 0, ry1 };
Point( 13 ) = { 91.65341172, 25, 0, ry1 };
Point( 14 ) = { yoke_x1, 25, 0, ry1 };
Point( 15 ) = { yoke_x2, 52, 0, ry2 };
Point( 16 ) = { yoke_x1, 52, 0, ry2 };
Point( 17 ) = { yoke_x0, coil_y5, 0, ry2 };
Point( 18 ) = { coil_x1, 52, 0, ry2 };
Point( 19 ) = { coil_x1+lc, 52, 0, ry2 };
Point( 20 ) = { coil_x2, 52, 0, ry2 };
Point( 21 ) = { coil_x2+lc, 52, 0, ry2 };

// anchor point
Point( 22 ) = { coil_x1, 0, 0, rs };

// point for cut at coil 4
Point( 23 ) = { coil_x1-coil_y4, 0, 0, rs };

// point for cut at coil 3
Point( 24 ) = { coil_x1-coil_y3, 0, 0, rs };

// point for cut at coil 2
Point( 25 ) = { coil_x1-coil_y2, 0, 0, rs };

// point for cut at coil 1
Point( 26 ) = { coil_x1-coil_y1, 0, 0, rs };

// anchor point
Point( 27 ) = { coil_x2, 0, 0, rs };

// point for cut at coil 5
Point( 28 ) = { coil_x2-coil_y5, 0, 0, rs };

// point for cut at coil 6
Point( 29 ) = { coil_x2-coil_y6, 0, 0, rs };

// point for cut at coil 7
Point( 30 ) = { coil_x2-coil_y7, 0, 0, rs };

// point for cut at coil 8
Point( 31 ) = { coil_x2-coil_y8, 0, 0, rs };

//------------------------------------------
// rounded yoke corners
//------------------------------------------
// Each entry is a corner P with its two neighbours A and B, in yoke
// traversal order. Point 15 is deliberately ABSENT: it lies at
// ( yoke_x2, 52 ) between two segments that are both at y = 52, so it is a
// waypoint on the straight top face, not a corner - its fillet would be
// degenerate ( uA + uB = 0 ).
//
// Corners 16, 12 and 11 are re-entrant ( interior 270, 193.6, 270 deg ) and
// 5, 14, 13 are convex. The construction below needs no branch for that: the
// centre always sits on the bisector of the wedge between the two edges, and
// is exactly tangent to both either way. Convexity only decides whether the
// arc trims iron ( convex ) or fills the notch ( re-entrant ), which is the
// wanted behaviour in both cases.

p = 31 ;

fax[] = { 30         , coil_x1+lc , yoke_x1     , yoke_x1     , 91.65341172 , coil_x2+lc };
fay[] = { 228.936375 , 52         , 52          , 25          , 25          , 52         };
fpx[] = { 30         , yoke_x1    , yoke_x1     , 91.65341172 , 203.318856  , yoke_x0    };
fpy[] = { 52         , 52         , 25          , 25          , 52          , 52         };
fbx[] = { coil_x1    , yoke_x1    , 91.65341172 , 203.318856  , yoke_x2     , yoke_x0    };
fby[] = { 52         , 25         , 25          , 52          , 52          , coil_y5    };

ft1[] = {} ;   // tangent point on the A side
ft2[] = {} ;   // tangent point on the B side
fcc[] = {} ;   // arc centre

For i In {0:5}
    fuax = fax[i] - fpx[i] ;   fuay = fay[i] - fpy[i] ;
    fubx = fbx[i] - fpx[i] ;   fuby = fby[i] - fpy[i] ;

    fla = Sqrt( fuax*fuax + fuay*fuay ) ;
    flb = Sqrt( fubx*fubx + fuby*fuby ) ;

    fuax /= fla ;  fuay /= fla ;
    fubx /= flb ;  fuby /= flb ;

    ftheta = Acos( fuax*fubx + fuay*fuby ) ;

    fd = r_fillet / Tan( 0.5*ftheta ) ;      // corner to tangent point
    fe = r_fillet / Sin( 0.5*ftheta ) ;      // corner to arc centre

    If ( fd >= fla || fd >= flb )
        Error( "r_fillet = %g is too large for corner %g: tangent offset %g exceeds an adjacent segment ( %g, %g )",
               r_fillet, i, fd, fla, flb ) ;
        Abort ;
    EndIf

    fmx = fuax + fubx ;   fmy = fuay + fuby ;
    flm = Sqrt( fmx*fmx + fmy*fmy ) ;
    fmx /= flm ;  fmy /= flm ;

    p += 1 ;  Point( p ) = { fpx[i] + fd*fuax, fpy[i] + fd*fuay, 0, ry1 } ;  ft1[] += p ;
    p += 1 ;  Point( p ) = { fpx[i] + fd*fubx, fpy[i] + fd*fuby, 0, ry1 } ;  ft2[] += p ;
    p += 1 ;  Point( p ) = { fpx[i] + fe*fmx,  fpy[i] + fe*fmy,  0, ry1 } ;  fcc[] += p ;
EndFor


// Symmetry Line
Line(3) = {4, 23};
Line(4) = {23, 24};
Line(5) = {24, 25};
Line(6) = {25, 26};
Line(7) = {26, 22};
Line(8) = {22, 31};
Line(9) = {31, 30};
Line(10) = {30, 29};
Line(11) = {29, 28};
Line(12) = {28, 27};
Line(13) = {27, 1};
Line(14) = {1, 10};
Line(15) = {10, 2};

// yoke
Line(16) = {10, 9};
Line(17) = {9, 8};
Line(18) = {8, 7};
Line(19) = {7, 6};
// the corner points 5, 16, 14, 13, 12 and 11 are rounded, so these segments
// stop at the tangent points instead. Point 15 is untouched ( collinear ).
Line(20) = {6, ft1[0]};
Line(21) = {ft2[0], 18};
Line(22) = {18, 19};
Line(23) = {19, ft1[1]};
Line(24) = {ft2[1], ft1[2]};
Line(25) = {ft2[2], ft1[3]};
Line(26) = {ft2[3], ft1[4]};
Line(27) = {ft2[4], 15};
Line(28) = {15, 20};
Line(29) = {20, 21};
Line(30) = {21, ft1[5]};
Line(31) = {ft2[5], 17};
Line(32) = {17, 1};
l = 32 ;

// the corner arcs. Appended here so that every downstream coil curve id
// shifts by exactly the number of arcs, automatically, via the same counter.
l += 1 ;  Circle(l) = { ft1[0], fcc[0], ft2[0] } ;  a5  = l ;
l += 1 ;  Circle(l) = { ft1[1], fcc[1], ft2[1] } ;  a16 = l ;
l += 1 ;  Circle(l) = { ft1[2], fcc[2], ft2[2] } ;  a14 = l ;
l += 1 ;  Circle(l) = { ft1[3], fcc[3], ft2[3] } ;  a13 = l ;
l += 1 ;  Circle(l) = { ft1[4], fcc[4], ft2[4] } ;  a12 = l ;
l += 1 ;  Circle(l) = { ft1[5], fcc[5], ft2[5] } ;  a11 = l ;

// the original corner points are now referenced by nothing. Left in place
// they would still be written as isolated 0-D elements, i.e. stray vertices
// with no curve through them, which the mesh reader has to filter. Drop
// them. Point 15 stays: it is still an endpoint of lines 27 and 28. The six
// arc centres cannot be dropped ( each Circle owns its centre ).
Delete { Point{ 5, 16, 14, 13, 12, 11 }; }


//------------------------------------------
// COIL 4
//------------------------------------------

x = coil_x1 ;
y = coil_y4 ;

p4a = p+1 ;
p4b = p+2 ;
p4c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p4d = p-2 ;
p4e = p-1 ;
p4f = p ;

// horizontal lines upper
q = p4a ;
l4hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {23, 22, p4b};
cut4 = l ;

// horizontal lines cut
q = p4b ;
l4hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p4c ;
l4hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l4v0 = l ;
Line( l ) = { 18, p4a };
l += 1 ;
Line( l ) = { p4d, 19 };

// vertical lines coil
q = p4a ;
l4v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l4v2 = l-1 ;

//------------------------------------------
// COIL 3
//------------------------------------------

x = coil_x1 ;
y = coil_y3 ;

p3a = p+1 ;
p3b = p+2 ;
p3c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p3d = p-2 ;
p3e = p-1 ;
p3f = p ;

// horizontal lines upper
q = p3a ;
l3hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {24, 22, p3b};
cut3 = l ;

// horizontal lines cut
q = p3b ;
l3hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p3c ;
l3hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l3v0 = l ;
Line( l ) = { p4c, p3a };
l += 1 ;
Line( l ) = { p4f, p3d };

// vertical lines coil
q = p3a ;
l3v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l3v2 = l-1 ;

//------------------------------------------
// COIL 2
//------------------------------------------

x = coil_x1 ;
y = coil_y2 ;

// points for air block between coils 2 and 3
p2m = p+1 ;
p+=1;
Point( p ) = {x, 0.5*(coil_y2+coil_y3), 0, rcm };
p+=1;
Point( p ) = {x+lc, 0.5*(coil_y2+coil_y3), 0, rcm };
l+=1 ;
Line( l ) = {p2m, p2m+1};
l2ha = l ;

p2a = p+1 ;
p2b = p+2 ;
p2c = p+3 ;

// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p2d = p-2 ;
p2e = p-1 ;
p2f = p ;

// horizontal lines upper
q = p2a ;
l2hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {25, 22, p2b};
cut2 = l ;

// horizontal lines cut
q = p2b ;
l2hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p2c ;
l2hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l2v0 = l ;
Line( l ) = { p3c, p2m };
l += 1 ;
Line( l ) = { p2m, p2a };

l += 1 ;
Line( l ) = { p3f, p2m+1 };
l += 1 ;
Line( l ) = { p2m+1, p2d };

// vertical lines coil
q = p2a ;
l2v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l2v2 = l-1 ;

//------------------------------------------
// COIL 1
//------------------------------------------

// coil
x = coil_x1;
y = coil_y1;

p1a = p+1 ;
p1b = p+2 ;
p1c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p1d = p-2 ;
p1e = p-1 ;
p1f = p ;

// horizontal lines upper
q = p1a ;
l1hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {26, 22, p1b};
cut1 = l ;

// horizontal lines cut
q = p1b ;
l1hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p1c ;
l1hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l1v0 = l ;
Line( l ) = { p2c, p1a };
l += 1 ;
Line( l ) = { p2f, p1d };

// vertical lines coil
q = p1a ;
l1v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l1v2 = l-1 ;

//------------------------------------------
// COIL 8
//------------------------------------------

x = coil_x2 ;
y = coil_y8 ; 

p8a = p+1 ;
p8b = p+2 ;
p8c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p8d = p-2 ;
p8e = p-1 ;
p8f = p ;

// horizontal lines upper
q = p8a ;
l8hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {31, 27, p8b};
cut8 = l ;

// horizontal lines cut
q = p8b ;
l8hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p8c ;
l8hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l8v0 = l ;
Line( l ) = { 20, p8a };
l += 1 ;
Line( l ) = { p8d, 21 };

// vertical lines coil
q = p8a ;
l8v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l8v2 = l-1 ;

//------------------------------------------
// COIL 7
//------------------------------------------

x = coil_x2 ;
y = coil_y7 ; 

p7a = p+1 ;
p7b = p+2 ;
p7c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p7d = p-2 ;
p7e = p-1 ;
p7f = p ;

// horizontal lines upper
q = p7a ;
l7hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {30, 27, p7b};
cut7 = l ;

// horizontal lines cut
q = p7b ;
l7hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p7c ;
l7hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l7v0 = l ;
Line( l ) = { p8c, p7a };
l += 1 ;
Line( l ) = { p8f, p7d };

// vertical lines coil
q = p7a ;
l7v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l7v2 = l-1 ;

//------------------------------------------
// COIL 6
//------------------------------------------

x = coil_x2 ;
y = coil_y6 ; 

// points for air block between coils 2 and 3
p6m = p+1 ;
p+=1;
Point( p ) = {x, 0.5*(coil_y6+coil_y7), 0, rcm };
p+=1;
Point( p ) = {x+lc, 0.5*(coil_y6+coil_y7), 0, rcm };
l+=1 ;
Line( l ) = {p6m, p6m+1};
l6ha = l ;

p2a = p+1 ;
p2b = p+2 ;
p2c = p+3 ;

p6a = p+1 ;
p6b = p+2 ;
p6c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p6d = p-2 ;
p6e = p-1 ;
p6f = p ;

// horizontal lines upper
q = p6a ;
l6hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {29, 27, p6b};
cut6 = l ;

// horizontal lines cut
q = p6b ;
l6hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p6c ;
l6hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l6v0 = l +1 ;
l += 1 ;
Line( l ) = { p7c, p6m };
l += 1 ;
Line( l ) = { p6m, p6a };
l += 1 ;
Line( l ) = { p7f, p6m+1 };
l += 1 ;
Line( l ) = { p6m+1, p6d };

// vertical lines coil
q = p6a ;
l6v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l6v2 = l-1 ;

//------------------------------------------
// COIL 5
//------------------------------------------

x = coil_x2 ;
y = coil_y5 ; 

p5a = p+1 ;
p5b = p+2 ;
p5c = p+3 ;


// point numbering scheme
//

//  a    d
//  b    e
//  c    f

// points for coil
For k In {1:n+1}
	p += 1 ;
	Point( p ) = {x, y+0.5*w, 0, rc };
	p += 1 ;
	Point( p ) = {x, y, 0, rc };
	p += 1 ;
	Point( p ) = {x, y-0.5*w, 0, rc };
	x += d;
EndFor

p5d = p-2 ;
p5e = p-1 ;
p5f = p ;

// horizontal lines upper
q = p5a ;
l5hu = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// round cut
l+= 1 ;
Circle(l) = {28, 27, p5b};
cut5 = l ;

// horizontal lines cut
q = p5b ;
l5hm = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// horizontal lines lower
q = p5c ;
l5hl = l+1 ;
For k In {1:n}
	l += 1 ;
	Line( l ) = { q, q+3 };
	q += 3 ;
EndFor

// vertical lines air block
l += 1 ;
l5v0 = l ;
Line( l ) = { p6c, p5a };
l += 1 ;
Line( l ) = { p6f, p5d };

// vertical lines coil
q = p5a ;
l5v1 = l+1 ;
For k In {0:n}
	l += 1 ;
	Line( l ) = { q, q+1 };
	l += 1 ;
	Line( l ) = { q+1, q+2 };
	q += 3 ;
EndFor

// remember last line of air block line
l5v2 = l-1 ;

//------------------------------------------
// Blocks
//------------------------------------------

b = 0 ;

b+=1;
Curve Loop(b) = {14, 16, 17, 18, 19, 20, a5, 21, 22, 23, a16, 24, a14, 25, a13, 26, a12, 27, 28, 29, 30, a11, 31, 32};
Plane Surface(b) = {b};

b += 1 ;
Curve Loop(b) = {15, 1, 2, 3, cut4, -l4v0-2, -l4v0, -21, -a5, -20, -19, -18, -17, -16};
Plane Surface(b) = {b};

//------------------------------------------
// air blocks with round cuts
//------------------------------------------

b+=1;
Curve Loop(b) = {4, cut3, -l3v0-2, -l3v0, -l4v0-3, -cut4};
Plane Surface(b) = {b};


b+=1;
Curve Loop(b) = {5, cut2, -l2v0-4, -l2v0-1, -l2v0, -l3v0-3, -cut3};
Plane Surface(b) = {b};

b+=1;
Curve Loop(b) = {6, cut1, -l1v0-2, -l1v0, -l2v0-5, -cut2};
Plane Surface(b) = {b};

b+=1;
Curve Loop(b) = {7, 8, cut8, -l8v0-2, -l8v0, -28, -27, -a12, -26, -a13, -25, -a14, -24, -a16, -23, -l4v0-1, l4v2, l4v2+1, l3v0+1, l3v2, l3v2+1, l2v0+2, l2v0+3, l2v2, l2v2+1, l1v0+1, l1v2, l1v2+1, -l1hl-n+1:-l1hl, -l1v0-3, -cut1};
Plane Surface(b) = {b};

b+=1;
Curve Loop(b) = {9, cut7, -l7v0-2, -l7v0, -l8v0-3, -cut8};
Plane Surface(b) = {b};


b+=1;
Curve Loop(b) = {10, cut6, -l6v0-4, -l6v0-1, -l6v0, -l7v0-3, -cut7};
Plane Surface(b) = {b};

b+=1;
Curve Loop(b) = {11, cut5, -l5v0-2, -l5v0, -l6v0-5, -cut6};
Plane Surface(b) = {b};

b+=1;
Curve Loop(b) = {12, 13, -32, -31, -a11, -30, -l8v0-1, l8v2, l8v2+1, l7v0+1, l7v2, l7v2+1, l6v0+2, l6v0+3, l6v2, l6v2+1, l5v0+1, l5v2, l5v2+1, -l5hl-n+1:-l5hl, -l5v0-3, -cut5};
Plane Surface(b) = {b};

//------------------------------------------

// air block between coil 4 and yoke
b += 1 ;
Curve Loop(b) = {l4v0, l4hu:l4hu+n-1, l4v0+1, -22};
Plane Surface(b) = {b};

// blocks coil 4
lv = l4v0+2 ;
lhu = l4hu ;
lhm = l4hm ;
lhl = l4hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air blocks between coil 3 and 4
b += 1 ;
Curve Loop(b) = {l3v0, l3hu:l3hu+n-1, -l3v0-1, -l4hl-n+1:-l4hl};
Plane Surface(b) = {b};

// blocks coil 3
lv = l3v0+2 ;
lhu = l3hu ;
lhm = l3hm ;
lhl = l3hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 2 and 3
b += 1 ;
Curve Loop(b) = {l2v0, l2ha, -l2v0-2, -l3hl-n+1:-l3hl};
Plane Surface(b) = {b};

b += 1 ;
Curve Loop(b) = {l2v0+1, l2hu:l2hu+n-1, -l2v0-3, -l2ha};
Plane Surface(b) = {b};

// blocks coil 2
lv = l2v0+4 ;
lhu = l2hu ;
lhm = l2hm ;
lhl = l2hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 1 and 2
b += 1 ;
Curve Loop(b) = {l1v0, l1hu:l1hu+n-1, -l1v0-1, -l2hl-n+1:-l2hl};
Plane Surface(b) = {b};

// blocks coil 2
lv = l1v0+2 ;
lhu = l1hu ;
lhm = l1hm ;
lhl = l1hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 8 and yoke
b += 1 ;
Curve Loop(b) = {l8v0, l8hu:l8hu+n-1, l8v0+1, -29};
Plane Surface(b) = {b};

// blocks coil 8
lv = l8v0+2 ;
lhu = l8hu ;
lhm = l8hm ;
lhl = l8hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 7 and 8
b += 1 ;
Curve Loop(b) = {l7v0, l7hu:l7hu+n-1, -l7v0-1, -l8hl-n+1:-l8hl};
Plane Surface(b) = {b};

// blocks coil 7
lv = l7v0+2 ;
lhu = l7hu ;
lhm = l7hm ;
lhl = l7hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 6 and 7
// b += 1 ;
// Curve Loop(b) = {l6v0, l6hu:l6hu+n-1, -l6v0-1, -l7hl-n+1:-l7hl};
// Plane Surface(b) = {b};

// air block between coil 2 and 3
b += 1 ;
Curve Loop(b) = {l6v0, l6ha, -l6v0-2, -l7hl-n+1:-l7hl};
Plane Surface(b) = {b};

b += 1 ;
Curve Loop(b) = {l6v0+1, l6hu:l6hu+n-1, -l6v0-3, -l6ha};
Plane Surface(b) = {b};

// blocks coil 6
lv = l6v0+4 ;
lhu = l6hu ;
lhm = l6hm ;
lhl = l6hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor

// air block between coil 5 and 6
b += 1 ;
Curve Loop(b) = {l5v0, l5hu:l5hu+n-1, -l5v0-1, -l6hl-n+1:-l6hl};
Plane Surface(b) = {b};

// blocks coil 5
lv = l5v0+2 ;
lhu = l5hu ;
lhm = l5hm ;
lhl = l5hl ;
For k In {1:n}
	b += 1 ;
	Curve Loop(b) = {lv, lhm, -lv-2, -lhu};
	Plane Surface(b) = {b};
	
	b += 1 ;
	Curve Loop(b) = {lv+1, lhl, -lv-3, -lhm};
	Plane Surface(b) = {b};
	
	lv += 2 ;
	lhu += 1 ;
	lhm += 1 ;
	lhl += 1 ;
EndFor


//------------------------------------------
// topology printout
//------------------------------------------
// The deck defines no Physical entities, so BELFEM consumes raw gmsh
// geometry tags and input.conf has to name them literally. Any inserted
// curve therefore renumbers the whole contract. Rather than tracking that
// by hand, the block below writes the current topology straight out of the
// deck's own counters into topology.conf, ready to paste into input.conf.

tapes = n + 1 ;                  // 58 turns per coil
halves = 2 * tapes ;             // two 2.5 mm halves per turn

// coil order as input.conf numbers them ( the deck BUILDS 4,3,2,1,8,7,6,5 )
cv1[] = { l1v1, l2v1, l3v1, l4v1, l5v1, l6v1, l7v1, l8v1 } ;
cx0[] = { coil_x1, coil_x1, coil_x1, coil_x1, coil_x2, coil_x2, coil_x2, coil_x2 } ;
cy0[] = { coil_y1, coil_y2, coil_y3, coil_y4, coil_y5, coil_y6, coil_y7, coil_y8 } ;

Printf("// generated by gantry.geo -- r_fillet = %g mm", r_fillet) > "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    iron") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
Printf("        blocks   : 1 ;") >> "topology.conf" ;
Printf("        material : iron ;") >> "topology.conf" ;
Printf("    }") >> "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    air") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
Printf("        blocks : 2:%g ;", b) >> "topology.conf" ;
Printf("    }") >> "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    thinshell : tape") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
Printf("        sidesets : %g:%g, %g:%g, %g:%g, %g:%g, %g:%g, %g:%g, %g:%g, %g:%g ;",
       cv1[0], cv1[0]+halves-1, cv1[1], cv1[1]+halves-1,
       cv1[2], cv1[2]+halves-1, cv1[3], cv1[3]+halves-1,
       cv1[4], cv1[4]+halves-1, cv1[5], cv1[5]+halves-1,
       cv1[6], cv1[6]+halves-1, cv1[7], cv1[7]+halves-1) >> "topology.conf" ;
Printf("    }") >> "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    // y = 0 midplane, current preserved under the mirror -> B x n = 0") >> "topology.conf" ;
Printf("    air symmetry") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
Printf("        sidesets : 3:13, 15 ;") >> "topology.conf" ;
Printf("    }") >> "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    ferro symmetry") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
Printf("        sidesets : 14 ;") >> "topology.conf" ;
Printf("    }") >> "topology.conf" ;
Printf("") >> "topology.conf" ;
Printf("    // one curve per tape: the two 2.5 mm halves of its 5 mm cross section") >> "topology.conf" ;
Printf("    curves") >> "topology.conf" ;
Printf("    {") >> "topology.conf" ;
kk = 0 ;
For c In {0:7}
    For t In {0:tapes-1}
        kk += 1 ;
        Printf("        %4g : %4g, %4g ;   // coil %g, turn %2g  ( x = %.4f   y = %.1f )",
               kk, cv1[c]+2*t, cv1[c]+2*t+1, c+1, t+1,
               cx0[c] + t*d, cy0[c]) >> "topology.conf" ;
    EndFor
EndFor
Printf("    }") >> "topology.conf" ;
