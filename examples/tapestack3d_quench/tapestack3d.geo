numTapes = 8 ;

tapeWidth = 4 ;
tapeDistance = 0.1 ;

tapeResolutionTip = 0.1;
tapeResolutionCenter = 0.5 ;
domainResolution = 5 ;

domainRadius = 50 ;
domainLength = 10 ;

y = -(numTapes-1)*0.5*tapeDistance ;

p = 0 ;

// domain
Point(p+1) = { 0,  0,  0,  tapeResolutionCenter };
Point(p+2) = { 0,  domainRadius,  0,  domainResolution };
Point(p+3) = { domainRadius,  0,  0,  domainResolution };
Point(p+4) = { 0,  -domainRadius,  0,  domainResolution };
Point(p+5) = { -domainRadius,  0,  0,  domainResolution };

Point(p+6) = { 0,  0,  domainLength,  tapeResolutionCenter };
Point(p+7) = { 0,  domainRadius,  domainLength,  domainResolution };
Point(p+8) = { domainRadius,  0,  domainLength,  domainResolution };
Point(p+9) = { 0,  -domainRadius,  domainLength,  domainResolution };
Point(p+10) = { -domainRadius,  0,  domainLength,  domainResolution };
p0 = 10 ;

l = 0 ;
Circle(l+1) = {p+2, p+1, p+3};
Circle(l+2) = {p+3, p+1, p+4};
Circle(l+3) = {p+4, p+1, p+5};
Circle(l+4) = {p+5, p+1, p+2};
Circle(l+5) = {p+7, p+6, p+8};
Circle(l+6) = {p+8, p+6, p+9};
Circle(l+7) = {p+9, p+6, p+10};
Circle(l+8) = {p+10, p+6, p+7};
l = l+8 ;

Line(l+1) = {p+2, p+7};
Line(l+2) = {p+3, p+8};
Line(l+3) = {p+4, p+9};
Line(l+4) = {p+5, p+10};
l = l+4 ;

l0 = l;

p=p0 ;

For t In {1:numTapes}

	Point(p+1) = { -0.5*tapeWidth,  y,  0,  tapeResolutionTip };
	Point(p+2) = { 0,  y,  0,  tapeResolutionCenter };
	Point(p+3) = { 0.5*tapeWidth,  y,  0,  tapeResolutionTip };
	Point(p+4) = { -0.5*tapeWidth,  y,  domainLength,  tapeResolutionTip };
	Point(p+5) = { 0,  y,  domainLength,  tapeResolutionCenter };
	Point(p+6) = { 0.5*tapeWidth,  y,  domainLength,  tapeResolutionTip };
	
	p = p + 6 ;
	y = y + tapeDistance ;
EndFor


p=p0 ;
For t In {1:numTapes}

	Line(l+1) = {p+1, p+2};
	Line(l+2) = {p+2, p+3};
	Line(l+3) = {p+6, p+5};
	Line(l+4) = {p+5, p+4};
	Line(l+5) = {p+1, p+4};
	Line(l+6) = {p+2, p+5};
	Line(l+7) = {p+3, p+6};
	
	p = p + 6 ;	
	l = l + 7 ;
EndFor

l1 = l ;
p=p0 ;

For t In {1:numTapes-1}
	l = l + 1 ;
	Line(l) = {p+1, p+7};
	p = p + 6 ;
EndFor

l2 = l ;
p=p0 ;
For t In {1:numTapes-1}
	l = l + 1 ;
	Line(l) = {p+3, p+9};
	p = p + 6 ;
EndFor

l3 = l ;
p=p0 ;
For t In {1:numTapes-1}
	l = l + 1 ;
	Line(l) = {p+4, p+10};
	p = p + 6 ;
EndFor

l4 = l ;
p=p0 ;
For t In {1:numTapes-1}
	l = l + 1 ;
	Line(l) = {p+6, p+12};
	p = p + 6 ;
EndFor
l5 = l ;


// Domain boundaries
Curve Loop(1) = {8, -9, -4, 12};
Surface(1) = {1};
Curve Loop(2) = {7, -12, -3, 11};
Surface(2) = {2};
Curve Loop(3) = {6, -11, -2, 10};
Surface(3) = {3};
Curve Loop(4) = {5, -10, -1, 9};
Surface(4) = {4};

// tapes
s = 4 ;
l = l0 + 1 ;

s0 = s ;
For t In {1:numTapes}
	// left side
	s = s + 1 ;
	Curve Loop(s) = {l, l+5, l+3, -l-4};
	Plane Surface(s) = {s};
	
	// right side
	s = s + 1 ;
	Curve Loop(s) = {l+1, l+6, l+2, -l-5};
	Plane Surface(s) = {s};
	
	l = l + 7 ;
EndFor

l = l0 ;
m = l1 ;

s1 = s ;

a = l0 ;
b = l1 ;
c = l2 ;
d = l3 ;
e = l4 ;

For t In {1:numTapes-1}
	s = s + 1 ;
	Curve Loop(s) = {b+1, a+12, -d-1, -a-5};
	Plane Surface(s) = {s};

	a = a + 7 ;
	b = b + 1 ;
	c = c + 1 ;
	d = d + 1 ;
	e = e + 1 ;	
EndFor

s2 = s ;

a = l0 ;
b = l1 ;
c = l2 ;
d = l3 ;
e = l4 ;

For t In {1:numTapes-1}
	s = s + 1 ;
	Curve Loop(s) = {a+4, d+1, -a-11, -a-10, -e-1, a+3};
	Plane Surface(s) = {s};
	a = a + 7 ;
	b = b + 1 ;
	c = c + 1 ;
	d = d + 1 ;
	e = e + 1 ;	
EndFor
s3 = s ;

a = l0 ;
b = l1 ;
c = l2 ;
d = l3 ;
e = l4 ;

For t In {1:numTapes-1}

	s = s + 1 ;
	Curve Loop(s) = {e+1, -a-14, -c-1, a+7};
	Plane Surface(s) = {s};
	
	a = a + 7 ;
	b = b + 1 ;
	c = c + 1 ;
	d = d + 1 ;
	e = e + 1 ;	
EndFor

s4 = s ;

a = l0 ;
b = l1 ;
c = l2 ;
d = l3 ;
e = l4 ;

For t In {1:numTapes-1}
	s = s + 1 ;
	Curve Loop(s) = {a+2, c+1, -a-9, -a-8, -b-1, a+1};
	Plane Surface(s) = {s};
	
	a = a + 7 ;
	b = b + 1 ;
	c = c + 1 ;
	d = d + 1 ;
	e = e + 1 ;	
EndFor

s5 = s ;

a = s0 + 1 ;
b = s1 + 1 ;
c = s2 + 1 ;
d = s3 + 1 ;
e = s4 + 1 ;

// volume counter
v = 1 ;
For t In {1:numTapes-1}
	v = v + 1 ;
	
	Surface Loop(v) = {a+1, e, d, c, a, b, a+3, a+2};
	Volume(v) = {v};

	a = a + 2 ;
	b = b + 1 ;
	c = c + 1 ;
	d = d + 1 ;
	e = e + 1 ;
	
EndFor

// loop counter
t = s ;

// boundary at front
t = t + 1 ;
Curve Loop(t) = {1, 2, 3, 4};

t = t + 1 ;
Curve Loop(t) = {l0+1, l0+2, l2+1:l3, -l1+5, -l1+6,-l2:-l1-1};

s = s + 1 ;
Plane Surface(s) = {t-1, t};


// boundary at back
t = t + 1 ;
Curve Loop(t) = {5, 6, 7, 8};

t = t + 1 ;
Curve Loop(t) = {l0+3, l0+4, l3+1:l4, -l1+3, -l1+4, -l5:-l4-1};

s = s + 1 ;
Plane Surface(s) = {t-1, t};


// The air region wraps the whole stack, so its boundary needs BOTH outer tape
// faces: tape 1 at the bottom ( s0+1, s0+2 ) and tape numTapes at the top
// ( s1-1, s1 ).  This read "s0-1, s0" = surfaces 3 and 4, i.e. a second copy of
// two cylinder quadrants already listed four terms earlier, while the top tape
// appeared nowhere -- so the air was bounded below by tape 1 and, at the top,
// by nothing.  The two literals "6, 5" are tape 1 and are now written through
// the counters like every other entry, so the loop no longer depends on s0
// happening to be 4.
//   1..4          outer cylinder
//   s-1, s        air end faces at z = 0 and z = domainLength
//   s0+1, s0+2    tape 1        ( bottom of the stack )
//   s1-1, s1      tape numTapes ( top of the stack )
//   s3+1..s4      tip walls at x = +tapeWidth/2
//   s1+1..s2      tip walls at x = -tapeWidth/2
Surface Loop(1) = {1, 2, 3, 4, s, s-1, s1-1, s1, s4:s3+1, s0+2, s0+1, s1+1:s2 };
Volume(1) = {1};

// --- mesh periodicity ---------------------------------------------------
// The solver pairs the z = 0 and z = domainLength faces facet by facet, so
// gmsh must mesh the back faces as translated COPIES of the front faces.
// Without these constraints both planes are meshed independently and the
// facet counts do not even agree.
//   front (z = 0)            back (z = domainLength)
//   s4+1 .. s5   solder      s2+1 .. s3   solder
//   s-1          air         s             air

For k In {1:numTapes-1}
	Periodic Surface { s2+k } = { s4+k } Translate { 0, 0, domainLength };
EndFor

Periodic Surface { s } = { s-1 } Translate { 0, 0, domainLength };

// ========================================================================
//  input.conf snippet
// ========================================================================
//  Everything below is derived from the counters above, so it stays correct
//  when numTapes, tapeWidth or domainLength change.  Run
//
//      gmsh -0 tapestack3d.geo
//
//  and copy the two printed blocks into input.conf.  Only the ids are
//  printed -- the layer stack, the material names and the source amplitude
//  are user choices that the geometry cannot know.
//
//  Surface groups produced above:
//     s0+1 .. s1   the 2*numTapes tape half surfaces (thin shells)
//     s1+1 .. s2   tip walls at x = -tapeWidth/2
//     s2+1 .. s3   solder end faces at z = domainLength
//     s3+1 .. s4   tip walls at x = +tapeWidth/2
//     s4+1 .. s5   solder end faces at z = 0
//     s5+1 / s5+2  air end faces at z = 0 / z = domainLength
//
//  Point ids double as node ids in the mesh, because every geometric point
//  is meshed and gmsh writes the point entities first.

Printf("") ;
Printf("topology") ;
Printf("{") ;
Printf("    thinshell : tape") ;
Printf("    {") ;
Printf("        sidesets : %g:%g ;", s0+1, s1 ) ;
Printf("    }") ;
Printf("") ;
Printf("    air") ;
Printf("    {") ;
Printf("        block : 1 ;") ;
Printf("    }") ;
Printf("    conductor") ;
Printf("    {") ;
Printf("        blocks : 2:%g ;", numTapes ) ;
Printf("    }") ;
Printf("") ;

//  A terminal curve is a tape sideset intersected with the end face that
//  actually carries that tape's edge.  Only tape 1 borders the air end
//  faces -- for every other tape the end face belongs to the solder volume
//  below it.  Intersecting an inner tape with the air face yields a single
//  corner node and aborts in CurveFactory::sort_end_nodes.

Printf("    curves") ;
Printf("    {") ;
For k In {1:numTapes}
	kLeft  = s0 + 2*k - 1 ;
	kRight = s0 + 2*k ;
	If ( k == 1 )
		kFront = s5 + 1 ;
	Else
		kFront = s4 + k - 1 ;
	EndIf
	Printf("        %2g : %2g @ %2g ;   //  tape %g, left",  2*k-1, kLeft,  kFront, k ) ;
	Printf("        %2g : %2g @ %2g ;   //  tape %g, right", 2*k,   kRight, kFront, k ) ;
EndFor
For k In {1:numTapes}
	kLeft  = s0 + 2*k - 1 ;
	kRight = s0 + 2*k ;
	If ( k == 1 )
		kBack = s5 + 2 ;
	Else
		kBack = s2 + k - 1 ;
	EndIf
	Printf("        %2g : %2g @ %2g ;   //  tape %g, left",  2*numTapes+2*k-1, kLeft,  kBack, k ) ;
	Printf("        %2g : %2g @ %2g ;   //  tape %g, right", 2*numTapes+2*k,   kRight, kBack, k ) ;
EndFor
Printf("    }") ;
Printf("") ;

//  Three points spanning each end plane.  They must not be collinear, so
//  both tips of tape 1 plus the left tip of the last tape are used.

Printf("    periodic") ;
Printf("    {") ;
Printf("        source : %g, %g, %g ;", p0+1, p0+3, p0+6*(numTapes-1)+1 ) ;
Printf("        target : %g, %g, %g ;", p0+4, p0+6, p0+6*(numTapes-1)+4 ) ;
Printf("    }") ;
Printf("}") ;
Printf("") ;

//  The current terminal.  One bracket = one condition, and the number of
//  conditions may not exceed the number of cohomology generators of the air
//  region minus the free axial generator contributed by the periodicity.
//  Which form is correct depends on whether the volumes between the tapes
//  conduct:
//
//    solder between the tapes ( blocks 2..numTapes are conductor )
//        the stack is ONE conductor, so there is exactly one condition and
//        the cycle has to encircle the whole cross-section.  Give the solder
//        end faces as SIDESETS -- the bulk branch of suggest_Homology takes
//        their boundary, which is the envelope of the stack -- and set the
//        amplitude to the TOTAL current.  This is what is printed below.
//
//    insulation between the tapes ( those blocks are air )
//        every tape is its own conductor, so use the per-tape curves with
//        one bracket per tape and the per-tape current:
//            input curves  : [1,2],[3,4],...
//            output curves : [2n+1,2n+2],...
//
//  Do not use the curve form for the soldered stack: its generator is the
//  plain sum of the listed curves, which misses the solder current, and the
//  inner tape curves lie inside the conductor.  Asking for more conditions
//  than the topology allows aborts in
//  Cohomology::updatekGeneratorsFromHomology.

Printf("boundary conditions") ;
Printf("{") ;
Printf("    current") ;
Printf("    {") ;
Printf("        input terminals  : [%g:%g] ;   //  solder end faces at z = 0", s4+1, s5 ) ;
Printf("        output terminals : [%g:%g] ;   //  solder end faces at z = domainLength", s2+1, s3 ) ;
Printf("    }") ;
Printf("}") ;
Printf("") ;
