numTapes = 2 ;

tapeWidth = 12 ;
tapeDistance = 1 ;

tapeResolutionTip = 1; //0.9
tapeResolutionCenter = 1 ; //1
domainResolution = 5 ;

domainRadius = 40 ;
domainLength = 50 ;

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


Surface Loop(1) = {1, 2, 3, 4, s, s-1, s0-1,s0, s4:s3+1, 6, 5, s1+1:s2 };
Volume(1) = {1};

a = 5;

// periodicity 1
Printf("Tapes :");
For t In {1:numTapes}
    Printf("%g : %g, %g ", t, a, a+1 );
    a = a + 2 ;
EndFor

a = s4 +1;
Printf(" ");


Printf("Boundary");
Printf("0 : 1, 2, 3, 4");
Printf(" ");


a = s4 +1;
Printf(" ");
// periodicity 1
Printf("Periodicity 1 :");
    Printf("0 : %g ", s-1 );
For t In {1:numTapes-1}
    Printf("%g : %g ", t, a );
    a = a + 1 ;
EndFor

Printf(" ");

a = s2 +1;
// periodicity 2
Printf("Periodicity 2 :");
    Printf("0 : %g ", s );
For t In {1:numTapes-1}
    Printf("%g : %g ", t, a );
    a = a + 1 ;
EndFor
