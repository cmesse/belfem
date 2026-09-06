numTapes = 11 ;

tapeWidth = 4 ;
tapeDistance = 0.2 ;

tapeResolutionTip = 0.01  ;
tapeResolutionCenter = 0.1 ;
domainResolution = 1 ;

domainRadius = 25 ;

y = -(numTapes-1)*0.5*tapeDistance ;

p = 0 ;


Point(p+1) = { 0,  0,  0,  tapeResolutionCenter };
Point(p+2) = { 0,  domainRadius,  0,  domainResolution };
Point(p+3) = { domainRadius,  0,  0,  domainResolution };
Point(p+4) = { 0,  -domainRadius,  0,  domainResolution };
Point(p+5) = { -domainRadius,  0,  0,  domainResolution };

l = 0 ;
Circle(l+1) = {p+2, p+1, p+3};
Circle(l+2) = {p+3, p+1, p+4};
Circle(l+3) = {p+4, p+1, p+5};
Circle(l+4) = {p+5, p+1, p+2};
l = l+4 ;


p=5 ;

For t In {1:numTapes}

	Point(p+1) = { -0.5*tapeWidth,  y,  0,  tapeResolutionTip };
	Point(p+2) = { 0,  y,  0,  tapeResolutionCenter };
	Point(p+3) = { 0.5*tapeWidth,  y,  0,  tapeResolutionTip };
	
	p = p + 3 ;	
	y = y + tapeDistance ;
EndFor


p=5 ;
For t In {1:numTapes}

	Line(l+1) = {p+1, p+2};
	Line(l+2) = {p+2, p+3};

	
	p = p + 3 ;	
	l = l + 2 ;
EndFor

a = l ;
p = 5 ;
For t In {1:numTapes-1}
	Line(l+1) = {p+1, p+4};
	l = l + 1 ;
	p = p + 3 ;
EndFor
p = p - 3 ;
b = l ;
For t In {1:numTapes-1}
	Line(l+1) = {p+6, p+3};
	l = l + 1 ;
	p = p - 3 ;
EndFor
d = l ;

p = 3*numTapes + 5 ;


c = 0 ;

l=4 ;
m = a + 1;
n = b+numTapes-1;

For t In {1:numTapes-1}
	c = c + 1 ;
    Curve Loop(c) = {l+1,l+2,-n,-l-4,-l-3,-m};
	l = l + 2 ;
	m = m + 1 ;
	n = n - 1 ;
	Plane Surface(c) = {c};
EndFor

l = 4;
p = 0 ;

m = a+numTapes-1;
n = b+numTapes-1;
Printf( "a %g", a );
Printf( "b %g", b );
Printf( "m %g", m );
Printf( "n %g", n );
Printf( "c %g", c );

Curve Loop(c+1) = {1,2,3,4};
Curve Loop(c+2) = {a+1:m, a-1, a, b+1:n, -l-2, -l-1};
Plane Surface(c+1) = {c+1,c+2};
