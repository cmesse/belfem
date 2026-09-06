numTapes = 2 ;

tapeWidth = 4 ;
tapeDistance = 1 ;

phi = 0.5*Sqrt( 5 ) - 0.5 ;
tapeResolution       = 0.3;
tapeResolutionCenter = tapeResolution * phi  * phi;
tapeResolutionQuench = tapeResolutionCenter * phi * phi;
domainResolution = 5 ;
blastWidth   = 1.8 ;
blastLength  = 1 + phi  ;

domainRadius = 20 ;
domainLength = 12.5 ;


// quench area
Point(1) = { 0,  0,  0,  tapeResolutionQuench };
Point(2) = { 0,  blastWidth,  0,  tapeResolutionQuench };
Point(3) = { blastLength,  0,  0,  tapeResolutionQuench };
Point(4) = { 0,  -blastWidth,  0,  tapeResolutionQuench };
Point(5) = {  -blastLength,  0, 0, tapeResolutionQuench };


Point(6) = { 0.5*domainLength,  -0.5*tapeWidth,  0,  tapeResolution };

Point(7) = {  0,  -0.5*tapeWidth,  0,  tapeResolutionCenter };
Point(8) = { 0.5*domainLength,   0.5*tapeWidth,  0,  tapeResolution };
Point(9) = {-0.5*domainLength,   0.5*tapeWidth,  0,  tapeResolution };
Point(10) = {  0,  0.5*tapeWidth,  0,  tapeResolutionCenter };
Point(11) = {-0.5*domainLength,  -0.5*tapeWidth,  0,  tapeResolution };


Ellipse(1) = {3, 1, 2, 2};
Ellipse(2) = {2, 1, 5, 5};
Ellipse(3) = {5, 1, 4, 4};
Ellipse(4) = {4, 1, 3, 3};
Line(5) = {11, 7};
Line(6) = {7, 6};
Line(7) = {8, 8};
Line(8) = {6, 8};
Line(9) = {8, 10};
Line(10) = {10, 9};
Line(11) = {9, 11};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};
Curve Loop(2) = {6, 8, 9, 10, 11, 5};
Plane Surface(2) = {1, 2};

Point(12) = { -0.5*domainLength,  0,  0,  domainResolution };
Point(13) = { -0.5*domainLength,  domainRadius,  0,  domainResolution };
Point(14) = { -0.5*domainLength,  0,  domainRadius,  domainResolution };
Point(15) = { -0.5*domainLength,  -domainRadius,  0,  domainResolution };
Point(16) = { -0.5*domainLength,  0,  -domainRadius,  domainResolution };


Point(17) = { 0.5*domainLength,  0,  0,  domainResolution };
Point(18) = { 0.5*domainLength,  domainRadius,  0,  domainResolution };
Point(19) = { 0.5*domainLength,  0,  domainRadius,  domainResolution };
Point(20) = { 0.5*domainLength,  -domainRadius,  0,  domainResolution };
Point(21) = { 0.5*domainLength,  0,  -domainRadius,  domainResolution };
Circle(12) = {16, 12, 15};
Circle(13) = {15, 12, 14};
Circle(14) = {14, 12, 13};
Circle(15) = {13, 12, 16};
Circle(16) = {20, 17, 19};
Circle(17) = {19, 17, 18};
Circle(18) = {18, 17, 21};
Circle(19) = {21, 17, 20};
Line(20) = {16, 21};
Line(21) = {15, 20};
Line(22) = {14, 19};
Line(23) = {13, 18};
Line(24) = {13, 9};
Line(25) = {11, 15};
Line(26) = {20, 6};
Line(27) = {8, 18};
Curve Loop(3) = {14, 24, 11, 25, 13};
Plane Surface(3) = {3};
Curve Loop(4) = {24, 11, 25, -12, -15};
Plane Surface(4) = {4};
Curve Loop(5) = {26, 8, 27, 18, 19};
Plane Surface(5) = {5};
Curve Loop(6) = {17, -27, -8, -26, 16};
Plane Surface(6) = {6};
Curve Loop(7) = {15, 20, -18, -23};
Surface(7) = {7};
Curve Loop(8) = {12, 21, -19, -20};
Surface(8) = {8};
Curve Loop(9) = {13, 22, -16, -21};
Surface(9) = {9};
Curve Loop(10) = {14, 23, -17, -22};
Surface(10) = {10};
Curve Loop(11) = {24, -10, -9, 27, -23};
Plane Surface(11) = {11};
Curve Loop(12) = {25, 21, 26, -6, -5};
Plane Surface(12) = {12};
Surface Loop(1) = {4, 8, 5, 7, 11, 2, 1, 12};
Volume(1) = {1};
Surface Loop(2) = {6, 10, 3, 9, 11, 2, 1, 12};
Volume(2) = {2};
