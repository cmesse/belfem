numTapes = 2 ;

discDiameter = 76.2 ;
discThickness = 12.2 ;

diskResolution = 4 ;

domainResolution = 25 ;
blastWidth   = 1.8 ;
blastLength  = 1 + phi  ;

domainRadius = 200 ;

// Center
Point(1) = { 0,  0,  0,  diskResolution };

// Disk
Point(2) = { 0, 0 , -0.5*discThickness,  diskResolution };
Point(3) = { 0.5*discDiameter, 0 , -0.5*discThickness,  diskResolution };
Point(4) = { 0, 0.5*discDiameter, -0.5*discThickness,  diskResolution };
Point(5) = { 0, -0.5*discDiameter, -0.5*discThickness,  diskResolution };
Point(6) = { -0.5*discDiameter, 0 , -0.5*discThickness,  diskResolution };

Circle(1) = {3, 2, 5};
Circle(2) = {5, 2, 6};
Circle(3) = {6, 2, 4};
Circle(4) = {4, 2, 3};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Point(7) = { 0, 0 , 0.5*discThickness,  diskResolution };
Point(8) = { 0.5*discDiameter, 0 , 0.5*discThickness,  diskResolution };
Point(9) = { 0, 0.5*discDiameter, 0.5*discThickness,  diskResolution };
Point(10) = { 0, -0.5*discDiameter, 0.5*discThickness,  diskResolution };
Point(11) = { -0.5*discDiameter, 0 , 0.5*discThickness,  diskResolution };
Circle(5) = {8, 7, 10};
Circle(6) = {10, 7, 11};
Circle(7) = {11, 7, 9};
Circle(8) = {9, 7, 8};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Line(9) = {3, 8};
Line(10) = {4, 9};
Line(11) = {6, 11};
Line(12) = {5, 10};
Curve Loop(3) = {10, -7, -11, 3};
Surface(3) = {3};
Curve Loop(4) = {2, 11, -6, -12};
Surface(4) = {4};
Curve Loop(5) = {12, -5, -9, 1};
Surface(5) = {5};
Curve Loop(6) = {4, 9, -8, -10};
Surface(6) = {6};
Surface Loop(1) = {2, 6, 1, 5, 4, 3};
Volume(1) = {1};


// Domain
Point(12) = { domainRadius, 0, 0, domainResolution };
Point(13) = { 0, domainRadius, 0, domainResolution };
Point(14) = { -domainRadius, 0, 0, domainResolution };
Point(15) = { 0, -domainRadius, 0, domainResolution };
Point(16) = { 0, 0, -domainRadius, domainResolution };
Point(17) = { 0, 0, domainRadius, domainResolution };



Circle(13) = {12, 1, 13};
Circle(14) = {13, 1, 14};
Circle(15) = {14, 1, 15};
Circle(16) = {15, 1, 12};
Circle(17) = {12, 1, 17};
Circle(18) = {13, 1, 17};
Circle(19) = {17, 1, 14};
Circle(20) = {15, 1, 17};
Circle(21) = {12, 1, 16};
Circle(22) = {13, 1, 16};
Circle(23) = {14, 1, 16};
Circle(24) = {15, 1, 16};

Curve Loop(7) = {13, 18, -17};
Surface(7) = {7};
Curve Loop(8) = {14, -19, -18};
Surface(8) = {8};
Curve Loop(9) = {15, 20, 19};
Surface(9) = {9};
Curve Loop(10) = {16, 17, -20};
Surface(10) = {10};
Curve Loop(11) = {16, 21, -24};
Surface(11) = {11};
Curve Loop(12) = {15, 24, -23};
Surface(12) = {12};
Curve Loop(13) = {23, -22, 14};
Surface(13) = {13};
Curve Loop(14) = {21, -22, -13};
Surface(14) = {14};
Surface Loop(2) = {9, 12, 11, 10, 7, 14, 13, 8};
Volume(2) = {1, 2};
