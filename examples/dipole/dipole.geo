
coilDiameter = 1.0 ;
coilDistance = 4.0 ;
domainRadius = 10 ;

domainResolution = 0.5 ;
coilResolution = 0.05 ;


Point(1) = {0, 0, 0, coilResolution } ;

// Left Coil
Point(2) = {-0.5*coilDistance, 0, 0, coilResolution } ;
Point(3) = {-0.5*coilDistance+0.5*coilDiameter, 0, 0, coilResolution } ;
Point(4) = {-0.5*coilDistance,0.5*coilDiameter, 0, coilResolution } ;
Point(5) = {-0.5*coilDistance-0.5*coilDiameter, 0, 0, coilResolution } ;
Point(6) = {-0.5*coilDistance,-0.5*coilDiameter, 0, coilResolution } ;
Circle(1) = {3, 2, 4};
Circle(2) = {4, 2, 5};
Circle(3) = {5, 2, 6};
Circle(4) = {6, 2, 3};
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

// Right Coil
Point(7) = {0.5*coilDistance, 0, 0, coilResolution } ;
Point(8) = {0.5*coilDistance+0.5*coilDiameter, 0, 0, coilResolution } ;
Point(9) = {0.5*coilDistance,0.5*coilDiameter, 0, coilResolution } ;
Point(10) = {0.5*coilDistance-0.5*coilDiameter, 0, 0, coilResolution } ;
Point(11) = {0.5*coilDistance,-0.5*coilDiameter, 0, coilResolution } ;
Circle(5) = {8, 7, 9};
Circle(6) = {9, 7, 10};
Circle(7) = {10, 7, 11};
Circle(8) = {11, 7, 8};
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

// Domain Edge
Point(12) = { domainRadius, 0, 0, domainResolution } ;
Point(13) = { 0, domainRadius, 0, domainResolution } ;
Point(14) = { -domainRadius, 0, 0, domainResolution } ;
Point(15) = { 0, -domainRadius, 0, domainResolution } ;

Circle(9) = {12, 1, 13};
Circle(10) = {13, 1, 14};
Circle(11) = {14, 1, 15};
Circle(12) = {15, 1, 12};
Curve Loop(3) = {10, 11, 12, 9};
Plane Surface(3) = {3, 1, 2};
