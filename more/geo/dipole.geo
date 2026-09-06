
coilDiameter = 1 ;
coilDistance = 6 ;
domainRadius = 50 ;

coilResolution = 0.1 ;
edgeResolution = 3 ;

// Coil 1
Point(1)  = { -0.5*coilDistance, 0, 0, coilResolution};
Point(2)  = { -0.5*coilDistance+0.5*coilDiameter, 0, 0, coilResolution};
Point(3)  = { -0.5*coilDistance, 0.5*coilDiameter, 0, coilResolution};
Point(4)  = { -0.5*coilDistance-0.5*coilDiameter, 0, 0, coilResolution};
Point(5)  = { -0.5*coilDistance,-0.5*coilDiameter, 0, coilResolution};

Circle(1) = {2, 1, 3};
Circle(2) = {3, 1, 4};
Circle(3) = {4, 1, 5};
Circle(4) = {5, 1, 2};


// Coil 2
Point(6)  = {  0.5*coilDistance, 0, 0, coilResolution};
Point(7)  = {  0.5*coilDistance+0.5*coilDiameter, 0, 0, coilResolution};
Point(8)  = {  0.5*coilDistance, 0.5*coilDiameter, 0, coilResolution};
Point(9)  = {  0.5*coilDistance-0.5*coilDiameter, 0, 0, coilResolution};
Point(10) = {  0.5*coilDistance,-0.5*coilDiameter, 0, coilResolution};

Circle(5) = {7, 6, 8};
Circle(6) = {8, 6, 9};
Circle(7) = {9, 6, 10};
Circle(8) = {10, 6, 7};

// Domain
Point(11)  = {  0, 0, 0, coilResolution };
Point(12)  = {  domainRadius, 0, 0, edgeResolution };
Point(13)  = {  0, domainRadius, 0, edgeResolution };
Point(14)  = {  -domainRadius, 0, 0, edgeResolution };
Point(15)  = {  0, -domainRadius, 0, edgeResolution };

Circle(9) = {12, 11, 13};
Circle(10) = {13, 11, 14};
Circle(11) = {14, 11, 15};
Circle(12) = {15, 11, 12};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(2) = {2};

Curve Loop(3) = {9, 10, 11, 12};
Plane Surface(3) = {3, 1, 2};
