cl1 = 0.0005;
cl2 = 0.001;
Point(1) = {-0.05, -0.01, 0, cl2};
Point(2) = {0, -0.01, 0, cl1};
Point(3) = {0.05, -0.01, 0, cl2};
Point(4) = {0.05, 0.01, 0, cl2};
Point(5) = {0, 0.01, 0, cl1};
Point(6) = {-0.05, 0.01, 0, cl2};
Point(7) = {0.0, 0.0, 0, cl1};
Point(8) = {-0.0, -0.005, 0, cl1};
Point(9) = {-0.005, -0.0, 0, cl1};
Point(10) = {-0.0, 0.005, 0, cl1};
Point(11) = {0.005, 0.0, 0, cl1};
//+
Line(1) = {1, 2};
//+
Line(2) = {2, 3};
//+
Line(3) = {3, 4};
//+
Line(4) = {4, 6};
//+
Line(5) = {6, 1};
//+
Circle(6) = {8, 7, 9};
//+
Circle(7) = {9, 7, 10};
//+
Circle(8) = {10, 7, 11};
//+
Circle(9) = {11, 7, 8};
//+
Curve Loop(1) = {5, 1, 2, 3, 4};
//+
Curve Loop(2) = {7, 8, 9, 6};
//+
Plane Surface(1) = {1, 2};
