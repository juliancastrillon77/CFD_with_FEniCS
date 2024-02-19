// Julian Castrillon
// 3D Mesh

h = 0.3;

//+
Point(1)  = {0, 0, 0, h};
Point(2)  = {5, 0, 0, h};
Point(3)  = {5, 1, 0, h};
Point(4)  = {0, 1, 0, h};
Point(5)  = {0, 0, 1, h};
Point(6)  = {5, 0, 1, h};
Point(7)  = {5, 1, 1, h};
Point(8)  = {0, 1, 1, h};
Point(9)   = {1.5, 0.2, 0.5, h};
Point(10)  = {1.3, 0.6, 0.5, h};
Point(11)  = {1.5, 0.6, 0.5, h};
Point(12)  = {1.5, 0.8, 0.5, h};

//+
Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};
Line(5) = {5, 6};
Line(6) = {6, 7};
Line(7) = {7, 8};
Line(8) = {8, 5};
Line(9) = {1, 5};
Line(10) = {8, 4};
Line(11) = {2, 6};
Line(12) = {7, 3};
Line(13) = {9, 10};
Circle(14) = {10, 11, 12};

//+
Extrude {{0, 1, 0}, {1.5, 0.6, 0.5}, Pi/2} {
  Curve{13}; Curve{14}; 
}
Extrude {{0, 1, 0}, {1.5, 0.6, 0.5}, Pi/2} {
  Curve{15}; Curve{18}; 
}
Extrude {{0, 1, 0}, {1.5, 0.6, 0.5}, Pi/2} {
  Curve{21}; Curve{24}; 
}
Extrude {{0, 1, 0}, {1.5, 0.6, 0.5}, Pi/2} {
  Curve{27}; Curve{30}; 
}

//+
Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(40) = {1}; // Left Wall
Curve Loop(2) = {5, 6, 7, 8};
Plane Surface(41) = {2}; // Right Wall
Curve Loop(3) = {9, -8, 10, 4};
Plane Surface(42) = {3}; // Entry Wall
Curve Loop(4) = {11, 6, 12, -2};
Plane Surface(43) = {4}; // Exit Wall
Curve Loop(5) = {1, 11, -5, -9};
Plane Surface(44) = {5}; // Bottom Wall
Curve Loop(6) = {3, -10, -7, 12};
Plane Surface(45) = {6}; // Top Wall

//+
Surface Loop(1) = {40, 44, 43, 41, 45, 42};
Surface Loop(2) = {20, 26, 32, 38, 35, 29, 23, 17};
Volume(1) = {1, 2};

//+
Physical Surface("LeftWall", 46) = {40};
Physical Surface("RightWall", 47) = {41};
Physical Surface("Entry", 48) = {42};
Physical Surface("Exit", 49) = {43};
Physical Surface("BottomWall", 50) = {44};
Physical Surface("TopWall", 51) = {45};

//+
Physical Surface("SphereForwardRight", 52) = {26};
Physical Surface("SphereForwardLeft", 53) = {32};
Physical Surface("SphereBackwardLeft", 54) = {38};
Physical Surface("SphereBackwardRight", 55) = {20};
Physical Surface("ConeForwardRight", 56) = {23};
Physical Surface("ConeForwardLeft", 57) = {29};
Physical Surface("ConeBackwardLeft", 58) = {35};
Physical Surface("ConeBackwardRight", 59) = {17};
Physical Volume("Vol", 60) = {1};


