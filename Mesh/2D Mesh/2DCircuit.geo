// Julian Castrillon 
// 2D Mesh

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;

h = 0.003;

Point(1)  = {0, 0, 0, h};
Point(2)  = {0.45, 0, 0, h};
Point(3)  = {0.45, 0.1, 0, h};
Point(4)  = {0, 0.1, 0, h};

Point(5)    = {0.15, 0.08, 0, h};
Point(51)   = {0.14, 0.08, 0, h};
Point(52)   = {0.16, 0.08, 0, h};
Circle(55)  = {51, 5, 52};
Circle(550) = {52, 5, 51};

Point(6)    = {0.15, 0.05, 0, h};
Point(61)   = {0.14, 0.05, 0, h};
Point(62)   = {0.16, 0.05, 0, h};
Circle(66)  = {61, 6, 62};
Circle(660) = {62, 6, 61};

Point(7)    = {0.15, 0.02, 0, h};
Point(71)   = {0.14, 0.02, 0, h};
Point(72)   = {0.16, 0.02, 0, h};
Circle(77)  = {71, 7, 72};
Circle(770) = {72, 7, 71};

Point(8)    = {0.25, 0.065, 0, h};
Point(81)   = {0.24, 0.065, 0, h};
Point(82)   = {0.26, 0.065, 0, h};
Circle(88)  = {81, 8, 82};
Circle(880) = {82, 8, 81};

Point(9)    = {0.25, 0.035, 0, h};
Point(91)   = {0.24, 0.035, 0, h};
Point(92)   = {0.26, 0.035, 0, h};
Circle(99)  = {91, 9, 92};
Circle(990) = {92, 9, 91};

Point(10) = {0.35, 0.08, 0, h};
Point(11) = {0.35, 0.05, 0, h};
Point(12) = {0.35, 0.02, 0, h};

Point(10)     = {0.35, 0.08, 0, h};
Point(101)    = {0.34, 0.08, 0, h};
Point(102)    = {0.36, 0.08, 0, h};
Circle(1010)  = {101, 10, 102};
Circle(10100) = {102, 10, 101};

Point(11)     = {0.35, 0.05, 0, h};
Point(111)    = {0.34, 0.05, 0, h};
Point(112)    = {0.36, 0.05, 0, h};
Circle(1111)  = {111, 11, 112};
Circle(11110) = {112, 11, 111};

Point(12)     = {0.35, 0.02, 0, h};
Point(121)    = {0.34, 0.02, 0, h};
Point(122)    = {0.36, 0.02, 0, h};
Circle(1212)  = {121, 12, 122};
Circle(12120) = {122, 12, 121};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Curve Loop(1) = {1, 2, 3, 4};
Curve Loop(2) = {550, 55};
Curve Loop(3) = {660, 66};
Curve Loop(4) = {770, 77};
Curve Loop(5) = {880, 88};
Curve Loop(6) = {990, 99};
Curve Loop(7) = {10100, 1010};
Curve Loop(8) = {11110, 1111};
Curve Loop(9) = {12120, 1212};
Plane Surface(1) = {1, 2, 3, 4, 5, 6, 7, 8, 9};

Physical Curve("Inlet", 1) = {4};
Physical Curve("BottomWallEntry", 2) = {1};
Physical Curve("Outlet", 3) = {2};
Physical Curve("TopWall", 4) = {3};

Physical Curve("Circle1T", 5) = {55};
Physical Curve("Circle1B", 6) = {550};

Physical Curve("Circle2T", 7) = {66};
Physical Curve("Circle2B", 8) = {660};

Physical Curve("Circle3T", 9) = {77};
Physical Curve("Circle3B", 10) = {770};

Physical Curve("Circle4T", 11) = {88};
Physical Curve("Circle4B", 12) = {880};

Physical Curve("Circle5T", 13) = {99};
Physical Curve("Circle5B", 14) = {990};

Physical Curve("Circle6T", 15) = {1010};
Physical Curve("Circle6B", 16) = {10100};

Physical Curve("Circle7T", 17) = {1111};
Physical Curve("Circle7B", 18) = {11110};

Physical Curve("Circle8T", 19) = {1212};
Physical Curve("Circle8B", 20) = {12120};

Physical Surface("Surface", 1) = {1};

Field[1] = BoundaryLayer;
Field[1].CurvesList = {1, 3, 55, 550, 66, 660, 77, 770, 88, 880, 99, 990, 1010, 10100, 1111, 11110, 1212, 12120};  // Apply to edges 1, 2, and 3
Field[1].hwall_n = 0.0001;           // Sets the first layer thickness to 0.01 units
Field[1].ratio = 1.3;                // Sets a growth rate of 10% per layer
Field[1].thickness = 0.002;          // Sets the total boundary layer thickness to 0.1 units
Field[1].Quads = 0;                  // Use quadrilateral elements
Field[1].IntersectMetrics = 1;       // Smoothly transition into the rest of the mesh
Field[1].PointsList = {1, 2, 3, 4};
BoundaryLayer Field = 1;

Field[2] = Box;
Field[2].Thickness = 0.1;
Field[2].VIn = h/3;
Field[2].VOut = h;
Field[2].XMax = 0.4;
Field[2].XMin = 0.13;
Field[2].YMax = 0.1;
Background Field = 2;

Mesh 2;
Save "2DCircuit.msh";






