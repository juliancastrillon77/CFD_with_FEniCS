// Julian Castrillon 
// 2D Mesh

h = 0.01;

Point(1)  = {0, 0, 0, h};
Point(2)  = {5, 0, 0, h};
Point(3)  = {5, 1, 0, h};
Point(4)  = {0, 1, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

//+
Curve Loop(1) = {4, 1, 2, 3};
Plane Surface(1) = {1};

//+
Physical Curve("Entry", 8)          = {4};
Physical Curve("BottomWall", 9)     = {1};
Physical Curve("Exit", 10)          = {2};
Physical Curve("TopWall", 11)       = {3};
Physical Surface("Surface", 2) = {1};

