// Julian Castrillon 
// 2D Mesh

h = 0.02;

Point(1)  = {0, 0, 0, h};
Point(2)  = {5, 0, 0, h};
Point(3)  = {5, 1, 0, h};
Point(4)  = {0, 1, 0, h};
Point(5)  = {0, 0.65, 0, h};
Point(6)  = {0, 0.35, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 5};
Line(5) = {5, 6};
Line(6) = {6, 1};

//+
Curve Loop(1) = {1, 2, 3, 4, 5, 6};
Plane Surface(1) = {1};

//+
Physical Curve("EntryTop", 8)       = {4};
Physical Curve("Nozzle", 12)        = {5};
Physical Curve("EntryBottom", 13)   = {6};
Physical Curve("BottomWall", 9)     = {1};
Physical Curve("Exit", 10)          = {2};
Physical Curve("TopWall", 11)       = {3};

Physical Surface("Surface", 2) = {1};


