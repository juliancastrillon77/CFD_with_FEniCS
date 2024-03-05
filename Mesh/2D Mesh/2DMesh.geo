// Julian Castrillon 
// 2D Mesh

h = 0.02;

Point(1)  = {0, 0, 0, h};
Point(2)  = {5, 0, 0, h};
Point(3)  = {5, 1, 0, h};
Point(4)  = {0, 1, 0, h};
Point(5)  = {1.5, 0.25, 0, h};

Line(1) = {1, 2};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Point(6)  = {1.35, 0.6, 0, h};
Point(7)  = {1.5, 0.6, 0, h};
Point(8)  = {1.65, 0.6, 0, h};

Circle(5) = {8, 7, 6};
Line(6) = {6, 5};
Line(7) = {5, 8};

//+
Curve Loop(1) = {4, 1, 2, 3};
Curve Loop(2) = {6, 7, 5};
Plane Surface(1) = {1, 2};

//+
Physical Curve("Entry", 8)          = {4};
Physical Curve("BottomWall", 9)     = {1};
Physical Curve("Exit", 10)          = {2};
Physical Curve("TopWall", 11)       = {3};
Physical Curve("Circle", 12)        = {5};
Physical Curve("TriangleLeft", 13)  = {6};
Physical Curve("TriangleRight", 14) = {7};
Physical Surface("Surface", 15) = {1};

//+
Field[1] = Box;
Field[1].Thickness = 0.5;
Field[1].VIn       = 0.004;
Field[1].VOut      = 0.02;
Field[1].XMax      = 1.8;
Field[1].XMin      = 1.3;
Field[1].YMax      = 0.9;
Field[1].YMin      = 0.1;

//+
//Field[2] = Box;
//Field[2].Thickness = 0.5;
//Field[2].VIn       = 0.01;
//Field[2].VOut      = 0.02;
//Field[2].XMax      = 5;
//Field[2].XMin      = 3;
//Field[2].YMax      = 0.9;
//Field[2].YMin      = 0.1;

//+
//Field[3] = Min;
//Field[3].FieldsList = {1, 2};
Background Field = 1; //3;
