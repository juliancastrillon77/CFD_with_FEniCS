// Julian Castrillon 
// 2D Mesh

Mesh.MshFileVersion = 2.2;
Mesh.Binary = 0;

h = 0.02;

Point(1)  = {0, 0, 0, h};
Point(2)  = {0.75, 0, 0, h};
Point(3)  = {0.8, 0, 0, h};
Point(4)  = {0.85, 0, 0, h};
Point(5)  = {1.8, 0, 0, h};
Point(6)  = {1.8, 0.1, 0, h};
Point(7) = {0, 0.1, 0, h};

Line(1)  = {1, 2};
Line(2)  = {4, 5};
Line(3)  = {5, 6};
Line(4)  = {6, 7};
Line(5)  = {7, 1};
Circle(6) = {4, 3, 2};

Physical Curve("Inlet", 1) = {5};
Physical Curve("BottomWall", 2) = {1, -6, 2};
Physical Curve("Outlet", 3) = {3};
Physical Curve("TopWall", 4) = {4};

Curve Loop(1) = {1, 2, -6, 3, 4, 5};
Plane Surface(1) = {1};

Physical Surface("Surface", 1) = {1};

Field[1] = BoundaryLayer;
Field[1].CurvesList = {1, 2, 4, 6};  // Apply to edges 1, 2, and 3
Field[1].hwall_n = 0.0001;            // Sets the first layer thickness to 0.01 units
Field[1].ratio = 1.4;                  // Sets a growth rate of 10% per layer
Field[1].thickness = 0.005;          // Sets the total boundary layer thickness to 0.1 units
Field[1].Quads = 0;                  // Use quadrilateral elements
Field[1].IntersectMetrics = 1;       // Smoothly transition into the rest of the mesh
Field[1].PointsList = {1, 5, 6, 7};
BoundaryLayer Field = 1;

Field[2] = Box;
Field[2].Thickness = 0.1;
Field[2].VIn = h/4;
Field[2].VOut = h;
Field[2].XMax = 1.2;
Field[2].XMin = 0.6;
Field[2].YMax = 0.1;
Background Field = 2;

Mesh.MeshSizeFactor = 1;

Mesh 2;
Save "2DChannelSemiCircle.msh";





