lc = 0.02;

Point(1) = {0, 0, 0, lc};
Point(2) = {1, 0, 0, lc};
Point(3) = {1, 1, 0, lc};
Point(4) = {0, 1, 0, lc};

Line(1) = {1,2};
Line(2) = {2,3};
Line(3) = {3,4};
Line(4) = {4,1};

Curve Loop(1) = {1, 2, 3, 4};
Plane Surface(1) = {1};

//Extrude {0,0,1} { Point{1} Layers{10}; }
//Extrude { 0,0.5,0} { Curve Loop{1} }

//+
Field[1] = Distance;
Field[1].CurvesList = {1, 2};
Field[1].Sampling = 90;
Field[1].PointsList = {1, 2, 3};

Field[2] = Threshold;
Field[2].InField = 1;
Field[2].SizeMin = 0.02;
Field[2].SizeMax = lc;
Field[2].DistMin = 0.1;
Field[2].DistMax = 0.11;


Field[7] = Min;
Field[7].FieldsList = {1, 2};
Background Field = 7;


Mesh 2;//+