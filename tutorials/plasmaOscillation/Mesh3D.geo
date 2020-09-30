//Gmsh mesh file https://gmsh.info
L1=6;
nCells1=120;
L2=0.8;
nCells2=16;
L3=0.8;
nCells3=16;

Point(1) = {0,L2/2,L3/2};
Point(2) = {0,-L2/2,L3/2};
Point(3) = {0,-L2/2,-L3/2};
Point(4) = {0,L2/2,-L3/2};

Point(5) = {L1,L2/2,L3/2};
Point(6) = {L1,-L2/2,L3/2};
Point(7) = {L1,-L2/2,-L3/2};
Point(8) = {L1,L2/2,-L3/2};


//y-dir
Line(1) = {1, 2};
Line(2) = {4, 3};
Line(3) = {5, 6};
Line(4) = {8, 7};

//z-dir
Line(5) = {4, 1};
Line(6) = {3, 2};
Line(7) = {8, 5};
Line(8) = {7, 6};

//x-dir
Line(9) = {1, 5};
Line(10) = {2, 6};
Line(11) = {3, 7};
Line(12) = {4, 8};


//X
Transfinite Line {9:12} = nCells1+1 Using Progression 1.0;
//Y
Transfinite Line {1:4} = nCells2+1 Using Progression 1.0;
//Z
Transfinite Line {5:8} = nCells3+1 Using Progression 1.0;


//Kathode
Line Loop(1) = {5, 1, -6, -2};
Surface(1) = {1};

//Anode
Line Loop(2) = {3, -8, -4, 7};
Surface(2) = {2};

//Empty Y-Richtung
Line Loop(3) = {1, 10, -3, -9};
Surface(3) = {3};
Line Loop(4) = {4, -11, -2, 12};
Surface(4) = {4};

//Empty Z-Richtung
Line Loop(5) = {5, 9, -7, -12};
Surface(5) = {5};
Line Loop(6) = {8, -10, -6, 11};
Surface(6) = {6};

For j In {1:6}
Transfinite Surface {j} = {};
EndFor
For j In {1:6}
Recombine Surface {j};
EndFor


//volume
Surface Loop(1) = {3, 1, 5, 2, 6, 4};
Volume(1) = {1};
Transfinite Volume {1} = {};



Physical Surface("leftBoundary") = {1};
Physical Surface("rightBoundary") = {2};
Physical Surface("upperBoundary") = {5};
Physical Surface("lowerBoundary") = {6};
Physical Surface("frontBoundary") = {3};
Physical Surface("backBoundary") = {4};
Physical Volume("internalField") = {1};


