//Gmsh mesh file https://gmsh.info
nCells = 200;
L1=0.05;
L2=0.024;


Point(1) = {-L1/2,L2/2,L2/2};
Point(2) = {-L1/2,-L2/2,L2/2};
Point(3) = {-L1/2,-L2/2,-L2/2};
Point(4) = {-L1/2,L2/2,-L2/2};

Point(5) = {L1/2,L2/2,L2/2};
Point(6) = {L1/2,-L2/2,L2/2};
Point(7) = {L1/2,-L2/2,-L2/2};
Point(8) = {L1/2,L2/2,-L2/2};

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
Transfinite Line {9:12} = nCells+1 Using Progression 1.0;
//Y+Z
Transfinite Line {1:8} = 1 Using Progression 1.0;


//left boundary
Line Loop(1) = {5, 1, -6, -2};
Surface(1) = {1};

//right boundary
Line Loop(2) = {3, -8, -4, 7};
Surface(2) = {2};

//empty y-dir
Line Loop(3) = {1, 10, -3, -9};
Surface(3) = {3};
Line Loop(4) = {4, -11, -2, 12};
Surface(4) = {4};

//empty z-dir
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

//volumen
Surface Loop(1) = {3, 1, 5, 2, 6, 4};
Volume(1) = {1};
Transfinite Volume {1} = {};


Physical Surface("source") = {1};
Physical Surface("collector") = {2};
Physical Surface("emptyBoundary") = {3,4,5,6};
Physical Volume("internalField") = {1};


