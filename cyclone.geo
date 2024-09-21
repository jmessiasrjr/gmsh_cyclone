//      Setting OpenCASCADE -- Gmsh 3.0.6
SetFactory("OpenCASCADE");

/*  Dimensions of the hydrocyclone
---------------------------------------------------------
    D      -> Diameter between 2 conic section [m]
    Dc     -> Diameter of cylindrical part [m]
---------------------------------------------------------
    Di     -> Diameter of feed pipe [m]
    Do     -> Diameter of overflow pipe [m]
    Du     -> Diameter of underflow pipe [m]
    VF     -> Vortex finder length [m]
    L1     -> length of cylindrical section [mm]
    theta1 -> Angle of first conic section [degree]
    theta2 -> Angle of second conic section [degree]
---------------------------------------------------------*/
//     Constants
D = 0.040;
Dc = 0.070;
//     Variables
Di = 0.0125;
Do = 0.0125;
Du = 0.0175;
VF = 0.022;
L1 = 0.035;
theta1 = 10;
theta2 = 8.5;

//     Admitted values
ho = 0.02; // height of the overflow pipe that is going out of the hydrocyclone
hi = 0.04; // height of the inlet pipe that is going out of the hydrocyclone
mes = 1.0; // mesh element size
Nzb = 30; Nxb = 30; // hydrocyclone body -- number of divisions
Nzo = 30; Nxo = 30; // overflow pipe -- number of divisions
Rzb = 1.2; Rxb = 1.; // hydrocyclone body -- geometric progression interpolation ratio
Rzo = 1.2; Rxo = 1.; // overflow pipe -- geometric progression interpolation ratio
nolb = 50;  // number of layers in hydrocyclone body extrusion
nolo = 50;  // number of layers in overflow pipe extrusion
npi = 20;   // number of points in Spline line 

// =========== POINTS AND SURFACES ===========  

//     Origin point
Point(1) = {0, 0, 0, mes};
//     Cylindrical section
Point(2) = {0, Dc/2, 0, mes};
Point(3) = {0, Dc/2, L1, mes};
Point(4) = {0, 0, L1, mes};

Line(1) = {2, 1};
Line(2) = {2, 3};
Line(3) = {3, 4};
Line(4) = {4, 1};

Transfinite Line {1, 3} = Nzb Using Progression Rzb;
Transfinite Line {2, 4} = Nxb Using Progression Rxb;

Line Loop(5) = {-1, 2, 3, 4};
Plane Surface(6) = {5};
Transfinite Surface {6};
Recombine Surface {6};

//     First conic section
b1 = (Dc-D)/(2*Tan(theta1*Pi/180));

Point(5) = {0, D/2, -b1, mes};
Point(6) = {0, 0, -b1, mes};

Line(5) = {2, 1};
Line(6) = {2, 5};
Line(7) = {5, 6};
Line(8) = {6, 1};

Transfinite Line {5, 7} = Nzb Using Progression Rzb;
Transfinite Line {6, 8} = Nxb Using Progression Rxb;

Line Loop(10) = {-5, 6, 7, 8};
Plane Surface(12) = {10};
Transfinite Surface {12};
Recombine Surface {12};

//     Second conic section
b2 = (D-Du)/(2*Tan(theta2*Pi/180));

Point(7) = {0, Du/2, -(b1+b2), mes};
Point(8) = {0, 0, -(b1+b2), mes};

Line(9) = {5, 6};
Line(10) = {5, 7};
Line(11) = {7, 8};
Line(12) = {8, 6};

Transfinite Line {9, 11} = Nzb Using Progression Rzb;
Transfinite Line {10, 12} = Nxb Using Progression Rxb;

Line Loop(15) = {-9, 10, 11, 12};
Plane Surface(18) = {15};
Transfinite Surface {18};
Recombine Surface {18};

//     Overflow
Point(9) = {0, 0, L1-VF, mes};
Point(10) = {0, Do/2, L1-VF, mes};
Point(11) = {0, Do/2, L1+ho, mes};
Point(12) = {0, 0, L1+ho, mes};

Line(13) = {10, 9};
Line(14) = {10, 11};
Line(15) = {11, 12};
Line(16) = {12, 9};

Transfinite Line {13, 15} = Nzo Using Progression Rzo;
Transfinite Line {14, 16} = Nxo Using Progression Rxo;

Line Loop(20) = {-13, 14, 15, 16};
Plane Surface(24) = {20};
Transfinite Surface {24};
Recombine Surface {24};
/*
//     First inlet pipe
x1 = Sqrt(Dc*Dc/4 - (Dc-Di)*(Dc-Di)/4);
x2 = Sqrt(Dc*Dc/4 - (Dc/2-Di)*(Dc/2-Di));

Point(13) = {-x1, -(Dc-Di)/2, L1-Di/2, mes};
Point(14) = {-x1, -(Dc-Di)/2, L1, mes};
Point(15) = {0, -Dc/2, L1-Di/2, mes};
Point(16) = {-x1, -(Dc-Di)/2, L1-Di, mes};
Point(17) = {-x2, -(Dc/2-Di), L1-Di/2, mes};
Point(18) = {0, 0, L1-Di/2, mes};

pList1[0] = 15;
pList2[0] = 14;
pList3[0] = 16;
pList4[0] = 15;

For i In {1:npi}

    yi1 = 0 + i*Di/(2*(npi+1));
    xi1 = Sqrt(Dc*Dc/4 - (Dc/2-yi1)*(Dc/2-yi1));
    zi1 = Sqrt(Di*Di/4 - (Di/2-yi1)*(Di/2-yi1));

    yi2 = Di/2 + i*Di/(2*(npi+1));
    xi2 = Sqrt(Dc*Dc/4 - (Dc/2-yi2)*(Dc/2-yi2));
    zi2 = Sqrt(Di*Di/4 - (Di/2-yi2)*(Di/2-yi2));

    pList1[i] = newp;
    Point(pList1[i]) = {-xi1, -(Dc/2-yi1), L1-Di/2+zi1, mes};
    
    pList2[i] = newp;
    Point(pList2[i]) = {-xi2, -(Dc/2-yi2), L1-Di/2+zi2, mes};

    pList3[i] = newp;
    Point(pList3[i]) = {-xi2, -(Dc/2-yi2), L1-Di/2-zi2, mes};
    
    pList4[i] = newp;
    Point(pList4[i]) = {-xi1, -(Dc/2-yi1), L1-Di/2-zi1, mes};
    
EndFor

pList1[npi+1] = 14;
pList2[npi+1] = 17;
pList3[npi+1] = 17;
pList4[npi+1] = 16;

np = newp;

Spline(17) = pList1[];
Spline(18) = pList2[];
Spline(19) = pList3[];
Spline(20) = pList4[];

Line(21) = {14, 13};
Circle(22) = {17, 18, 13};
Line(23) = {16, 13};
Circle(24) = {15, 18, 13};

Point(np) = {-hi, -(Dc-Di)/2, L1-Di/2, mes};
Point(np+1) = {-hi, -(Dc-Di)/2, L1-Di, mes};
Point(np+2) = {-hi, -(Dc/2-Di), L1-Di/2, mes};
Point(np+3) = {-hi, -(Dc-Di)/2, L1, mes};
Point(np+4) = {-hi, -Dc/2, L1-Di/2, mes};

Circle(25) = {np+2, np, np+3};
Circle(26) = {np+3, np, np+4};
Circle(27) = {np+1, np, np+4};
Circle(28) = {np+2, np, np+1};

Line(29) = {np, np+3};
Line(30) = {np, np+2};
Line(31) = {np, np+1};
Line(32) = {np, np+4};

Line(33) = {np+1, 16};
Line(34) = {np+2, 17};
Line(35) = {np+3, 14};
Line(36) = {np+4, 15};

Transfinite Line {21, 22, 23, 24, 29, 30, 31, 32} = Nzo Using Progression Rzo;

Transfinite Line {33, 34, 35, 36} = Nxo Using Progression Rxo;

Line Loop(40) = {33, 19, -34, 28};

Line Loop(45) = {34, -18, -35, -25};

Line Loop(50) = {35, -17, -36, -26};

Line Loop(55) = {36, 20, -33, 27};

Line Loop(60) = {19, -18, -17, 20};

Line Loop(65) = {-28, 25, 26, -27};

Surface(41) = {40};
Transfinite Surface {41};
Recombine Surface {41};

Surface(46) = {45};
Transfinite Surface {46};
Recombine Surface {46};

Surface(51) = {50};
Transfinite Surface {51};
Recombine Surface {51};

Surface(56) = {55};
Transfinite Surface {56};
Recombine Surface {56};

Surface(61) = {60};
Transfinite Surface {61};
Recombine Surface {61};



Surface(66) = {65};
Transfinite Surface {66};
Recombine Surface {66};

voli1[] = Extrude {-hi, 0, 0} {
    Surface{61};
    Layers{nolo};
    Recombine;
};



Line Loop(p+5) = {-(p+2), -(p+3), p, p+1};
//Line Loop(p+5) = {p, p+4, -(p+3)};
//Line Loop(p+5) = {p, p+1, -(p+2)};
//Line Loop(p+4) = {-(p+1), cn, -(cn+1), p};
Surface(p+6) = {p+5};
Transfinite Surface {p+6};
Recombine Surface {p+6};
*/

x1 = Sqrt(Dc*Dc/4 - (Dc-Di)*(Dc-Di)/4);
x2 = Sqrt(Dc*Dc/4 - (Dc/2-Di)*(Dc/2-Di));

Point(17) = {-x1, -(Dc-Di)/2, L1-Di/2, mes};
Point(18) = {-x1, -(Dc-Di)/2, L1, mes};
Point(19) = {0, -Dc/2, L1-Di/2, mes};
Point(20) = {-x1, -(Dc-Di)/2, L1-Di, mes};
Point(21) = {-x2, -(Dc/2-Di), L1-Di/2, mes};
/*
Point(22) = {x1, (Dc-Di)/2, L1-Di/2, mes};
Point(23) = {x1, (Dc-Di)/2, L1, mes};
Point(24) = {0, Dc/2, L1-Di/2, mes};
Point(25) = {x1, (Dc-Di)/2, L1-Di, mes};
Point(26) = {x2, (Dc/2-Di), L1-Di/2, mes};


Point(13) = {-hi, -(Dc-Di)/2, L1-Di, mes};
Point(14) = {-hi, -(Dc/2-Di), L1-Di/2, mes};
Point(15) = {-hi, -(Dc-Di)/2, L1, mes};
Point(16) = {-hi, -Dc/2, L1-Di/2, mes};



Line Loop(25) = {13, 14, 15, 16};//create circle first
Plane Surface(30) = {25};
Transfinite Surface {30};
Recombine Surface {30};
*/


/*
Line(17) = {28, 27};
Line(18) = {28, 29};
Line(19) = {29, 30};
Line(20) = {30, 27};

Transfinite Line {17} = Nzo Using Progression Rzo;
Transfinite Line {18} = Nxo Using Progression Rxo;
Transfinite Line {19} = Nzo Using Progression Rzo;
Transfinite Line {20} = Nxo Using Progression Rxo;

Line Loop(25) = {-17, 18, 19, 20};
Plane Surface(30) = {25};
Transfinite Surface {30};
Recombine Surface {30};

suri11[] = Extrude {{1, 0, 0}, {0, -(Dc-Di)/2, L1-Di/2}, Pi/2} {
    Line{17};
    Layers{nolo};
    Recombine;
};

suri12[] = Extrude {{1, 0, 0}, {0, -(Dc-Di)/2, L1-Di/2}, Pi/2} {
    Line{suri11[0]};
    Layers{nolo};
    Recombine;
};

suri13[] = Extrude {{1, 0, 0}, {0, -(Dc-Di)/2, L1-Di/2}, Pi/2} {
    Line{suri12[0]};
    Layers{nolo};
    Recombine;
};

suri14[] = Extrude {{1, 0, 0}, {0, -(Dc-Di)/2, L1-Di/2}, Pi/2} {
    Line{suri13[0]};
    Layers{nolo};
    Recombine;
};
*/
// ========= EXTRUSIONS (VOLUMES) ===============

// cylindrical extrusion
volcy1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {
    Surface{6};
    Layers{nolb};
    Recombine;
};

volcy2[] = Extrude {{0, 0, 1}, {0, 0, 0}, -Pi} {
    Surface{6};
    Layers{nolb};
    Recombine;
};

// first conic section
volc11[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {
    Surface{12};
    Layers{nolb};
    Recombine;
};

volc12[] = Extrude {{0, 0, 1}, {0, 0, 0}, -Pi} {
    Surface{12};
    Layers{nolb};
    Recombine;
};

// second conic section
volc21[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {
    Surface{18};
    Layers{nolb};
    Recombine;
};

volc22[] = Extrude {{0, 0, 1}, {0, 0, 0}, -Pi} {
    Surface{18};
    Layers{nolb};
    Recombine;
};

// overflow
volo1[] = Extrude {{0, 0, 1}, {0, 0, 0}, Pi} {
    Surface{24};
    Layers{nolo};
    Recombine;
};

volo2[] = Extrude {{0, 0, 1}, {0, 0, 0}, -Pi} {
    Surface{24};
    Layers{nolo};
    Recombine;
};
/*
// first inlet pipe
voli12[] = Extrude {hi, 0, 0} {
    Surface{suri11[1]};
    Layers{nolo};
    Recombine;
};

voli12[] = Extrude {hi, 0, 0} {
    Surface{suri12[1]};
    Layers{nolo};
    Recombine;
};

voli13[] = Extrude {hi, 0, 0} {
    Surface{suri13[1]};
    Layers{nolo};
    Recombine;
};

voli14[] = Extrude {hi, 0, 0} {
    Surface{suri14[1]};
    Layers{nolo};
    Recombine;
};
*/
//+
//vol[] = BooleanDifference{ Volume{9}; Volume{10}; Volume{11}; Volume{12}; Delete; }{ Volume{1};};


Cylinder(20) = {0, -(Dc-Di)/2, L1-Di/2, -hi, 0, 0, Di/2, 2*Pi};

BooleanDifference{ Volume{20}; Delete; }{ Volume{1}; }
/*
Line Loop(59) = {60, 64, -61};

Plane Surface(62) = {59};
Transfinite Surface {62};
Recombine Surface {62};

Recursive Delete {
  Volume{20}; 
}

np = newp;

Point(np) = {-h1, -(Dc-Di)/2, L1-Di/2, mes};
Point(np+1) = {-h1, -(Dc-Di)/2, L1, mes};
Point(np+2) = {-h1, -(Dc/2), L1-Di/2, mes};
Point(np+3) = {-h1, -(Dc-Di)/2, L1-Di, mes};
Point(np+4) = {-h1, -(Dc/2-Di), L1-Di/2, mes};

Circle(80) = {np+1, np, np+2};
Circle(81) = {np+2, np, np+3};
Circle(82) = {np+3, np, np+4};
Circle(83) = {np+4, np, np+1};

Line(29) = {np, np+3};
Line(30) = {np, np+2};
Line(31) = {np, np+1};
Line(32) = {np, np+4};

Line(84) = {np+1, 18};
Line(85) = {np+2, 19};
Line(86) = {np+3, 20};
//Line(87) = {np+4, 21};

Line Loop(88) = {80, 81, 82, 83};
Line Loop(89) = {80, 85, 64, -84};
Line Loop(90) = {81, 86, 61, -85};
Line Loop(91) = {82, 83, 84, 60, -86};

Plane Surface(92) = {88};
Transfinite Surface {92};
Recombine Surface {92};

Surface(93) = {89};
Transfinite Surface {93};
Recombine Surface {93};

Surface(94) = {90};
Transfinite Surface {94};
Recombine Surface {94};

Surface(95) = {91};
Transfinite Surface {95};
Recombine Surface {95};


// Inlets
//Cylinder(20) = {0, -(Dc-Di)/2, L1-Di/2, -hi, 0, 0, Di/2, 2*Pi};
//Cylinder(21) = {0, (Dc-Di)/2, L1-Di/2, hi, 0, 0, Di/2, 2*Pi};
/*
BooleanDifference{
    Volume{voli11[1]}; Volume{voli12[1]};
    Volume{voli13[1]}; Volume{voli14[1]};
    Delete; }{ Volume{1}; 
}
*/
//BooleanDifference{ Volume{20}; Delete; }{ Volume{1}; }
//BooleanIntersection{ Volume{20}; Delete; }{ Volume{1}; }
/*
Transfinite Surface {64} = {21, 18, 20, 19};
Recombine Surface {64};

Transfinite Surface {65} = {26, 23, 24, 25};
Recombine Surface {65};

Extrude {-hi, 0, 0} {
    Surface{64};
    Layers{nolo};
    Recombine;
}

BooleanIntersection{ Volume{20}; Delete; }{ Volume{22}; Delete; }


Extrude {hi, 0, 0} {
    Surface{65};
    Layers{nolo};
    Recombine;
}

//+
Recursive Delete {
  Volume{23}; Volume{21}; 
}
//+
Recursive Delete {
  Volume{20}; 
}
*/
//+
/*
Transfinite Surface {64} = {20, 21, 18, 19};
//+
//Transfinite Surface {59} = {13, 14, 15, 16};
//+
Transfinite Volume{11} = {16, 27, 13, 20, 19, 17};
//+
Transfinite Volume{12} = {14, 27, 13, 21, 17, 20};
//+
Transfinite Volume{9} = {14, 27, 15, 21, 17, 18};
*/
//+
Surface{57, 58, 59} In Volume{20};

Transfinite Line {57} = Nzo Using Progression Rzo;
Transfinite Line {58, 59, 60, 61} = Nxo Using Progression Rxo;
