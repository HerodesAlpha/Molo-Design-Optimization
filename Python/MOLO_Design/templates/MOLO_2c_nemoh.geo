
// ----------------------------------------------------------------------------
//	MOLO MESH
// 	PREPARED BY EIVIND SØNJU
// ----------------------------------------------------------------------------

// SetFactory("OpenCASCADE");
Geometry.CopyMeshingMethod = 1;

// ----------------------------------------------------------------------------
// 	INPUT
// ----------------------------------------------------------------------------

lc = 1e-2;	        //
drc = #dia_rc#;     // Diameter of radial columns
dcc = #dia_hc#;     // Diameter of central column
gaf = #gap#; 	    // Gap factor (ratio of drc)
hgt = #hgt#;	    // Height of columns
nel_rrc = #nel#;    // Number of elements around cylinder circ.

// ----------------------------------------------------------------------------
// 	DERIVED UNITS
// ----------------------------------------------------------------------------

xO = 0;		        // Model origin x-axis
yO = 0;		        // Model origin y-axis
zO = -hgt;	        // Model origin z-axis


rrc = drc/2;
rcc = dcc/2;
dx  = (1 + gaf)*rrc;
dz = hgt;
nop_flat1 = nel_rrc/8 + 1;
nop_flat2 = nel_rrc/4 + 1;
nop_height1 = Floor(dz/(Pi*drc/nel_rrc))+1;
nop_height2 = Floor(dz/(Pi*drc/(2*nel_rrc)))+1;

// ----------------------------------------------------------------------------
// 	CENTRAL HUB LOWER FLANGE
// ----------------------------------------------------------------------------

p01x = xO + rcc;
p02x = xO + dx;
p03x = xO + rcc/2;
p03y = yO + Sqrt(3)/2 * rcc;
// Need p05 to solve for p04
alpha=1.0;
p05x = xO + (1 + Sqrt(5))/2 * rrc*alpha; 
p05y = yO + rrc*alpha; 
// Set up line equations
m4 = Sqrt(3);
m5 = -1/2;
b4 = 0;
b5 = p05y -m5*p05x;
// Find intersection between line 4 and line 5
p04x = xO + (-b5/(m5-m4)); 
p04y = yO + p04x*m4;
p09x = xO + Sqrt(3)/2 * rcc;
p09y = yO + rcc/2;
p10x = (p01x/2 + p03x/2 + p09x)/3;
p10y = (yO + p03y/2 + p09y)/3;
p25y = yO + rrc/2; 

Point(0) = {  xO,    yO,  zO, lc};
Point(1) = {p01x,    yO,  zO, lc};
Point(2) = {p02x,    yO,  zO, lc};
Point(3) = {p03x,  p03y,  zO, lc};
Point(4) = {p04x,  p04y,  zO, lc};
Point(5) = {p05x,  p05y,  zO, lc};
Point(6) = {p02x,  p05y,  zO, lc};
Point(7) = {p01x/2,    yO,  zO, lc};
Point(8) = {p03x/2,  p03y/2,  zO, lc};
Point(9) = {p09x,  p09y,  zO, lc};
Point(10) = {p10x,  p10y,  zO, lc};
Point(11) = {p02x,  p25y,  zO, lc};

Printf("p1 = (%f, %f)",p01x,yO);
Printf("p2 = (%f, %f)",p02x,yO);
Printf("p3 = (%f, %f)",p03x,p03y);
Printf("p4 = (%f, %f)",p04x,p04y);
Printf("p5 = (%f, %f)",p05x,p05y);
Printf("p6 = (%f, %f)",p02x,p05y);

Line(1) = {0, 7};
Line(2) = {7, 1};
Line(3) = {1, 2};
Line(4) = {2, 11};
Line(5) = {11, 6};
Line(6) = {6, 5};
Line(7) = {5, 4};
Line(8) = {4, 3};
Line(9) = {3, 8};
Line(10) = {8, 0};
Line(11) = {7, 10};
Line(12) = {8, 8};
Line(13) = {10, 8};

Circle(14) = {1, 0, 9};
Circle(15) = {9, 0, 3};

Line(16) = {10, 9};
Line(17) = {9, 11};
Line(18) = {9, 5};

Curve Loop(1) = {10, 1, 11, 13};
Plane Surface(1) = {1};
Curve Loop(2) = {2, 14, -16, -11};
Plane Surface(2) = {2};
Curve Loop(3) = {16, 15, 9, -13};
Plane Surface(3) = {3};
Curve Loop(4) = {3, 4, -17, -14};
Plane Surface(4) = {4};
Curve Loop(5) = {17, 5, 6, -18};
Plane Surface(5) = {5};
Curve Loop(6) = {18, 7, 8, -15};
Plane Surface(6) = {6};

Transfinite Curve {16, 8, 9, 10, 11, 12, 13, 14, 15, 7, 17, 18, 1, 2, 3, 4, 5, 6} = nop_flat1 Using Progression 1;
Transfinite Surface {6};
Transfinite Surface {5};
Transfinite Surface {4};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {1};

Recombine Surface {1, 2, 3, 4, 5, 6};

// ----------------------------------------------------------------------------
// 	FIRST RADIAL COLUMN LOWER FLANGE
// ----------------------------------------------------------------------------

p11x = xO + 2*dx - rrc;
p12x = xO + 2*dx;
p13x = p12x;
p13y = yO + rrc;
p14x = p12x - rrc/2;
p14y = yO + rrc/2;
p15x = p12x - Sqrt(2)/2*rrc;
p15y = yO + Sqrt(2)/2*rrc;
p18x = xO + 2*dx - rrc;

Point(12) = {  p12x,    yO,  zO, lc};
Point(13) = {  p13x,    p13y,  zO, lc};
Point(14) = {  p14x,    p14y,  zO, lc};
Point(15) = {  p15x,    p15y,  zO, lc};
Point(16) = {  p14x,    yO,  zO, lc};
Point(17) = {  p12x,    p14y,  zO, lc};
Point(18) = {  p18x,    yO,  zO, lc};

Line(19) = {2, 18};
Line(20) = {18, 16};
Line(21) = {16, 12};
Line(22) = {12, 17};
Line(23) = {17, 13};
Line(24) = {13, 6};
Line(25) = {16, 14};
Line(26) = {14, 17};
Line(27) = {14, 15};
Line(28) = {15, 11};
Circle(29) = {18, 12, 15};
Circle(30) = {15, 12, 13};

Curve Loop(7) = {28, 5, -24, -30};
Plane Surface(7) = {7};
Curve Loop(8) = {4, -28, -29, -19};
Plane Surface(8) = {8};
Curve Loop(9) = {20, 25, 27, -29};
Plane Surface(9) = {9};
Curve Loop(10) = {21, 22, -26, -25};
Plane Surface(10) = {10};
Curve Loop(11) = {26, 23, -30, -27};
Plane Surface(11) = {11};

Transfinite Curve {21, 30, 29, 27, 26, 25, 23, 22, 20} = nop_flat1 Using Progression 1;
Transfinite Curve {28, 24, 19} = nop_flat2 Using Progression 1;

Transfinite Surface {7};
Transfinite Surface {11};
Transfinite Surface {10};
Transfinite Surface {9};
Transfinite Surface {8};

Recombine Surface {7, 11, 10, 9, 8};

Extrude {0, 0, dz} {
  Curve{30}; Curve{29};
}
Extrude {0, 0, dz} {
  Curve{14}; Curve{15};
}
Transfinite Curve {33, 32, 36} = nop_height1 Using Progression 1;
Transfinite Curve {45, 41, 40} = nop_height2 Using Progression 1;
Transfinite Surface {46};
Transfinite Surface {42};
Transfinite Surface {38};
Transfinite Surface {34};
Recombine Surface {46, 42, 38, 34};
