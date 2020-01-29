
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
drc = 7.5;     	// Diameter of radial columns
dcc = 7;     	// Diameter of central column
wlf = 13; 	//Width lower flange/plate
gaf = 0.8; 	    // Gap factor (ratio of drc)
hgt = 6.35;	    // Height of columns
xO = 0;		// Model origin x-axis
yO = 0;		// Model origin y-axis
zO = -6.356260514543988;	// Model origin z-axis
nel_rrc = 16;    // Number of elements around cylinder circ.
vdist = 0.5;

// ----------------------------------------------------------------------------
// 	DERIVED UNITS
// ----------------------------------------------------------------------------

rrc = drc/2;
rcc = dcc/2;
dx  = (1 + gaf)*rrc;
dz = hgt;
nop_flat1 = nel_rrc/4 + 1;
nop_flat2 = nel_rrc/2 + 1;
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
Point(4) = {p02x,  p05y,  zO, lc};

// ----------------------------------------------------------------------------
// 	FIRST RADIAL COLUMN LOWER FLANGE (Z-POSITIVE FACE)
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

Point(5) = {  p12x,    yO,  zO, lc};
Point(13) = {  p13x,    p13y,  zO, lc};
Point(18) = {  p18x,    yO,  zO, lc};

// ----------------------------------------------------------------------------
// 	Z-NEGATIVE FACE LOWER FLANGE
// ----------------------------------------------------------------------------

Point(20) = {  xO,    yO,  zO - vdist, lc};
Point(21) = {p03x,  p03y,  zO - vdist, lc};
Point(22) = {  p12x,    yO,  zO - vdist, lc};
Point(23) = {  p13x,    p13y,  zO - vdist, lc};
Point(24) = {p02x,    yO,  zO - vdist, lc};
Point(45) = {p02x,  p05y,  zO - vdist, lc};


// ----------------------------------------------------------------------------
// 	FLANGES; LINES >> SURFACE
// ----------------------------------------------------------------------------

Line(1) = {20, 21};
Line(2) = {21, 45};
Line(3) = {45, 24};
Line(4) = {24, 20};
Line(5) = {45, 23};
Line(6) = {23, 22};
Line(7) = {22, 24};
Line(8) = {24, 45};
Circle(9) = {3, 0, 1};
Line(10) = {1, 2};
Line(11) = {2, 4};
Line(12) = {4, 3};
Line(13) = {2, 18};
Circle(14) = {18, 5, 13};
Line(15) = {13, 4};
Curve Loop(1) = {9, 10, 11, 12};
Plane Surface(1) = {1};
Curve Loop(2) = {11, -15, -14, -13};
Plane Surface(2) = {2};
Curve Loop(3) = {1, 2, 3, 4};
Plane Surface(3) = {3};
Curve Loop(4) = {7, -3, 5, 6};
Plane Surface(4) = {4};

// ----------------------------------------------------------------------------
// 	FLANGES; MESH PROPERTIES
// ----------------------------------------------------------------------------

//+
Transfinite Curve {13, 10, 14, 15, 11, 12, 9, 7, 4, 6, 5, 8, 3, 2, 1} = nop_flat1 Using Progression 1;
//Transfinite Curve {15} = nop_flat1 Using Progression .7;
//Transfinite Curve {14} = nop_flat1 Using Progression 1.3;

Transfinite Surface {1};
Transfinite Surface {2};
Transfinite Surface {3};
Transfinite Surface {4};

Recombine Surface {1, 2, 3, 4};


// ----------------------------------------------------------------------------
// 	CYLINDERS; EXTRUDE AND APPLY MESH PROPERTIES
// ----------------------------------------------------------------------------


Extrude {0, 0, dz} {
  Curve{9}; Curve{14};
}
Transfinite Curve {21,22} = nop_height1 Using Progression 1;
Transfinite Curve {17,18} = nop_height2 Using Progression 1;
Transfinite Surface {19};
Transfinite Surface {23};
Recombine Surface {19, 23};

// ----------------------------------------------------------------------------
// 	COPY COLUMNS ALONG RADIAL
// ----------------------------------------------------------------------------

//+
Symmetry {1, 0, 0, -p12x} {
  Duplicata { Surface{23}; Surface{2}; Surface{4}; }
}
//+
Symmetry {1, 0, 0, -1.5*p12x} {
  Duplicata { Surface{24}; Surface{2}; Surface{23}; Surface{29}; Surface{34}; Surface{4}; }
}

