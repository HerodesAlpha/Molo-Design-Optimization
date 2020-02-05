
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
drc = 8.6;     	// Diameter of radial columns
dcc = 8.6;     	// Diameter of central column
wlf = 10; 	//Width lower flange/plate
overlength = 1; // Lower flange overlength
gaf = 0.8; 	    // Gap factor (ratio of drc)
hgt = 6.35;	    // Height of columns
x0 = 0;		// Model origin x-axis
y0 = 0;		// Model origin y-axis
z0 = -6.356260514543988;	// Model origin z-axis
nel_rrc = 16;    // Number of elements around cylinder circ.
vdist = 5;

// ----------------------------------------------------------------------------
// 	DERIVED UNITS
// ----------------------------------------------------------------------------

rrc = drc/2;
rcc = dcc/2;
dx  = (1 + gaf)*drc; // Distance between column centers
dz = hgt;
Printf("dz = %f",dz);
l_el = Pi*drc/nel_rrc;
Printf("l_el = %f",l_el);

val=Floor((wlf/2-rrc)/l_el);
nop_flat0 =  (val < 2) ? 2 : val;
nop_flat1 = nel_rrc/8 + 1;
//nop_flat2 = nel_rrc/2 + 1;
nop_flat3 = nop_flat1 + nop_flat0 - 1;

nop_x_edge = Floor(((1 + 2*gaf)*drc/2)/l_el) + 1;


nop_height_radial_column = Floor(dz/(l_el)) + 1;
Printf("dz/(l_el) = %f",dz/(l_el));
nop_height_center_column = Floor(dz/(l_el*dcc/drc)) + 1;
nop_lower_long_edge = Floor(2.5*dx/l_el) + 1;
nop_lower_short_edge = Floor(wlf/l_el) + 1;

// ----------------------------------------------------------------------------
// 	CENTRAL HUB LOWER FLANGE
// ----------------------------------------------------------------------------
Point(1) = {  x0,    y0,  z0, lc};

p02x = x0 + rcc;
Point(2) = {  p02x,  y0,  z0, lc};

p03x =  x0 + dx/2;
Point(3) = {  p03x,  y0,  z0, lc};

p04x = x0 + dx - rrc;
Point(4) = {  p04x,  y0,  z0, lc};

p05x = x0 + dx;
Point(5) = {  p05x,  y0,  z0, lc};

p06x = x0 + rcc/2;
p06y = y0 + rcc *Sqrt(3)/2;
Point(6) = {  p06x,  p06y,  z0, lc};

p07x = p05x;
p07y = y0 + rrc;
Point(7) = {  p07x,  p07y,  z0, lc};

p08x = x0 + (wlf/2)/Sqrt(3);
p08y = y0 + wlf/2;
Point(8) = {  p08x,  p08y,  z0, lc};

p09x = p03x;
p09y = y0 + wlf/2;
Point(9) = {  p09x,  p09y,  z0, lc};

p10x = p05x;
p10y = y0 + wlf/2;
Point(10) = {  p10x,  p10y,  z0, lc};

p11x = x0 + rcc *Sqrt(3)/2;
p11y = y0 + rcc/2;
Point(11) = {  p11x,  p11y,  z0, lc};

p12x = x0 + p03x;
p12y = y0 + rrc*1/2;
Point(12) = {  p12x,  p12y,  z0, lc};//+

//p13x = x0 + dx - rrc*Sqrt(2)/2;
//p13y = y0 + rrc*Sqrt(2)/2;
p13x = x0 + dx - rrc*Sqrt(3)/2;
p13y = y0 + rrc*1/2;
Point(13) = {  p13x,  p13y,  z0, lc};//+
//+
Line(1) = {2, 3};
//+
Line(2) = {3, 12};
//+
Line(3) = {12, 11};
//+
Circle(4) = {11, 1, 2};
//+
Line(5) = {12, 9};
//+
Line(6) = {9, 8};
//+
Line(7) = {8, 6};
//+
Circle(8) = {6, 1, 11};
//+
Line(9) = {3, 4};
//+
Circle(10) = {4, 5, 13};
//+
Line(11) = {13, 12};
//+
Circle(12) = {13, 5, 7};
//+
Line(13) = {7, 10};
//+
Line(14) = {10, 9};
//+
Curve Loop(1) = {1, 2, 3, 4};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {3, -8, -7, -6, -5};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {9, 10, 11, -2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {11, 5, -14, -13, -12};
//+
Plane Surface(4) = {4};
//+

// ----------------------------------------------------------------------------
// 	FLANGES; MESH PROPERTIES
// ----------------------------------------------------------------------------

Transfinite Curve {13, 12, 10, 14, 5, 11, 9, 2, 1, 4, 3, 6, 7, 8} = nop_flat1 Using Progression 1;
//+
Transfinite Surface {2} = {8, 11, 12, 9};
//+
Transfinite Surface {4} = {9, 12, 13, 10};
//+
Transfinite Surface {1};
//+
Transfinite Surface {3};
//+
Transfinite Curve {7, 13} = nop_flat0 Using Progression 1;
//+
Transfinite Curve {5} = nop_flat3 Using Progression 1;
//+
Transfinite Curve {6, 14, 3, 11, 1, 9} = nop_x_edge Using Progression 1;
//+
Transfinite Curve {14} = nop_x_edge Using Progression 0.9;
//+
Recombine Surface {4, 3, 1, 2};


