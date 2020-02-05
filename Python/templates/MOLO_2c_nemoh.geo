
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
drc = 8.0;     	// Diameter of radial columns
dcc = 8.0;     	// Diameter of central column
wlf = 10; 	//Width lower flange/plate
overlength = 1; // Lower flange overlength
gaf = 3.5; 	    // Gap factor (ratio of drc)
hgt = 12;	    // Height of columns
x0 = 0;		// Model origin x-axis
y0 = 0;		// Model origin y-axis
z0 = -12;	// Model origin z-axis
nel_rrc = 16;    // Number of elements around cylinder circ.
vdist = 0.5;

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
nop_flat0 =  (val < 2) ? 2 : val; // Y-directoon overwidth
nop_flat1 = nel_rrc/8 + 1; // 1/8 cirumferences
nop_flat2 = nel_rrc/16 + 1;// 1/16 cirumferences
nop_flat3 = nop_flat1 + nop_flat0 - 1;
val2 = Floor(overlength / l_el)+1;
nop_flat4 =   (val2 < 2) ? 2 : val2;
Printf("nop_flat4 = %f",nop_flat4);
nop_x_edge = Floor(((1 + gaf)*drc/4)/l_el) + 1;


nop_height_radial_column = Floor(dz/(l_el)) + 1;
Printf("dz/(l_el) = %f",dz/(l_el));
nop_height_center_column = Floor(dz/(l_el*dcc/drc)) + 1;
nop_lower_long_edge = Floor(2*dx/l_el) + 1;
nop_lower_short_edge = Floor((wlf/2)/l_el) + 1;

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

p11x = x0 + rcc*Sqrt(2)/2;
p11y = y0 + rcc*Sqrt(2)/2;
Point(11) = {  p11x,  p11y,  z0, lc};

p12x = x0 + p03x;
p12y = y0 + rrc*Sqrt(2)/2;
Point(12) = {  p12x,  p12y,  z0, lc};//+

p13x = x0 + dx - rrc*Sqrt(2)/2;
p13y = y0 + rrc*Sqrt(2)/2;
//p13x = x0 + dx - rrc*Sqrt(3)/2;
//p13y = y0 + rrc*1/2;
Point(13) = {  p13x,  p13y,  z0, lc};//+

p14x = p11x + (wlf/2)/Sqrt(3) -rcc/2;
p14y = p08y;
Point(14) = {  p14x,  p14y,  z0, lc};

p15x = p13x;
p15y = p08y;
Point(15) = {  p15x,  p15y,  z0, lc};





//+
Line(1) = {11, 14};
//+
Line(2) = {14, 8};
//+
Line(3) = {8, 6};
//+
Circle(4) = {6, 1, 11};
//+
Circle(5) = {11, 1, 2};
//+
Line(6) = {2, 3};
//+
Line(7) = {3, 12};
//+
Line(8) = {12, 11};
//+
Line(9) = {12, 9};
//+
Line(10) = {9, 14};
//+
Line(11) = {3, 4};
//+
Circle(12) = {4, 5, 13};
//+
Line(13) = {13, 12};
//+
Circle(14) = {13, 5, 7};
//+
Line(15) = {7, 10};
//+
Line(16) = {10, 15};
//+
Line(17) = {15, 13};
//+
Line(18) = {15, 9};
//+
Curve Loop(1) = {7, 8, 5, 6};
//+
Plane Surface(1) = {1};
//+
Curve Loop(2) = {11, 12, 13, -7};
//+
Plane Surface(2) = {2};
//+
Curve Loop(3) = {3, 4, 1, 2};
//+
Plane Surface(3) = {3};
//+
Curve Loop(4) = {1, -10, -9, 8};
//+
Plane Surface(4) = {4};
//+
Curve Loop(5) = {9, -18, 17, 13};
//+
Plane Surface(5) = {5};
//+
Curve Loop(6) = {17, 14, 15, 16};
//+
Plane Surface(6) = {6};
//+
Transfinite Curve {3, 1, 9, 17, 15} = nop_flat0 Using Progression 1;
//+
Transfinite Curve {4, 2} = nop_flat2 Using Progression 1;
//+
Transfinite Curve {5, 7, 12, 14, 16} = nop_flat1 Using Progression 1;
//+
Transfinite Curve {6, 8, 10, 11, 13, 18} = nop_x_edge Using Progression 1;
//+
Transfinite Surface {3};
//+
Transfinite Surface {4};
//+
Transfinite Surface {5};
//+
Transfinite Surface {6};
//+
Transfinite Surface {1};
//+
Transfinite Surface {2};
//+
Recombine Surface {6, 2, 5, 1, 4, 3};


// ----------------------------------------------------------------------------
// 	CYLINDERS; EXTRUDE AND APPLY MESH PROPERTIES
// ----------------------------------------------------------------------------
//+
Extrude {0, 0, dz} {
  Curve{4}; Curve{5}; Curve{12}; Curve{14};
}
//+
Transfinite Curve {28, 29, 33} = nop_height_radial_column Using Progression 1;
//+
Transfinite Curve {20, 21, 25} = nop_height_center_column Using Progression 1;
//+
Transfinite Surface {22};
//+
Transfinite Surface {26};

// ----------------------------------------------------------------------------
// 	Order nodes on surfaces
// ----------------------------------------------------------------------------
//+
//+
Transfinite Surface {3} = {6, 11, 14, 8};
//+
Transfinite Surface {4} = {14, 11, 12, 9};
//+
Transfinite Surface {5} = {12, 13, 15, 9};
//+
Transfinite Surface {6} = {15, 13, 7, 10};
//+
Transfinite Surface {2} = {12, 3, 4, 13};
//+
Transfinite Surface {1} = {3, 12, 11, 2};
//+
Transfinite Surface {22} = {6, 16, 18, 11};
//+
Transfinite Surface {26} = {11, 18, 21, 2};
//+
Transfinite Surface {30} = {4, 22, 24, 13};
//+
Transfinite Surface {34} = {24, 27, 7, 13};

Reverse Surface { 1,2,3,4,5,6,22,26,30,34 };

Recombine Surface {22, 26, 30, 34};

// ----------------------------------------------------------------------------
// 	COPY COLUMNS ALONG RADIAL
// ----------------------------------------------------------------------------//+

//+
Symmetry {1, 0, 0, -p05x} {
  Duplicata { Surface{5}; Surface{2}; Surface{6}; Surface{34}; Surface{30}; }
}
//+
Symmetry {1, 0, 0, -1.5*p05x} {
  Duplicata { Surface{2}; Surface{5}; Surface{6}; Surface{45}; Surface{35}; Surface{40}; Surface{30}; Surface{34}; Surface{50}; Surface{55}; }
}
//+
Translate {rcc + overlength -dx/2, 0, 0} {
  Point{108}; Point{118}; Point{126};
}

Transfinite Curve {61, 63, 67} = nop_flat4 Using Progression 1;

// ----------------------------------------------------------------------------
// 	Z-NEGATIVE FACE LOWER FLANGE
// ----------------------------------------------------------------------------

//+
Translate {0, 0, -vdist} {
  Duplicata { Point{1}; Point{8}; Point{126}; Point{108}; }
}

//+
Line(110) = {274, 277};
//+
Line(111) = {277, 276};
//+
Line(112) = {276, 275};
//+
Line(113) = {275, 274};
//+
Curve Loop(7) = {110, 111, 112, 113};
//+
Plane Surface(106) = {7};
//+
Transfinite Surface {106} = {274, 275, 276, 277};
//+
Transfinite Curve {110, 112} = nop_lower_long_edge Using Progression 1;
//+
Transfinite Curve {113, 111} = nop_lower_short_edge Using Progression 1;
//+
Recombine Surface {106};
//+
Symmetry {-Sqrt(3)/2, 1/2, 0, 0} {
  Duplicata { Surface{26}; Surface{22}; Surface{1}; Surface{3}; Surface{30}; Surface{4}; Surface{2}; Surface{5}; Surface{34}; Surface{55}; Surface{6}; Surface{50}; Surface{45}; Surface{40}; Surface{105}; Surface{35}; Surface{85}; Surface{80}; Surface{100}; Surface{90}; Surface{75}; Surface{95}; Surface{60}; Surface{70}; Surface{65}; Surface{106}; }
}
//+
Symmetry {-Sqrt(3)/2, -1/2, 0, 0} {
  Duplicata { Surface{114}; Surface{119}; Surface{129}; Surface{124}; Surface{134}; Surface{139}; Surface{144}; Surface{149}; Surface{154}; Surface{159}; Surface{164}; Surface{169}; Surface{174}; Surface{179}; Surface{184}; Surface{189}; Surface{194}; Surface{199}; Surface{204}; Surface{209}; Surface{214}; Surface{219}; Surface{224}; Surface{229}; Surface{234}; Surface{239}; }
}
//+

