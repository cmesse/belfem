Mesh.SaveAll = 1;

m = 1;
cm = 10^-2;
mm = 10^-3;
um = 10^-6;
mu0 = 4e-7*Pi;  // [Hm⁻¹]
deg = Pi/180.;

p_Flag = "Input/01Flags/";
p_Geo  = "Input/02Geometry/";
p_Ele  = "Input/03Eletric/";


DefineConstant[	
// Define air region geometrical parameters
	air_r = { 4*cm, Name StrCat[p_Geo, "01Air/01radius"], visible 0 }
// Flags
	flag_time =  { 1, Choices { 0 = "hamonic", 1 = "time"}, Name StrCat[p_Flag, "01Time regime"] }
	flag_NL =  { 1, Choices { 0, 1}, Name StrCat[p_Flag, "01Nonlinear"] }
// Mesh paramenters
	sf = { 4*mm, Name StrCat[p_Geo, "01Air/02Mesh/01sf"] }
// Geometrical parameters
	thin_dx = { 4*mm, Name StrCat[p_Geo, "02Thin/03Width"], visible 0 } 
	thin_dy = { 1*um, Choices { 10*um, 100*um}, Name StrCat[p_Geo, "02Thin/04Thickness"], visible 0 }
	thin_dz = { 1*mm, Name StrCat[p_Geo, "02Thin/05Lenght"], visible 0 }
//Mesh parameters
	thin_sf_in = { 4*mm/200, Name StrCat[p_Geo, "02Thin/06Mesh/01sf_thin_in"] }
	thin_sf_out = { 4*mm/200, Name StrCat[p_Geo, "02Thin/06Mesh/01sf_thin_out"] }
// Electrical parameters
	thin_J = { 90/(4*mm*1*um), Name StrCat[p_Geo, "02Thin/07Eletric/01J[A par m^2]"] }
	thin_sigma = {5.96e7, Name StrCat[p_Geo, "02Thin/07Eletric/02Sigma"] }
	thin_mu = {4e-7*Pi, Name StrCat[p_Geo, "02Thin/07Eletric/03Mu"] }
// Coil number
	n_coils= { 19, Name StrCat[p_Geo, "02Thin/08Number of coils"], visible 0 }
	coils_dist= { 250*um, Name StrCat[p_Geo, "02Thin/11Coils distance"], visible 0 }
	N_ele = { 1., Name StrCat[p_Geo, "02Thin/06Number elements in the TS model"] }
];
	
I_thin = thin_J*(thin_dx*thin_dy); // Transport current

If (flag_time ==0) // time-harmonic simulation parameters
	DefineConstant[
		Freq = { 50, Choices { 50, 1000000, 10000000, 100000000}, Name StrCat["Input/Solver/", "02Freq [Hz]"] }
	];
ElseIf (flag_time ==1)  // time-dependent simulation parameters
	DefineConstant[
	source_type = { 0, Choices { 0 = "sinus", 1 = "ramp"}, Name StrCat["Input/Solver/", "01Source type"] }
	];
	
	If (source_type==0) // sinus source
		DefineConstant[
			Freq = { 50, Name StrCat["Input/Solver/", "02Freq [Hz]"] }
			periods = {1, Min 0.1, Max 2.0, Step 0.05,
			Name "Input/Solver/0Periods to simulate"},
			time0 = 0, // initial time
			time1 = periods * (1 / Freq), // final time
			dt = {25e-6, Min 25e-6, Max 25e-6, Step 1e-6,
			  Name "Input/Solver/1Time step [s]"}
			dt_max = {25e-6,
			  Name "Input/Solver/2Maximum time step [s]"},
			tol_abs = {1e-6,
			  Name "Input/Solver/3Absolute tolerance on nonlinear residual"},
			tol_rel = {1e-5,
			  Name "Input/Solver/3Relative tolerance on nonlinear residual"},
			iter_max = {12,
			  Name "Input/Solver/Maximum number of nonlinear iterations"},
			visu = {1, Choices{0, 1}, AutoCheck 0,
			  Name "Input/Solver/Visu", Label "Real-time visualization"}
		];

	ElseIf (source_type==1) //ramp source
		DefineConstant[
			time1 = {1/50, Name StrCat["Input/Solver/", "0Final time"] },
			time0 = 0., // initial time
			time_N = {1000, Name StrCat["Input/Solver/", "1Time step [s]"]},
			dt = 25e-6//(time1-time0)/(time_N-1)
			tol_abs = {1e-10,
			  Name StrCat["Input/Solver/", "3Absolute tolerance on nonlinear residual"] },
			tol_rel = {1e-5,
			  Name StrCat["Input/Solver/", "3Relative tolerance on nonlinear residual"] },
			iter_max = {12,
			  Name StrCat["Input/Solver/", "Maximum number of nonlinear iterations"] },
			dt_max = {0.5e-3,
			  Name StrCat["Input/Solver/", "2Maximum time step [s]"] },
		];
	EndIf
EndIf

// REGIONS IDs
AIR = 100000;
AIR_BND = 100000+1;

For i In {1:n_coils}
	COIL~{i} = 200000+i*1000;
EndFor

// domain width
a = air_r*2 ;

// domain height
b = a ;

// tape width
c = thin_dx ;

// distance between tapes
d = coils_dist ;

// domain mesh size
e = sf ;

// tape mesh size
f = thin_sf_in ; //40 20
g = thin_sf_out ; //40 200
//f = 0.025 ;
//g = 0.005 ;


n = n_coils ; // must be an odd number! 
dn = 0.5 * ( n-1) - 1 ;

// domain corner nodes

// circle nodes
Point(1) = {      0,-0.5*b, 0, e } ;
Point(2) = {  0.5*a,0, 0, e } ;
Point(3) = {      0, 0.5*b, 0, e } ;
Point(4) = { -0.5*a, 0, 0, e } ;

// tapes
h = - (n-1) * d * 0.5 ;

p = 4 ;
For k In {1:n}
    p += 1 ;
    Point(p) = { -0.5*c, h, 0, g } ;
    p += 1 ;
    Point(p) = { 0, h, 0, f } ;
    p += 1 ;
    Point(p) = { 0.5*c,h, 0, g } ;
    h += d ;
EndFor


// center point
m = 33; //Generalize this

// domain edges
Circle(1) = {1, m, 2};
Circle(2) = {2, m, 3};
Circle(3) = {3, m, 4};
Circle(4) = {4, m, 1};

// tapes
l = 4 ;
p = 5 ;
ta = l + 1 ;
For k In {1:n}
    l+=1;
    Line(l) = {p, p+1};
    l+=1;
    Line(l) = {p+1, p+2};
    p+=3;
EndFor
tn = l ;

// cuts
l+=1;
ca = l ;
Line(l) = {1, 6};
p = 6 ;
For k In {1:n-1}
	l+=1;
	Line(l) = {p, p+3};
	p+=3;
EndFor

l+=1 ;
Line(l) = {m+1, 2};
hr = l ;

l+=1;
Line(l) = {p, 3};
cn = l ;

l+=1 ;
Line(l) = {m-1, 4};
hl = l ;


// left edges
la = l+1 ;
p = 5 ;
For k In {1:n-1}
	l+=1;
	Line(l) = {p, p+3};
	p+=3;
EndFor
ln = l ;
ra = l + 1 ;
p += 2 ;
For k In {1:n-1}
	l+=1;
	Line(l) = {p, p-3};
	p-=3;
EndFor
rn = l ;

// Air Domains

w = la + dn ;

Line Loop(1) = {4, -5, 43, 64, 65:73}; //Generalize this
Plane Surface(1) = {1};


Line Loop(2) = {1, -6, -43, -62, 92:100}; //Generalize this
Plane Surface(2) = {2};

B = 2 ;
bl = 5 ;
br = 6 ;
ll = la ;
rr = rn ;
cc = ca + 1 ;
tl = 7 ;
tr = 8 ;



For k In {1:n-1}
	B+=1;
	Line Loop(B) = {bl, cc, -tl, -ll};
	Plane Surface(B) = {B};
	B+=1;
	Line Loop(B) = {br, -rr, -tr, -cc};
	Plane Surface(B) = {B};
	bl+=2;
	br+=2;
	tl+=2;
	tr+=2;
	cc+=1;
	ll+=1;
	rr-=1;
EndFor

B+=1;
Line Loop(B) = {3, -hl, ln-dn:ln, tn-1, cn};
Plane Surface(B) = {B};

B+=1;
Line Loop(B) = {2, -cn, tn, ra:ra+dn, hr};
Plane Surface(B) = {B};

