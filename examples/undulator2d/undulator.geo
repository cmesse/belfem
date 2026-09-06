// Gregory Giard 07/19/2023
// Undulator benchmark geometric definition

Mesh.SaveAll = 1;

Include "ParamsGeo.dat";

// domain corner nodes
p=1;
Point(1) = {0,0,0,sf_core/2};
l = 0;

For k In {1:n_groove/2}
	p+=1;
	pd~{k} = p;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w/2,0,0,sf_core/2};
	l+=1;
	Line(l) = {p-1,p};
EndFor

p+=1;
Point(p) = {air_r,0,0,e};
l+=1;
Line(l) = {p-1,p};

For i In {1:n_groove/2}
	k = n_groove/2-(i-1);
	p+=1;
	pu~{k} = p;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w/2,Sqrt[air_r^2-((2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w/2)^2],0,e};
	l+=1;
	Circle(l) = {p-1,1,p};
EndFor

p+=1;
Point(p) = {0,air_r,0,e};
l+=1;
Circle(l) = {p-1,1,p};

// domain nodes in material
p+=1;
Point(p) = {0,gap/2+2*pole_h+core_h,0,sf_sc_center};
pul = p;
l+=1;
Line(l) = {p-1,p};
p+=1;
Point(p) = {0,gap/2+pole_h+core_h/2,0,sf_core};
pcl = p;
l+=1;
Line(l) = {p-1,p};
p+=1;
Point(p) = {0,gap/2,0,sf_sc_center};
pdl = p;
l+=1;
Line(l) = {p-1,p};

p+=1;
Point(p) = {pole_w/2,gap/2,0,sf_sc*2};

l+=1;
Line(l) = {p-1,1};

l+=1;
Line(l) = {p-1,p};
ldl = l;


// bottom pole lines
For k In {1:n_groove/2}
	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w,gap/2+pole_h,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+groove_w/2,gap/2+pole_h,0,sf_sc_center};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+k*groove_w,gap/2+pole_h,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+k*groove_w,gap/2,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k)*pole_w/2+k*groove_w,gap/2,0,sf_sc_center};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;

	If ( k == n_groove/2 )
		Point(p) = {(2*k+1)*pole_w/2+k*groove_w,gap/2,0,sf_sc_center};
	Else
		Point(p) = {(2*k+1)*pole_w/2+k*groove_w,gap/2,0,sf_sc*2};
	EndIf

	l+=1;
	Line(l) = {p-1,p};
EndFor
ldr = l;
pdr = p;

p+=1;
Point(p) = {pole_w/2,gap/2+2*pole_h+core_h,0,sf_sc*2};
l+=1;
Line(l) = {pul,p};
lul = l;

// top pole lines
For k In {1:n_groove/2}
	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w,gap/2+pole_h+core_h,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+0.5*groove_w,gap/2+pole_h+core_h,0,sf_sc_center};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+k*groove_w,gap/2+pole_h+core_h,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	Point(p) = {(2*k-1)*pole_w/2+k*groove_w,gap/2+2*pole_h+core_h,0,sf_sc*2};
	l+=1;
	Line(l) = {p-1,p};
	
	p+=1;
	Point(p) = {(2*k)*pole_w/2+k*groove_w,gap/2+2*pole_h+core_h,0,sf_sc_center};
	l+=1;
	Line(l) = {p-1,p};

	p+=1;
	If ( k == n_groove/2 )
		Point(p) = {(2*k+1)*pole_w/2+k*groove_w,gap/2+2*pole_h+core_h,0,sf_sc_center};
	Else
		Point(p) = {(2*k+1)*pole_w/2+k*groove_w,gap/2+2*pole_h+core_h,0,sf_sc*2};
	EndIf
	l+=1;
	Line(l) = {p-1,p};
EndFor

pur = p;
lur = l;

p+=1;
Point(p) = {(n_groove/2)*groove_w+(n_groove + 1)*pole_w/2,gap/2+pole_h+core_h/2,0,sf_core};
pcr = p;

l+=1;
lc1 = l;
l+=1;
lc2 = l;
Line(lc1)={pdr,pcr};
Line(lc2)={pcr,pur};

l+=1;
Line(l) = {pcl,pcr};
lc = l;

ll=2;

//Air surface
//Line Loop(1) = {1:n_groove + 3,lul:lur,-lc2,-lc1,-ldr:-ldl,(n_groove + 6)};
//Plane Surface(1) = {1};

//Core surface
Line Loop(2) = {n_groove + 4,lc,lc2,-lur:-lul};
Plane Surface(2) = {2};
Line Loop(3) = {n_groove + 5,ldl:ldr,lc1,-lc};
Plane Surface(3) = {3};

s = 3;
ll = 3;

// bottom conductors
For k In {1:n_groove/2}

	lsdl_[] = {};
	lsdr_[] = {};

	If (k == n_groove/2)
		n_tapes = n_tapes_ext;
		tape_space = tape_space_ext;
	Else
		n_tapes = n_tapes_int;
		tape_space = tape_space_int;
	EndIf

	For i In {1:n_tapes}

		p+=1;
		ptdl~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2,gap/2+pole_h-G10_h-i*tape_space,0,sf_sc};

		p+=1;
		ptdr~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w,gap/2+pole_h-G10_h-i*tape_space,0,sf_sc};

		// center point
		p+=1;
		ptdc~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w/2,gap/2+pole_h-G10_h-i*tape_space,0,sf_sc_center};

		l+=1;
		ltdl~{k}~{i} = l;
		Line(l) = {ptdl~{k}~{i},ptdc~{k}~{i}};

		l+=1;
		ltdr~{k}~{i} = l;
		Line(l) = {ptdc~{k}~{i},ptdr~{k}~{i}};


		If (i > 1)
			l+=1;
			Line(l) = {ptdl~{k}~{i-1},ptdl~{k}~{i}};
			lsdl~{k}~{i} = l;
			lsdl_[]+=lsdl~{k}~{i};
			
			l+=1;
			Line(l) = {ptdr~{k}~{i-1},ptdr~{k}~{i}};
			lsdr~{k}~{i} = l;
			lsdr_[]+=lsdr~{k}~{i};
			
			l+=1;
			Line(l) = {ptdc~{k}~{i-1},ptdc~{k}~{i}};
			lcd~{k}~{i} = l;
			
			ll+=1;
			Line Loop(ll) = {-ltdl~{k}~{i-1},lsdl~{k}~{i},
							ltdl~{k}~{i},-lcd~{k}~{i}};
			s+=1;
			Plane Surface(s) = {s};
			
			ll+=1;
			Line Loop(ll) = {-ltdr~{k}~{i-1},lcd~{k}~{i},
							ltdr~{k}~{i},-lsdr~{k}~{i}};
			s+=1;
			Plane Surface(s) = {s};
			
		EndIf

		If (i == n_tapes)
			l+=1;
			lgdl~{k} = l;
			Line(l) = {ptdl~{k}~{i},pdl+6*(k-1)+1};
			
			l+=1;
			lgdr~{k} = l;
			Line(l) = {ptdr~{k}~{i},pdl+6*(k-1)+5};
		
			l+=1;
			Line(l) = {ptdc~{k}~{i},pd~{k}};
			lcd~{k} = l;
			
			If (k == 1)
				ll+=1;
				Line Loop(ll) = {-ltdl~{k}~{i},lgdl~{k},-ldl,(ldl-1),k,-lcd~{k}};
				s+=1;
				Plane Surface(s) = {s};
			Else
				ll+=1;
				Line Loop(ll) = {-ltdl~{k}~{i},lgdl~{k},-(ldl+6*(k-1)),
								-(ldl+6*(k-1)-1),-lgdr~{k-1},-ltdr~{k-1}~{n_tapes_int}, 
								lcd~{k-1},k,-lcd~{k}};
				s+=1;
				Plane Surface(s) = {s};
			EndIf
			
			ll+=1;
			Line Loop(ll) = {lgdr~{k},-(ldl+6*(k-1)+4):-(ldl+6*(k-1)+1),
							-lgdl~{k},-lsdl_[],ltdl~{k}~{1},ltdr~{k}~{1},
							lsdr_[]};
			s+=1;
			Plane Surface(s) = {s};
			
		EndIf


	EndFor

EndFor

// top conductors
For k In {1:n_groove/2}

	lsul_[] = {};
	lsur_[] = {};

	If (k == n_groove/2)
		n_tapes = n_tapes_ext;
		tape_space = tape_space_ext;
	Else
		n_tapes = n_tapes_int;
		tape_space = tape_space_int;
	EndIf
	For i In {1:n_tapes}

		p+=1;
		ptul~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2,gap/2+pole_h+core_h+G10_h+i*tape_space,0,sf_sc};

		p+=1;
		ptur~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w,gap/2+pole_h+core_h+G10_h+i*tape_space,0,sf_sc};

		p+=1;
		ptuc~{k}~{i} = p;
		Point(p) = {(2*k-1)*pole_w/2+(k-1)*groove_w+(groove_w-tape_w)/2+tape_w/2,gap/2+pole_h+core_h+G10_h+i*tape_space,0,sf_sc_center};
		
		l+=1;
		ltul~{k}~{i} = l;
		Line(l) = {ptul~{k}~{i},ptuc~{k}~{i}};

		l+=1;
		ltur~{k}~{i} = l;
		Line(l) = {ptuc~{k}~{i},ptur~{k}~{i}};

		If (i > 1)
			l+=1;
			Line(l) = {ptul~{k}~{i-1},ptul~{k}~{i}};
			lsul~{k}~{i} = l;
			lsul_[]+=lsul~{k}~{i};
			
			l+=1;
			Line(l) = {ptur~{k}~{i-1},ptur~{k}~{i}};
			lsur~{k}~{i} = l;
			lsur_[]+=lsur~{k}~{i};
			
			l+=1;
			Line(l) = {ptuc~{k}~{i-1},ptuc~{k}~{i}};
			lcu~{k}~{i} = l;
			
			ll+=1;
			Line Loop(ll) = {ltul~{k}~{i-1},lcu~{k}~{i},
								-ltul~{k}~{i},-lsul~{k}~{i}};
			s+=1;
			Plane Surface(s) = {s};
			
			ll+=1;
			Line Loop(ll) = {ltur~{k}~{i-1},lsur~{k}~{i},
							-ltur~{k}~{i},-lcu~{k}~{i}};
			s+=1;
			Plane Surface(s) = {s};
		EndIf

		If (i == n_tapes)
			l+=1;
			lgul~{k} = l;
			Line(l) = {ptul~{k}~{i},pdr+6*(k-1)+1};
			
			l+=1;
			lgur~{k} = l;
			Line(l) = {ptur~{k}~{i},pdr+6*(k-1)+5};
		
			l+=1;
			Line(l) = {ptuc~{k}~{i},pu~{k}};
			lcu~{k} = l;
			
			If (k == 1)
				ll+=1;
				Line Loop(ll) = {ltul~{k}~{i},-lgul~{k},lul,(n_groove+4-k),(n_groove+3-k),lcu~{k}};
				s+=1;
				Plane Surface(s) = {s};
			Else
				ll+=1;
				Line Loop(ll) = {ltul~{k}~{i},-lgul~{k},(lul+6*(k-1)),
								(lul+6*(k-1)-1),lgur~{k-1},ltur~{k-1}~{n_tapes_int}, 
								-lcu~{k-1},(n_groove+3-k),lcu~{k}};
				s+=1;
				Plane Surface(s) = {s};
			EndIf
			
			ll+=1;
			Line Loop(ll) = {-lgur~{k},(lul+6*(k-1)+4):(lul+6*(k-1)+1),
							lgul~{k},lsul_[],-ltul~{k}~{1},-ltur~{k}~{1},
							-lsur_[]};
			s+=1;
			Plane Surface(s) = {s};
		EndIf


	EndFor

EndFor

// Air surface
Line Loop(1) = {n_groove/2+1,n_groove/2+2,-lcu~{n_groove/2},ltur~{n_groove/2}~{n_tapes_ext},
				lgur~{n_groove/2},(lul+6*(n_groove/2-1)+5),
				(lul+6*(n_groove/2-1)+6),-lc2,-lc1,
				-(ldl+6*(n_groove/2-1)+6),-(ldl+6*(n_groove/2-1)+5),
				-lgdr~{n_groove/2},-ltdr~{n_groove/2}~{n_tapes_ext},
				lcd~{n_groove/2}};
Plane Surface(1) = {1};


