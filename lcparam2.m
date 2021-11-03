(* rotinas para o calculo de parametros unitarios e da matriz de admitancia nodal
   de linhas de transmissao aereas e cabos SC *)
(* nov-2016  *)
(* Wolfram Language Raw Program *)
(* algumas funcoes de uso geral *)
(* rotinas para copiar o linspace e o logspace do matlab *)
logspac = Compile[{{i,_Real},{f,_Real},{np,_Integer}},
                 Table[10.0^(i + (x*(f - i))/(Floor[np] - 1)), {x, 0, np - 1}]
                 ]
     
linspac = Compile[{{i,_Real}, {f,_Real}, {np,_Integer}}, 
                 Table[i + (x*(f - i))/(Floor[np] - 1), {x, 0, np - 1}]
                 ]
                 
(* no caso de haver a eliminacao de blindagens e/ou armardura  *)        
(* eh o mesmo codigo para a eliminacao de cabos pararraios em LTs *)        
eliprc = Compile[{{m,_Complex,2},{nc,_Integer},{np,_Integer}},
	            Take[Inverse[m], nc - np, nc - np]
	            ]

(* nao sei qual a utilidade para sistemas de cabos mais deixei o codigo *)	            
(*  reduce bundle to equivalent phase conductors -- matrizes devem estar em formato de admitancia *)
elibndc = Compile[{{aux,_Complex,2},{nb,_Integer},{nf,_Integer}},
	              Table[Sum[aux[[i, j]], {i, nb*m - (nb - 1), nb*m}, {j, nb*n - (nb - 1), nb*n}],
                       {m, nf},{n, nf}]
                 ]
                 
                 
(* Elimina o Cruzamento de Autovetores *)
(* Metodo de Newton-Raphson *)
(* Elimina o Cruzamento de Autovetores *)
(* Metodo de Newton-Raphson \
*)
eigvNR[G_?MatrixQ, V0_, d0_, Tol_: 10^-12] := 
  Module[{A = G, nc, scale, v, d, vec, res, jac, i, itr, outd, outv},
   nc = Length[A];
   scale = Norm[A];
   A = A/scale;
   outd = ConstantArray[0., nc];
   outv = ConstantArray[0., {nc, nc}];
   itr = 0;
   jac = ConstantArray[0., {nc + 1, nc + 1}];
   Do[{v = V0[[1 ;;, i]];
     d = d0[[i]];
     res = Max[Abs[A.v - v*d]];
     While[(res > Tol) && (itr < 50000),
      itr = itr + 1;
      jac[[1 ;; nc, 1 ;; nc]] = A - d*IdentityMatrix[nc];
      jac[[1 ;; nc, nc + 1]] = -v;
      jac[[nc + 1, 1 ;; nc]] = 2*v;
      res = Join[A.v - v*d, {v.v - 1}];
      vec = LinearSolve[jac, res];
      res = Max[Abs[vec]];
      v = v - vec[[1 ;; nc]];
      d = d - vec[[nc + 1]];
      ];
     outd[[i]] = d*scale;
     outv[[1 ;;, i]] = Normalize[v];}, {i, 1, nc}];
   {outd, outv}
   ]
   
   
(* Rotacao da Matriz de Autovetores para minimizar parte Imaginaria *)
(* Complemento ao Metodo de Newton-Raphson *)
rot[s_?MatrixQ]:=Module[{S=s,Nc,SA,SB,scale,ang,scale1,col,scale2,err1,err2,numerator,denominator,aaa,bbb,ccc,ddd,eee,fff,j},
Nc=Length[S];
scale1=scale2=scale=ang=err1=err2=numerator=denominator=ConstantArray[0.,Nc];
SA=SB=ConstantArray[0.,{Nc,Nc}];

Table[

Table[numerator[[col]]=numerator[[col]]+Im[S[[j,col]]]*Re[S[[j,col]]];
denominator[[col]]=denominator[[col]]+Re[S[[j,col]]]^2-Im[S[[j,col]]]^2
,{j,Nc}];

numerator[[col]]=-2*numerator[[col]];
ang[[col]]=0.5*ArcTan[denominator[[col]],numerator[[col]]];
scale1[[col]]=Cos[ang[[col]]]+I*Sin[ang[[col]]];
scale2[[col]]=Cos[ang[[col]]+\[Pi]/2.]+I*Sin[ang[[col]]+\[Pi]/2.];
aaa=bbb=ccc=ddd=eee=fff=0.;

Table[SA[[j,col]]=S[[j,col]]*scale1[[col]];SB[[j,col]]=S[[j,col]]*scale2[[col]];
aaa=aaa+Im[SA[[j,col]]]^2;bbb=bbb+Re[SA[[j,col]]]*Im[SA[[j,col]]];ccc=ccc+Re[SA[[j,col]]]^2;
ddd=ddd+Im[SB[[j,col]]]^2;eee=eee+Re[SB[[j,col]]]*Im[SB[[j,col]]];fff=fff+Re[SB[[j,col]]]^2
,{j,Nc}];

err1[[col]]=aaa*Cos[ang[[col]]]^2+bbb*Sin[N[2*ang[[col]],16]]+ccc*Sin[ang[[col]]]^2;
err2[[col]]=ddd*Cos[ang[[col]]]^2+eee*Sin[N[2*ang[[col]],16]]+fff*Sin[ang[[col]]]^2;
If[err1[[col]]<err2[[col]],scale[[col]]=scale1[[col]],scale[[col]]=scale2[[col]]];
S[[All,col]]=S[[All,col]]*scale[[col]]

,{col,Nc}];

Return[S]
]


(* Produto Interno (Procedural) *)
(* Tentativa de traducao direta da versao em Matlab *)
xchEig[oldT_?MatrixQ,T_?MatrixQ,d1_?ListQ]:=
Module[{oldT0=oldT,T0=T,nc,ugh,ii,j,ilargest,rlargest,
	    dotprod,dot,ind,taken,d,hjelp,ihjelp,dum,dum2,l,eval, evect},

nc=Length[oldT0];
dot=ind=taken=hjelp=ConstantArray[0,nc];
d=DiagonalMatrix[d1];

ugh=Abs[Re[ConjugateTranspose[oldT0].T0]];

For[ii=1,ii<=nc,ii++,
ilargest=0;
rlargest=0;
For[j=1,j<=nc,j++,
dotprod=ugh[[ii,j]];
If[dotprod>rlargest,rlargest=Abs[Re[dotprod]];ilargest=j];
];
dot[[ii]]=rlargest;
ind[[ii]]=ii;
taken[[ii]]=0
];

(*sorting inner products in descending order*)
For[ii=1,ii<= nc,ii++,
For[j=1,j<=(nc-1),j++,
If[dot[[j]]    <dot[[j+1]],
    hjelp[[1]]=dot[[j+1]];
    ihjelp         =ind[[j+1]];
    dot[[j+1]]=dot[[j]];
    ind[[j+1]]=ind[[j]];
    dot[[j]]      =hjelp[[1]];
    ind[[j]]      =ihjelp
]
]
];

(*Doing the interchange in a prioritized sequence*)
For[l=1,l<=nc,l++,
ii=ind[[l]];
ilargest=0;
rlargest=0;

For[j=1,j<=nc,j++,
If[taken[[j]]==0,dotprod=ugh[[ii,j]];
If[dotprod>rlargest,
   rlargest=Abs[Re[dotprod]];
  ilargest=j
]
]
];

taken[[ii]]=1;

hjelp=T0[[1;;,ii]];
T0[[1;;,ii]]=T0[[1;;,ilargest]];
T0[[1;;,ilargest]]=hjelp;

hjelp=d[[ii,ii]];
d[[ii,ii]]=d[[ilargest,ilargest]];
d[[ilargest,ilargest]]=hjelp;

dum=ugh[[1;;,ii]];
ugh[[1;;,ii]]=ugh[[1;;,ilargest]];
ugh[[1;;,ilargest]]=dum];


dum2=DiagonalMatrix[Diagonal[Sign[Re[ConjugateTranspose[oldT0].T0]]]];
  
eval=Diagonal[d];
evect=T0.dum2;

{eval,evect}]


(* rotinas para evitar as oscilacoes de angulo de -Pi e +Pi *)

UnwrapPhase[data_?VectorQ, tol_: Pi, inc_: 2 Pi] := 
 FixedPoint[# + 
    inc*FoldList[Plus, 0., 
      Sign[Chop[ListCorrelate[{1, -1}, #], tol]   (*close Chop*)]
      (*close Sign*)] &,(*close FoldList*)data]
(*close FixedPoint and overall function*)

UnwrapPhase[list : {{_, _} ..}] := 
 Transpose[{list[[All, 1]], UnwrapPhase[list[[All, -1]]]}]


phase = Compile[{{l, _Complex, 1}}, 
   FoldList[
    Function[{prev, new}, # + Round[prev - #, 2 Pi] &@Arg@new], 
    Arg@First@l, Rest@l], CompilationTarget -> "C", 
   RuntimeOptions -> "Speed"]

(*                    *)
(* linhas aereas      *)
(*                    *)
(* impedancia interna de condutores tubulares *)
ZintTubo[Omega_, Rhoc_, rf_, rint_, Mur_: 1, Mu_: (4.*Pi)/10^7] := 
 With[{Etac = N[Sqrt[(I*Omega*Mur*Mu)/Rhoc]], ri = rint + 10^-6}, 
  With[{Den = 
     BesselK[1, Etac*ri]*BesselI[1, Etac*rf] - 
      BesselK[1, Etac*rf]*BesselI[1, Etac*ri], 
    Num = BesselK[1, Etac*ri]*BesselI[0, Etac*rf] + 
      BesselK[0, Etac*rf]*BesselI[1, Etac*ri]}, ((Rhoc*Etac)*
      Num)/((2*Pi*rf)*Den)]]
      
(* impedancia interna de condutor sem alma de aco *)
(* impedancia interna de condutores cilindricos *) 
(*Mur=90 para cables de aco*)    
Zin[Omega_, Rhopr_, rpr_, Mur_:90, Mu_:(4*Pi)/10^7] := 
    With[{Etapr = Sqrt[(I*Omega*Mu*Mur)/Rhopr]}, 
         ((Etapr*Rhopr)*BesselI[0, Etapr*rpr])/((2*Pi*rpr)*BesselI[1,Etapr*rpr])]
         
         (* internal impedance of cylindrical conductors *)
zintc= Compile[{{omega, _Complex}, {Rhoc, _Real}, {rf, _Real}, {ri,_Real}, {Mur, _Real}},
 ZintTubo[omega, Rhoc, rf, ri, Mur]]
 
zic= Compile[{{omega,_Complex},{rhoc, _Real},{rf, _Real}, {mur,_Real}},
	 Zin[omega, rhoc, rf, mur]]

(* external impedance of overhead lines *)
(* with ground wires *)
zextc = Compile[{{omega, _Complex}, {sigma, _Complex}, {npr,_Integer},
	             {x, _Real, 1}, {y, _Real, 1}, {rf, _Real}, {rpr, _Real}}, 
        Module[{Mu = 4.*10^-7*Pi, nc = Length[x], p, Deltaxij, yij, Deltayij},
        p = Sqrt[1.0/(I*omega*Mu*sigma)];
        I*omega*Mu/(2*Pi)*
        Table[
        	If[i != j, 
        	Deltaxij = (x[[i]] - x[[j]]);
        	yij = y[[i]] + y[[j]]; 
        	Deltayij = y[[i]] - y[[j]];
        	0.5*Log[(Deltaxij^2 + (2*p + yij)^2)/(Deltaxij^2 + (Deltayij)^2)],
        	If[i <= nc - npr,
        	Log[2.*(y[[i]] + p)/rf],
        	Log[2.*(y[[i]] + p)/rpr]]], {i, nc}, {j, nc}]]]
        	
(* 
  evaluate impedance and admittannce matrices per unit of length 
  using complex ground plane 
  as overloading in Mma does not work with compiled functions using uncompiled version
 *)  
 (* case 1: with ground wires and bundled conductors  *)
cZYlt[omega_,(x_)?VectorQ, (y_)?VectorQ, sigmas_, rdc_, rf_, rint_, npr_, rdcpr_, rpr_, nb_]:=
    Module[{mu=4*Pi*10^-7, eps=8.854*10^-12, nc=Length[x],  nf, rhoc, rhopr, mp, zin, ze, Z1, Y1},
    	
    nf= IntegerPart[(nc-npr)/nb];
    rhoc= rdc*Pi*(rf^2-rint^2);
    rhopr= rdcpr*Pi*rpr^2;
    
   zin= DiagonalMatrix[
   	    Join[zintc[omega, rhoc, rf, rint, 1] Table[1, {i, nc - npr}], 
             zic[omega, rhopr, rpr, 90] Table[1, {i, npr}]]
             ];
   ze= I omega mu/(2 Pi) With[{p = Sqrt[1/(I*omega*mu*sigmas)]}, 
         Table[If[i != j, 
         	  (Log[((x[[i]] - x[[j]])^2 + (2*p + y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 +(y[[i]] - y[[j]])^2)])/(2), 
              If[i <= nc - npr, 
               (Log[(2*(y[[i]] + p))/rf]), 
               (Log[(2*(y[[i]] + p))/rpr])]],
         {i, 1, nc}, 
         {j, 1, nc}]];
         
    Z1 = Inverse[elibndc[eliprc[zin+ze, nc, npr], nb, nf]];     
         
    mp= Table[If[i != j, 
   	    (1/2)*Log[((x[[i]] - x[[j]])^2 + (y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 + (y[[i]] - y[[j]])^2)], 
        If[i <= nc - npr, Log[(2*y[[i]])/rf], 
        Log[(2*y[[i]])/rpr]]],{i, 1, nc}, {j, 1, nc}];
 
 (* inclusao de condutancia shunt para evitar alguns problemas numericos *)         
    Y1= 3.0*10^-12 DiagonalMatrix[Table[1,{nf}]] +I omega*2*Pi*eps*elibndc[eliprc[mp, nc, npr], nb, nf];
          
    {Z1, Y1}                 

    ];
    
(* case  2: ground wires and unbundled conductors *)
cZYlt[omega_,(x_)?VectorQ, (y_)?VectorQ, sigmas_, rdc_,rf_,rint_, npr_, rdcpr_, rpr_]:=
    Module[{mu=4*Pi*10^-7, eps=8.854*10^-12, nc=Length[x],  nf, rhoc, rhopr, mp, zin, ze, Z1, Y1},
    	
    nf= IntegerPart[(nc-npr)];
    rhoc= rdc*Pi*(rf^2-rint^2);
    rhopr= rdcpr*Pi*rpr^2;
   
   If[ rint != 0, 
   zin= DiagonalMatrix[
   	    Join[zintc[omega, rhoc, rf, rint, 1] Table[1, {i, nc - npr}], 
             zic[omega, rhopr, rpr, 1] Table[1, {i, npr}]]
             ],
   zin= DiagonalMatrix[
   	    Join[zic[omega, rhoc, rf, 1] Table[1, {i, nc - npr}], 
             zic[omega, rhopr, rpr,1] Table[1, {i, npr}]]
             ];
   ];
   
             
   ze= I omega mu/(2 Pi) With[{p = Sqrt[1/(I*omega*mu*sigmas)]}, 
         Table[If[i != j, 
         	  (Log[((x[[i]] - x[[j]])^2 + (2*p + y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 +(y[[i]] - y[[j]])^2)])/(2), 
              If[i <= nc - npr, 
               (Log[(2*(y[[i]] + p))/rf]), 
               (Log[(2*(y[[i]] + p))/rpr])]],
         {i, 1, nc}, 
         {j, 1, nc}]];
         
    Z1 = Inverse[eliprc[zin+ze, nc, npr]];     
         
    mp= Table[If[i != j, 
   	    (1/2)*Log[((x[[i]] - x[[j]])^2 + (y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 + (y[[i]] - y[[j]])^2)], 
        If[i <= nc - npr, Log[(2*y[[i]])/rf], 
        Log[(2*y[[i]])/rpr]]],{i, 1, nc}, {j, 1, nc}];
          
    Y1=3.0*10^-12 DiagonalMatrix[Table[1,{nf}]] + I omega*2*Pi*eps* eliprc[mp, nc, npr];
          
    {Z1, Y1}       
];
     
     
(* case 3: no ground wires or bundled conductors *)   
cZYlt[omega_,(x_)?VectorQ, (y_)?VectorQ, sigmas_, rdc_,rf_,rint_]:=
    Module[{mu=4*Pi*10^-7, eps=8.854*10^-12, nc=Length[x],  nf, rhoc, mp, zin, ze, Z1, Y1},
    	
    nf= (nc);
    rhoc= rdc*Pi*(rf^2-rint^2);
    
   
   If[ rint != 0, 
   zin= DiagonalMatrix[zintc[omega, rhoc, rf, rint, 1] Table[1, {i, nc}]],
   zin= DiagonalMatrix[zic[omega, rhoc, rf, 1] Table[1, {i, nc}]]
   ];
   
             
   ze= I omega mu/(2 Pi) With[{p = Sqrt[1/(I*omega*mu*sigmas)]}, 
         Table[If[i != j, 
         	  (Log[((x[[i]] - x[[j]])^2 + (2*p + y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 +(y[[i]] - y[[j]])^2)])/(2), 
         	  (Log[(2*(y[[i]] + p))/rf])],
         {i, 1, nc}, 
         {j, 1, nc}]];
         
    Z1 =  zin+ze ;     
         
    mp= Table[If[i != j, 
   	    (1/2)*Log[((x[[i]] - x[[j]])^2 + (y[[i]] + y[[j]])^2)/((x[[i]] - x[[j]])^2 + (y[[i]] - y[[j]])^2)], 
   	    Log[(2*y[[i]])/rf]], {i, 1, nc}, {j, 1, nc}];
          
    Y1= 3.0*10^-12 DiagonalMatrix[Table[1,{nf}]] + I omega*2*Pi*eps*Inverse[mp]; 
        
    {Z1, Y1}      
];


(*                    *)
(* cabos subterraneos *)
(*                    *)

cZYci[omega_, r__, rhoc_, er1_, mur_: 1] := 
 Module[{r1 = r[[1]], r2 = r[[2]], z1, z2, y1, etac, 
   mu = 4.0*Pi*10^-7, eps = 8.854 10^-12},
   
  etac = Sqrt[(I*omega*mur*mu)/rhoc];
  
  z1 = etac rhoc/(2. Pi r1) BesselI[0, etac*r1]/BesselI[1, etac*r1];
  
  z2 = I omega mu 1./(2 \[Pi]) Log[r2/r1];
  
  y1 = I omega  2 \[Pi] (er1 eps)/Log[r2/r1];
  
  {z1 + z2, y1}
  ];

(* Z e Y cabo isolado  usando expressoes do Wedepohl *)
cZYciw[omega_, r__, rhoc_, er1_, mur_: 1] := 
 Module[{r1 = r[[1]], r2 = r[[2]], z1, z2, y1, etac, 
   mu = 4.0*Pi*10^-7, eps = 8.854 10^-12},
   
  etac = Sqrt[(I*omega*mur*mu)/rhoc];
  
  z1 = etac rhoc/(2. Pi r1) Coth[.7765*etac*r1] +.356 rhoc/(Pi*r1^2);
  
  z2 = I omega mu 1./(2 \[Pi]) Log[r2/r1];
  
  y1 = I omega  2 \[Pi] (er1 eps)/Log[r2/r1];
  
  {z1 + z2, y1}
  ] ;
  
(* Z e Y cabo com blindagem e isolacaoo *)
cZYcbi[omega_, r__, rhoc_, rhob_, er1_, er2_, mur2_: 1, mur1_: 1] := 
 Module[{r1 = r[[1]], r2 = r[[2]], r3 = r[[3]], r4 = r[[4]], Zi, Yi, 
   den, z1, z2, z3, z4, z5, z6, y1, y2, etac, etab, 
   mu = 4.0*Pi*10^-7, eps = 8.854 10^-12},
  etac = Sqrt[(I*omega*mur1*mu)/rhoc];
  etab = Sqrt[(I*omega*mur2*mu)/rhob];
  
  z1 = etac*rhoc/(2.*Pi*r1) BesselI[0, etac*r1]/BesselI[1, etac*r1];
  
  z2 = I*omega*mu*1./(2*Pi)* Log[r2/r1];
  
  den = BesselI[1, etab r3] BesselK[1, etab r2] - BesselI[1, etab  r2] BesselK[1, etab r3];
  
  z3 = (rhob etab)/(2 Pi r2) 1/den (BesselI[0, etab r2] BesselK[1, etab r3]
      + BesselK[0, etab  r2] BesselI[1, etab  r3]);
  
  z4 = rhob /(2 Pi r2 r3) 1/den;
  
  z5 = (rhob etab)/(2 Pi r3) 1/den (BesselI[0, etab r3] BesselK[1, etab r2] + BesselK[0, etab r3] BesselI[1, etab r2]);
  
  z6 = I*omega*mu/(2.*Pi) Log[r4/r3];
  
  Zi = {{z1 + z2 + z3 + z5 + z6 - 2*z4, z5 + z6 - z4},
  	    {z5 + z6 - z4, z5 + z6}};
  
  y1 = I*omega*2.*Pi*(er1 eps)/Log[r2/r1];
  
  y2 = I*omega*2.*Pi*(er1 eps)/Log[r4/r3];
   
   Yi = {{y1, -y1}, {-y1, y1 + y2}};
  
  
  	{Zi, Yi}
  	];


(* Line nodal admittance assembly *)
(* Calculo da matriz de admitancia nodal a partir de parametros unitarios e 
   comprimento do circuito *)
ynLT[Z_,Y_,length_]:=
     Module[{Z1=Z, Y1=Y, eval, evect, d, Tv, Tvi, hm, Am, Bm, y11, y12},
            {eval, evect} = Eigensystem[N[Z1.Y1]];
             d = Sqrt[eval];
             Tv = Transpose[evect];
             Tvi = Inverse[Transpose[evect]];
             hm = Exp[-d*length];
             Am = d*(1 + hm^2)/(1 - hm^2);
             Bm = -2.0*d*hm/(1 - hm^2);
             y11 = N[Inverse[Z1]].Tv.DiagonalMatrix[Am].Tvi;
             y12 = N[Inverse[Z1]].Tv.DiagonalMatrix[Bm].Tvi;
            {y11,y12} 
   
];
 	
  	


(* Z e Y cabo com blindagem e isolacaoo usando expressoes do Wedepohl *)
cZYcbiw[omega_, r__, rhoc_, rhob_, er1_, er2_, mur2_: 1, mur1_: 1] := 
 Module[{r1 = r[[1]], r2 = r[[2]], r3 = r[[3]], r4 = r[[4]], Zi, Yi, 
   z1, z2, z3, z4, z5, z6, y1, y2, etac, etab, delta, sum,
   mu = 4.0*Pi*10^-7, eps = 8.854 10^-12},
   delta = r3 - r2;
   sum = r2 + r3;
  etac = Sqrt[(I*omega*mur1*mu)/rhoc];
  etab = Sqrt[(I*omega*mur2*mu)/rhob];

(* condutor central *)  
  z1 = etac rhoc/(2. Pi r1) Coth[.7765*etac*r1] +.356 rhoc/(Pi*r1^2);
  z2 = I*omega*mu*1./(2*Pi)* Log[r2/r1];

(* blindagem *)  
  z3 = (rhob etab)/(2 Pi r2) Coth[etab* delta] - rhob/(2 Pi r2 sum);
  z4 = rhob etab/( Pi sum) Csch[etab*delta];
  z5 = (rhob etab)/(2 Pi r3) Coth[etab* delta] - rhob/(2 Pi r3 sum);
  z6 = I*omega*mu/(2.*Pi) Log[r4/r3];
  
  Zi = {{z1 + z2 + z3 + z5 + z6 - 2*z4, z5 + z6 - z4},
  	    {z5 + z6 - z4, z5 + z6}};
  
  y1 = I*omega*2.*Pi*(er1 eps)/Log[r2/r1];
  
  y2 = I*omega*2.*Pi*(er1 eps)/Log[r4/r3];
   
   Yi = {{y1, -y1}, {-y1, y1 + y2}};
  
  
  	{Zi, Yi}
  	];


(*  impedancia do solo usando Pollaczek *)    
 cZsolo[omega_, y1_, y2_, r_, s1_, er1_: 0.0] := 
 Module[{aux, d, di, eta, mu = 4.*Pi*10.0^-7, eps = 8.854 10^-12},
         eta = Sqrt[I omega mu (s1 + I omega er1 eps)];
         d = Sqrt[r^2 + (y1-y2)^2];
         di = Sqrt[r^2 + (y1+y2)^2];
  
  aux = NIntegrate[
    Exp[-(y1+y2)*Sqrt[lambda^2 + eta^2]]/
    (Abs[lambda] + Sqrt[lambda^2 + eta^2]) Exp[I r lambda],
    {lambda, -Infinity, Infinity}, 
     AccuracyGoal -> 9, PrecisionGoal -> 9, MaxRecursion -> 2000];
    
    I omega mu/(2*Pi)*(BesselK[0, eta*d] - BesselK[0, eta*di] +  aux)]   
    
(* impedancia de retorno de cabo isolado formula aproximada *)    
 cZsolo2[omega_, y1_, y2_, r_, s1_, er1_: 0.0] :=
    Module[{aux, d, di, eta, ell,mu = 4.*Pi*10.0^-7, eps = 8.854 10^-12},
         eta = Sqrt[I omega mu (s1 + I omega er1 eps)];
         ell = y1+y2;
         d = Sqrt[r^2 + (y1-y2)^2];
         di = Sqrt[r^2 + (y1+y2)^2];
         
          aux= (BesselK[0,eta*Sqrt[r^2+(y1-y2)^2] ]+
               ((y1+y2)^2-r^2)/(r^2+(y1+y2)^2)*(BesselK[2, eta Sqrt[r^2+(y1+y2)^2]]
               -2*( Exp[-(y1+y2)eta] (1+(y1+y2)eta ))/(eta^2 ( r^2+(y1+y2)^2)) ));
              
           I*omega*mu/(2*Pi)*aux    
        
         
        
         ];