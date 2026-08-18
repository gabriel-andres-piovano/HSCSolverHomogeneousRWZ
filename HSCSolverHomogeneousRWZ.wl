(* ::Package:: *)

(* ::Section:: *)
(*Begin package*)


BeginPackage["HSCSolverHomogeneousRWZ`"];


(* ::Subsubsection::Closed:: *)
(*Non-rescaled functions*)


RWsolverIn::usage = "RWsolverIn[l, \[Omega], r1g] gives the In solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r1g is the outermost integration radius.";


RWsolverUp::usage = "RWsolverUp[l, \[Omega], r2g] gives the Up solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r2g is the innermost integration radius.";


RWsolverUpNearHorizon::usage = "RWsolverUpNearHorizon[l, \[Omega], {\[Psi]up,d\[Psi]up}] provides the boundary condition near the horizon for the Up solution of the Regge-Wheeler equation."


ZsolverIn::usage = "RWsolverIn[l, \[Omega], r1g] gives the In solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r1g is the outermost integration radius.";


ZsolverUp::usage = "RWsolverUp[l, \[Omega], r2g] gives the Up solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r2g is the innermost integration radius.";


ZsolverUpNearHorizon::usage = "ZsolverUpNearHorizon[l, \[Omega], {\[Psi]up,d\[Psi]up}] provides the boundary condition near the horizon for the Up solution of the Zerilli equation."


(* ::Subsubsection::Closed:: *)
(*Rescaled functions*)


RWsolverInres::usage = "RWsolverInres[l, \[Omega], r1g] gives the rescaled In solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r1g is the outermost integration radius.";


RWsolverUpres::usage = "RWsolverUpres[l, \[Omega], r2g] gives the rescaled Up solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r2g is the innermost integration radius.";


RWsolverUpNearHorizonres::usage = "RWsolverUpNearHorizonres[l, \[Omega], {\[Psi]up,d\[Psi]up}] provides the rescaled boundary condition near the horizon for the Up solution of the Regge-Wheeler equation."


ZsolverInres::usage = "RWsolverInres[l, \[Omega], r1g] gives the rescaled In solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r1g is the outermost integration radius.";


ZsolverUpres::usage = "RWsolverUpres[l, \[Omega], r2g] gives the rescaled Up solution and its first-order derivative of the homogeneous Regge-Wheeler equation. r2g is the innermost integration radius.";


ZsolverUpNearHorizonres::usage = "ZsolverUpNearHorizonres[l, \[Omega], {\[Psi]up,d\[Psi]up}] provides the rescaled boundary condition near the horizon for the Up solution of the Zerilli equation."


(* ::Subsubsection::Closed:: *)
(*Start package*)


Begin["`Private`"];


(* ::Section:: *)
(*Functions for BCs and Regge-Wheeler solver in HS coordinates*)


(* ::Subsection::Closed:: *)
(*Boundary condition at the horizon - minus solution*)


bchorminRW[workingprecision_,l_,\[Omega]_]:=Module[{deltarp,A2n,rp,rin,err,errold,i,chor,p,q,phormin,qhor,dphormin,dqhor,\[Psi]hor,a2n,cInHor},
	rp=2;
	deltarp=rp/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor=4I*\[Omega];

	dqhor[n_]:=Which[n==0,0,n>0, (-1)^n/2^n (l+l^2-3 n)];
	dphormin[n_]:=Which[n==0,1-4I*\[Omega],n==1,-(1/2)-2I*\[Omega],n>1,(-1/2)^n];

	a2n[0]=1;
	A2n[n_]:= -(1/(n(n-chor)))Sum[(j*dphormin[n-j]+dqhor[n-j])a2n[j],{j,0,n-1}];

	err=1;
	errold=1;
	i=1;
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	phormin=p[rin,-1];
	qhor=q[rin];

	While[err > 10^(-workingprecision),
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=N[rp+deltarp,workingprecision];
			phormin=p[rin,-1];
			qhor=q[rin];
		];
		a2n[i]=Chop[A2n[i],10^(-workingprecision)];
		\[Psi]hor=Evaluate[1+Sum[a2n[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor''[rin]+phormin \[Psi]hor'[rin]+qhor \[Psi]hor[rin]];

		If[errold<=err&&errold> 10^(-workingprecision)&&i>5,
			Break[];
			,
			errold=err;
		];
		i++;
		If[i > 100, Break[]] (*Safeguard to avoid ruwaway computation*)   
	];
	cInHor=Table[a2n[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	{cInHor,rin}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at infinity - plus solution*)


bcinfplusRW[workingprecision_,l_,\[Omega]_]:=Module[{B1n,rp,rout,err,i,p,q,pinfplus,qinf,dpinfplus,dqinf,\[Psi]inf,b1n,cOutinf},
	rp=2;
	(*rout=2\[Pi] l(1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));*)
	rout=2\[Pi](1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));

	dqinf[n_]:=Which[n==0,0,n==1,0,n==2,-l(1+l),n>2,-(2^n/4)(-3+l+l^2)];
	dpinfplus[n_]:=Which[n==0,2I*\[Omega],n==1,4I*\[Omega],n>1, 1/2 (1+4I*\[Omega])2^n];
	
	err=1;
	i=1;
	b1n[0]=1;
	B1n[n_]:=(n-1)/(2I*\[Omega] ) b1n[n-1]+1/(2I*\[Omega]*n) Sum[(dqinf[j+1]-(n-j)dpinfplus[j])b1n[n-j],{j,1,n}];

	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=(6-l(1+l)r)/((-2+r)r^2);
	pinfplus=p[rout,1];
	qinf=q[rout];

	While[err > 10^(-workingprecision),
		If[Mod[i,15]==0,
			rout=2rout;
			pinfplus=p[rout,1];
			qinf=q[rout];
		];
        b1n[i]=B1n[i];
		\[Psi]inf=Evaluate[1+Sum[b1n[k](#)^(-k),{k,i}]]&;
		err=Abs[\[Psi]inf''[rout]+pinfplus \[Psi]inf'[rout]+qinf \[Psi]inf[rout]];
        i++;
		If[i > 150, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf=Table[b1n[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
	{cOutinf,rout}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at the horizon - plus solution*)


bchorplusRW[workingprecision_,l_,\[Omega]_]:=Module[{deltarp,A2n,rp,rin,err,errold,i,chor,p,q,phorplus,qhor,dphorplus,dqhor,\[Psi]hor,a2n,cInHor},
	rp=2;
	deltarp=rp/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor=-4I*\[Omega];

	dqhor[n_]:=Which[n==0,0,n>0,(-1)^n/2^n (l+l^2-3n)];
	dphorplus[n_]:=Which[n==0,1+4I*\[Omega],n==1,-(1/2)+2I*\[Omega],n>1,(-1)^n/2^n];

	a2n[0]=1;
	A2n[n_]:= -(1/(n(n-chor)))Sum[(j dphorplus[n-j]+dqhor[n-j])a2n[j],{j,0,n-1}];

	err=1;
	errold=1;
	i=1;
	p[r_,H_]:=2 (1+ I H r^2 \[Omega])/((-2+r)r);
	q[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	phorplus=p[rin,1];
	qhor=q[rin];

	While[err > 10^(-workingprecision),
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=N[rp+deltarp,workingprecision];
			phorplus=p[rin,-1];
			qhor=q[rin];
		];
        a2n[i]=Chop[A2n[i],10^(-workingprecision)];
		\[Psi]hor=Evaluate[1+Sum[a2n[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor''[rin]+phorplus \[Psi]hor'[rin]+qhor \[Psi]hor[rin]];

		If[errold<=err&&errold> 10^(-workingprecision)&&i>5,
			Break[];
			,
			errold=err;
		];
		i++;
		If[i > 100, Break[]] (*Safeguard to avoid ruwaway computation*)   
	];
	cInHor=Table[a2n[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	{cInHor,rin}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at infinity - minus solution*)


bcinfminRW[workingprecision_,l_,\[Omega]_]:=Module[{B1n,rp,rout,err,i,p,q,pinfmin,qinf,dpinfmin,dqinf,\[Psi]inf,b1n,cOutinf},
	rp=2;
	(*rout=2\[Pi] l(1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));*)
	rout=2\[Pi](1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));

	dqinf[n_]:=Which[n==0,0,n==1,0,n==2,-l(1+l),n>2,-(2^n/4)(-3+l+l^2)];
	 dpinfmin[n_]:=Which[n==0,-2I*\[Omega],n==1,-4I*\[Omega],n>1, 1/2 (1-4 I \[Omega])2^n];
	 
	err=1;
	i=1;
	b1n[0]=1;
	B1n[n_]:=-((n-1)/(2I*\[Omega] ))b1n[n-1]-1/(2I*\[Omega]*n) Sum[(dqinf[j+1]-(n-j)dpinfmin[j])b1n[n-j],{j,1,n}];

	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=(6-l(1+l)r)/((-2+r)r^2);
	pinfmin=p[rout,-1];
	qinf=q[rout];

	While[err > 10^(-workingprecision),
		If[Mod[i,15]==0,
			rout=2rout;
			pinfmin=p[rout,-1];
			qinf=q[rout];
		];
        b1n[i]=B1n[i];
		\[Psi]inf=Evaluate[1+Sum[b1n[k](#)^(-k),{k,i}]]&;
		err=Abs[\[Psi]inf''[rout]+pinfmin \[Psi]inf'[rout]+qinf \[Psi]inf[rout]];
        i++;
		If[i > 150, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf=Table[b1n[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
	{cOutinf,rout}
]


(* ::Subsection::Closed:: *)
(*Regge-Wheeler solver in hyperboloidal slicing coordinates*)


(* ::Subsubsection::Closed:: *)
(*In solution rescaled*)


RWsolverInres[l_,\[Omega]_,r1g_]:=Module[{workODE,precBC,precGoal,accGoal,rin,rp,p,qRW,cInH,\[Psi]hor,eqhor,X,Y,r,\[Psi]in,d\[Psi]in,nmaxhor},
	rp=2;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r1g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r1g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r1g]}];
	];
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	qRW[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	{cInH,rin}=bchorminRW[precBC,l,\[Omega]];
	nmaxhor=Length[cInH];
	\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;

	eqhor={
			X'[r]== Y[r],
			Y'[r]==-p[r,-1]Y[r]-qRW[r] X[r],
			X[rin]==\[Psi]hor[rin],Y[rin]== \[Psi]hor'[rin]
		};

	{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,r1g},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Function[{r},Evaluate[If[r<=rin,Evaluate[\[Psi]hor[r]],Evaluate[\[Psi]in[r]]]],Listable],Function[{r},Evaluate[If[r<=rin,Evaluate[\[Psi]hor'[r]],Evaluate[d\[Psi]in[r]]]],Listable],rin}
]


(* ::Subsubsection::Closed:: *)
(*Up solution rescaled*)


RWsolverUpres[l_,\[Omega]_,r2g_]:=Module[{workODE,precBC,precGoal,accGoal,rout,p,qRW,cOutinf,\[Psi]inf,eqinf,r,X,Y,\[Psi]up,d\[Psi]up,nmaxinf},
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r2g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r2g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r2g]}];
	];
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	qRW[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	{cOutinf,rout}=bcinfplusRW[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	eqinf={
			X'[r]== Y[r],
			Y'[r]==-p[r,1]Y[r]-qRW[r] X[r],
			X[rout]==\[Psi]inf[rout],Y[rout]==\[Psi]inf'[rout]
		};  

	{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,r2g,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];
	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Function[{r},Evaluate[If[r>=rout,Evaluate[\[Psi]inf[r]],Evaluate[\[Psi]up[r]]]],Listable],Function[{r},Evaluate[If[r>=rout,Evaluate[\[Psi]inf'[r]],Evaluate[d\[Psi]up[r]]]],Listable],rout}
]


(* ::Subsubsection::Closed:: *)
(*In solution*)


RWsolverIn[l_,\[Omega]_,r1g_]:=Module[{workODE,precBC,precGoal,accGoal,rin,rout,routplus,routmin,rp,p,qRW,cInH,\[Psi]hor,cOutinfplus,cOutinfmin,nmaxinfplus,nmaxinfmin,eqhor,r,rtor,X,Y,\[Psi]in,d\[Psi]in,nmaxhor,\[Psi]infplus,\[Psi]infmin,Bref,Binc,Rinfplus,dRinfplus,Rinfmin,dRinfmin,Rin,dRin,B1,B2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r1g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r1g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r1g]}];
	];
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	qRW[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	{cInH,rin}=bchorminRW[precBC,l,\[Omega]];
	nmaxhor=Length[cInH];
	\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;

	{cOutinfplus,routplus}=bcinfplusRW[precBC,l,\[Omega]];
	{cOutinfmin,routmin}=bcinfminRW[precBC,l,\[Omega]];
	nmaxinfplus=Length[cOutinfplus];
	nmaxinfmin=Length[cOutinfmin];
	rout=Max[{routplus,routmin}];

	\[Psi]infplus=Evaluate[Sum[cOutinfplus[[i]]#^(-i+1),{i,nmaxinfplus}]]&;
	\[Psi]infmin=Evaluate[Sum[cOutinfmin[[i]]#^(-i+1),{i,nmaxinfmin}]]&;

	Rinfplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]infplus[#]]&;
	dRinfplus=Evaluate[Rinfplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]infplus'[#]]&;
	Rinfmin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]infmin[#]]&;
	dRinfmin=Evaluate[Rinfmin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]infmin'[#]]&;

	eqhor={
			X'[r]== Y[r],
			Y'[r]==-p[r,-1]Y[r]-qRW[r] X[r],
			X[rin]==\[Psi]hor[rin],Y[rin]== \[Psi]hor'[rin]
		};

	If[r1g>=rout,
		{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

		{Bref,Binc}={B1,B2}/.NSolve[B1*Rinfplus[rout]+B2*Rinfmin[rout]==(Exp[-I*\[Omega]*rtor[rout]]\[Psi]in[rout])&&B1*dRinfplus[rout]+B2*dRinfmin[rout]==(Exp[-I*\[Omega]*rtor[rout]]((-((I*\[Omega]*rout)/(rout-rp)))\[Psi]in[rout]+d\[Psi]in[rout])),{B1,B2}][[1]];

		Rin=Function[{r},Evaluate[If[r<=rout,
									Evaluate[If[r>=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]in[r]],Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]hor[r]]]]
									,
									Evaluate[{Bref,Binc} . {Rinfplus[r],Rinfmin[r]}]
								]],Listable];
		dRin=Function[{r},Evaluate[If[r<=rout,
									Evaluate[If[r>=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]in[r]+d\[Psi]in[r])],Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]hor[r]+\[Psi]hor'[r])]]]
									,
									Evaluate[{Bref,Binc} . {dRinfplus[r],dRinfmin[r]}]
								]],Listable];
		,
		{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,r1g},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

		Rin=Function[{r},Evaluate[If[r<=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]hor[r]],Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]in[r]]]],Listable];
		dRin=Function[{r},Evaluate[If[r<=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]hor[r]+\[Psi]hor'[r])],Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]in[r]+d\[Psi]in[r])]]],Listable];
	];
	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Rin,dRin,{rin,rout}}
]


(* ::Subsubsection::Closed:: *)
(*Up solution*)


RWsolverUp[l_,\[Omega]_,r2g_]:=Module[{workODE,precBC,precGoal,accGoal,\[Lambda],rin,rout,rinplus,rinmin,rp,p,qRW,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,\[Psi]up,d\[Psi]up,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,\[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,Rup,dRup,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r2g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r2g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r2g]}];
	];
	\[Lambda]=1/2(-1+l)(2+l);
	p[r_,H_]:=2(1+I*H*r^2*\[Omega])/((-2+r)r);
	qRW[r_]:=(6-l(1+l)r)/((-2+r)r^2);

	{cOutinf,rout}=bcinfplusRW[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusRW[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminRW[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	eqinf={
			X'[r]== Y[r],
			Y'[r]==-p[r,1]Y[r]-qRW[r] X[r],
			X[rout]==\[Psi]inf[rout],Y[rout]==\[Psi]inf'[rout]
		};  

	If[r2g>=rout,
		Rup=Function[{r},Evaluate[Exp[I \[Omega] rtor[r]]\[Psi]inf[r]],Listable];
		dRup=Function[{r},Evaluate[Exp[I \[Omega] rtor[r]](((I \[Omega] r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])],Listable];
		,
		If[r2g<=rin,
			{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,rin,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];
			{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2 Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];
		
			Rup=Function[{r},Evaluate[If[r>=rin,
										Evaluate[If[r<=rout,Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]up[r]],Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]inf[r]]]]
										,
										Evaluate[{Cinc,Cref} . {Rhorplus[r],Rhormin[r]}]
									]],Listable];
			dRup=Function[{r},Evaluate[If[r>=rin,
										Evaluate[If[r<=rout,Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]up[r]+d\[Psi]up[r])],Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])]]]
										,
										Evaluate[{Cinc,Cref} . {dRhorplus[r],dRhormin[r]}]
									]],Listable];
			,
			{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,r2g,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

			Rup=Function[{r},Evaluate[If[r>rout,Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]inf[r]],Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]up[r]]]],Listable];
			dRup=Function[{r},Evaluate[If[r>rout,Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])],Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]up[r]+d\[Psi]up[r])]]],Listable];
		];
	];
	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Rup,dRup,{rin,rout}}
]


(* ::Subsubsection::Closed:: *)
(*Up solution near the horizon*)


RWsolverUpNearHorizon[l_,\[Omega]_,{\[Psi]up_,d\[Psi]up_}]:=Module[{precBC,rin,rout,rinplus,rinmin,rp,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,\[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;

	If[(Precision[\[Omega]]==MachinePrecision),
		precBC=13;
		,
		precBC=Precision[\[Omega]];
	];

	{cOutinf,rout}=bcinfplusRW[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusRW[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminRW[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2*Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];
	
	{
		{Function[{r},Evaluate[Cinc*Rhorplus[r]],Listable],Function[{r},Evaluate[Cinc*dRhorplus[r]],Listable]},
		{Function[{r},Evaluate[Cref*Rhormin[r]],Listable],Function[{r},Evaluate[Cref*dRhormin[r]],Listable]},
		{rin,rout},{Cinc,Cref}
	}
]


(* ::Subsubsection::Closed:: *)
(*Up solution near the horizon rescaled*)


RWsolverUpNearHorizonres[l_,\[Omega]_,{\[Psi]up_,d\[Psi]up_}]:=Module[
	{precBC,rin,rout,rinplus,rinmin,rp,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,
	 \[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,Rhorplusres,dRhorplusres,Rhorminres,dRhorminres,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;

	If[(Precision[\[Omega]]==MachinePrecision),
		precBC=13;
		,
		precBC=Precision[\[Omega]];
	];

	{cOutinf,rout}=bcinfplusRW[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusRW[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminRW[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2*Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];
	
	Rhorplusres=Evaluate[\[Psi]horplus[#]]&;
	dRhorplusres=Evaluate[\[Psi]horplus'[#]]&;
	Rhorminres=Evaluate[\[Psi]hormin[#]]&;
	dRhorminres=Evaluate[\[Psi]hormin'[#]]&;
	
	{
		{Function[{r},Evaluate[Cinc*Rhorplusres[r]],Listable],Function[{r},Evaluate[Cinc*dRhorplusres[r]],Listable]},
		{Function[{r},Evaluate[Cref*Rhorminres[r]],Listable],Function[{r},Evaluate[Cref*dRhorminres[r]],Listable]},
		{rin,rout},{Cinc,Cref}
	}
]


(* ::Section:: *)
(*Functions for BCs and Zerilli solver in HS coordinates*)


(* ::Subsection::Closed:: *)
(*Boundary condition at the horizon - minus solution*)


bchorminZ[workingprecision_,l_,\[Omega]_]:=Module[{deltarp,A2n,rp,rin,err,errold,i,chor,p,q,phormin,qhor,dphormin,dqhor,\[Psi]hor,a2n,cInHor,\[Lambda]},
	rp=2;
	deltarp=rp/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor=4I*\[Omega];
	\[Lambda]=1/2 (-1+l)(2+l);
	dqhor[n_]:=Which[n==0,0,n==1,(-3-4\[Lambda](1+\[Lambda]))/(6+4\[Lambda]),n>=2,2^-n ((-1)^n (2-3n+2\[Lambda])+2(1/(3+2\[Lambda]))^n (2^(1+n) (-\[Lambda])^n (1+\[Lambda])+(-3-2\[Lambda])^n (3n-2(1+\[Lambda]))))-(18(-1)^n (-1+n))/(2^n(3+2\[Lambda])^2) Hypergeometric2F1[2,2-n,1-n,(2\[Lambda])/(3+2\[Lambda])]];
	dphormin[n_]:=Which[n==0,1-4I*\[Omega],n==1,-(1/2)-2I*\[Omega],n>1,(-1)^n/2^n];

	a2n[0]=1;
	A2n[n_]:= -(1/(n(n-chor)))Sum[(j dphormin[n-j]+dqhor[n-j])a2n[j],{j,0,n-1}];

	err=1;
	errold=1;
	i=1;
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=-((2(9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));
	phormin=p[rin,-1];
	qhor=q[rin];

	While[err > 10^(-workingprecision),
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=N[rp+deltarp,workingprecision];
			phormin=p[rin,-1];
			qhor=q[rin];
		];
        a2n[i]=Chop[A2n[i],10^(-workingprecision)];
		\[Psi]hor=Evaluate[1+Sum[a2n[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor''[rin]+phormin \[Psi]hor'[rin]+qhor \[Psi]hor[rin]];

		If[errold<=err&&errold> 10^(-workingprecision)&&i>5,
			Break[];
			,
			errold=err;
		];
		i++;

		If[i > 100, Break[]] (*Safeguard to avoid ruwaway computation*)   
	];
	cInHor=Table[a2n[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	{cInHor,rin}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at infinity - plus solution*)


bcinfplusZ[workingprecision_,l_,\[Omega]_]:=Module[{M=1,B1n,rp,rout,err,i,p,q,pinfplus,qinf,dpinfplus,dqinf,\[Psi]inf,b1n,cOutinf,\[Lambda]},
	rp=2;
	(*rout=2\[Pi] l(1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));*)
	rout=2\[Pi](1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));

	\[Lambda]=1/2 (-1+l)(2+l);

	dqinf[n_]:=Which[n==0,0,n==1,0,n==2,-2(1+\[Lambda]),n==3,2+12/\[Lambda]-4\[Lambda],n>=4,((-\[Lambda])^-n (-27 2^n (-\[Lambda])^n+9 2^(2+n) (-\[Lambda])^(1+n)-12(-2 3^n+2 3^n*n+3 2^n (-\[Lambda])^n)\[Lambda]^2-16 3^n*n*\[Lambda]^3))/(36(3+2\[Lambda]))];
	dpinfplus[n_]:=Which[n==0,2I*\[Omega],n==1,4I*\[Omega],n>1,1/2 2^n(1+4I*\[Omega])];

	err=1;
	i=1;
	b1n[0]=1;
	B1n[n_]:=(n-1)/(2I*\[Omega] ) b1n[n-1]+1/(2I*\[Omega]*n) Sum[(dqinf[j+1]-(n-j)dpinfplus[j])b1n[n-j],{j,1,n}];

	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=-((2(9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));
	pinfplus=p[rout,1];
	qinf=q[rout];

	While[err > 10^(-workingprecision),
		If[Mod[i,15]==0,
			rout=2rout;
			pinfplus=p[rout,1];
			qinf=q[rout];
		];
        b1n[i]=B1n[i];
		\[Psi]inf=Evaluate[1+Sum[b1n[k](#)^-k,{k,i}]]&;
		err=Abs[\[Psi]inf''[rout]+pinfplus \[Psi]inf'[rout]+qinf \[Psi]inf[rout]];
        i++;
		If[i > 150, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf=Table[b1n[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
	{cOutinf,rout}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at the horizon - plus solution*)


bchorplusZ[workingprecision_,l_,\[Omega]_]:=Module[{deltarp,A2n,rp,rin,err,errold,i,chor,p,q,phorplus,qhor,dphorplus,dqhor,\[Psi]hor,a2n,cInHor,\[Lambda]},
	rp=2;
	deltarp=rp/50;
	rin=rp+deltarp; (*Closest r to the horizon *)
	chor=-4I*\[Omega];
	\[Lambda]=1/2(-1+l)(2+l);

	dqhor[n_]:=Which[n==0,0,n==1,(-3-4\[Lambda](1+\[Lambda]))/(6+4\[Lambda]),n>=2, 2^-n ((-1)^n (2-3n+2\[Lambda])+2(1/(3+2\[Lambda]))^n (2^(1+n) (-\[Lambda])^n (1+\[Lambda])+(-3-2\[Lambda])^n (3n-2(1+\[Lambda]))))-(18(-1)^n (-1+n))/(2^n(3+2\[Lambda])^2) Hypergeometric2F1[2,2-n,1-n,(2\[Lambda])/(3+2\[Lambda])]];
	dphorplus[n_]:=Which[n==0,1+4I*\[Omega],n==1,-(1/2)+2I*\[Omega],n>1,(-1)^n/2^n];

	a2n[0]=1;
	A2n[n_]:= -(1/(n(n-chor)))Sum[(j dphorplus[n-j]+dqhor[n-j])a2n[j],{j,0,n-1}];

	err=1;
	errold=1;
	i=1;
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=-((2 (9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));
	phorplus=p[rin,1];
	qhor=q[rin];

	While[err > 10^(-workingprecision),
		If[Mod[i,30]==0,
			deltarp=deltarp/2;
			rin=N[rp+deltarp,workingprecision];
			phorplus=p[rin,1];
			qhor=q[rin];
		];
        a2n[i]=Chop[A2n[i],10^(-workingprecision)];
		\[Psi]hor=Evaluate[1+Sum[a2n[k](#-rp)^k,{k,i}]]&;
		err=Abs[\[Psi]hor''[rin]+phorplus \[Psi]hor'[rin]+qhor \[Psi]hor[rin]];

		If[errold<=err&&errold> 10^(-workingprecision)&&i>5,
			Break[];
			,
			errold=err;
		];
		i++;
		If[i > 100, Break[]]  (*Safeguard to avoid ruwaway computation*)   
	];
	cInHor=Table[a2n[k],{k,0,i-1}]; (*Coefficients for ingoing waves near horizon (-)*)
	{cInHor,rin}
]


(* ::Subsection::Closed:: *)
(*Boundary condition at infinity - minus solution*)


bcinfminZ[workingprecision_,l_,\[Omega]_]:=Module[{M=1,B1n,rp,rout,err,i,p,q,pinfmin,qinf,dpinfmin,dqinf,\[Psi]inf,b1n,cOutinf,\[Lambda]},
	rp=2;
	(*rout=2\[Pi] l(1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));*)
	rout=2\[Pi](1/Abs[\[Omega]]+Abs[\[Omega]]/(1+Abs[\[Omega]]));

	\[Lambda]=1/2 (-1+l)(2+l);
	dqinf[n_]:=Which[n==0,0,n==1,0,n==2,-2(1+\[Lambda]),n==3,2+12/\[Lambda]-4\[Lambda],n>=4,((-\[Lambda])^-n (-27 2^n (-\[Lambda])^n+9 2^(2+n) (-\[Lambda])^(1+n)-12(-2 3^n+2 3^n n+3 2^n (-\[Lambda])^n)\[Lambda]^2-16 3^n n*\[Lambda]^3))/(36(3+2\[Lambda]))];
	 dpinfmin[n_]:=Which[n==0,-2I*\[Omega],n==1,-4I*\[Omega],n>1, 1/2 2^n(1-4I*\[Omega])];

	err=1;
	i=1;
	b1n[0]=1;
	B1n[n_]:=-((n-1)/(2I*\[Omega] ))b1n[n-1]-1/(2I*\[Omega]*n) Sum[(dqinf[j+1]-(n-j)dpinfmin[j])b1n[n-j],{j,1,n}];

	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	q[r_]:=-((2 (9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));
	pinfmin=p[rout,-1];
	qinf=q[rout];

	While[err > 10^(-workingprecision),
		If[Mod[i,15]==0,
			rout=2rout;
			pinfmin=p[rout,-1];
			qinf=q[rout];
		];
        b1n[i]=B1n[i];
		\[Psi]inf=Evaluate[1+Sum[b1n[k](#)^-k,{k,i}]]&;
		err=Abs[\[Psi]inf''[rout]+pinfmin \[Psi]inf'[rout]+qinf \[Psi]inf[rout]];
        i++;
		If[i > 150, Break[]]  (*Asymptotic expansions are not convergent*)    
	];
	cOutinf=Table[b1n[k],{k,0,i-1}]; (*Coefficients for outgoing waves at \[Infinity] (+)*)
	{cOutinf,rout}
]


(* ::Subsection::Closed:: *)
(*Zerilli solver in hyperboloidal slicing coordinates*)


(* ::Subsubsection::Closed:: *)
(*In solution rescaled*)


ZsolverInres[l_,\[Omega]_,r1g_]:=Module[{workODE,precBC,precGoal,accGoal,\[Lambda],rin,rp,p,qZ,cInH,\[Psi]hor,eqhor,r,X,Y,\[Psi]in,d\[Psi]in,nmaxhor},
	rp=2;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r1g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r1g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r1g]}];
	];
	\[Lambda]=1/2(-1+l)(2+l);
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	qZ[r_]:=-((2 (9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));

	{cInH,rin}=bchorminZ[precBC,l,\[Omega]];
	nmaxhor=Length[cInH];
	\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;

	eqhor={
			X'[r]== Y[r],
			Y'[r]==-p[r,-1]Y[r]-qZ[r] X[r],
			X[rin]==\[Psi]hor[rin],Y[rin]== \[Psi]hor'[rin]
		};

	{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,r1g},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Function[{r},Evaluate[If[r<=rin,Evaluate[\[Psi]hor[r]],Evaluate[\[Psi]in[r]]]],Listable],Function[{r},Evaluate[If[r<=rin,Evaluate[\[Psi]hor'[r]],Evaluate[d\[Psi]in[r]]]],Listable],rin}
]


(* ::Subsubsection::Closed:: *)
(*Up solution rescaled*)


ZsolverUpres[l_,\[Omega]_,r2g_]:=Module[{workODE,precBC,precGoal,accGoal,\[Lambda],rout,p,qZ,cOutinf,\[Psi]inf,eqinf,r,X,Y,\[Psi]up,d\[Psi]up,nmaxinf},
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r2g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r2g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r2g]}];
	];
	\[Lambda]=1/2(-1+l)(2+l);
	p[r_,H_]:=2 (1+I*H*r^2*\[Omega])/((-2+r)r);
	qZ[r_]:=-((2 (9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));

	{cOutinf,rout}=bcinfplusZ[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	eqinf={
			X'[r]== Y[r],
			Y'[r]==-p[r,1]Y[r]-qZ[r] X[r],
			X[rout]==\[Psi]inf[rout],Y[rout]==\[Psi]inf'[rout]
		};  

	{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,r2g,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Function[{r},Evaluate[If[r>=rout,Evaluate[\[Psi]inf[r]],Evaluate[\[Psi]up[r]]]],Listable],Function[{r},Evaluate[If[r>=rout,Evaluate[\[Psi]inf'[r]],Evaluate[d\[Psi]up[r]]]],Listable],rout}
]


(* ::Subsubsection::Closed:: *)
(*In solution*)


ZsolverIn[l_,\[Omega]_,r1g_]:=Module[{workODE,precBC,precGoal,accGoal,\[Lambda],rin,rout,routplus,routmin,rp,p,qZ,cInH,\[Psi]hor,cOutinfplus,cOutinfmin,nmaxinfplus,nmaxinfmin,eqhor,r,rtor,X,Y,\[Psi]in,d\[Psi]in,nmaxhor,\[Psi]infplus,\[Psi]infmin,Bref,Binc,Rinfplus,dRinfplus,Rinfmin,dRinfmin,Rin,dRin,B1,B2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r1g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r1g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r1g]}];
	];
	\[Lambda]=1/2(-1+l)(2+l);
	p[r_,H_]:=2(1+I*H*r^2*\[Omega])/((-2+r)r);
	qZ[r_]:=-((2(9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));

	{cInH,rin}=bchorminZ[precBC,l,\[Omega]];
	nmaxhor=Length[cInH];
	\[Psi]hor=Evaluate[Sum[cInH[[i]](#-rp)^(i-1),{i,nmaxhor}]]&;

	{cOutinfplus,routplus}=bcinfplusZ[precBC,l,\[Omega]];
	{cOutinfmin,routmin}=bcinfminZ[precBC,l,\[Omega]];
	nmaxinfplus=Length[cOutinfplus];
	nmaxinfmin=Length[cOutinfmin];
	rout=Max[{routplus,routmin}];

	\[Psi]infplus=Evaluate[Sum[cOutinfplus[[i]]#^(-i+1),{i,nmaxinfplus}]]&;
	\[Psi]infmin=Evaluate[Sum[cOutinfmin[[i]]#^(-i+1),{i,nmaxinfmin}]]&;

	Rinfplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]infplus[#]]&;
	dRinfplus=Evaluate[Rinfplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]infplus'[#]]&;
	Rinfmin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]infmin[#]]&;
	dRinfmin=Evaluate[Rinfmin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]infmin'[#]]&;

	eqhor={
			X'[r]== Y[r],
			Y'[r]==-p[r,-1]Y[r]-qZ[r] X[r],
			X[rin]==\[Psi]hor[rin],Y[rin]== \[Psi]hor'[rin]
		};

	If[r1g>=rout,
		{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];
		{Bref,Binc}={B1,B2}/.NSolve[B1*Rinfplus[rout]+B2*Rinfmin[rout]==(Exp[-I*\[Omega]*rtor[rout]]\[Psi]in[rout])&&B1 dRinfplus[rout]+B2*dRinfmin[rout]==(Exp[-I*\[Omega]*rtor[rout]]((-((I*\[Omega]*rout)/(rout-rp)))\[Psi]in[rout]+d\[Psi]in[rout])),{B1,B2}][[1]];

		Rin=Function[{r},Evaluate[If[r<=rout,
									Evaluate[If[r>=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]in[r]],Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]hor[r]]]]
									,
									Evaluate[{Bref,Binc} . {Rinfplus[r],Rinfmin[r]}]
								]],Listable];
		dRin=Function[{r},Evaluate[If[r<=rout,
									Evaluate[If[r>=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]in[r]+d\[Psi]in[r])],Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]hor[r]+\[Psi]hor'[r])]]]
									,
									Evaluate[{Bref,Binc} . {dRinfplus[r],dRinfmin[r]}]
								]],Listable];
		,
		{\[Psi]in,d\[Psi]in}={X,Y}/.First@NDSolve[eqhor,{X,Y},{r,rin,r1g},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

		Rin=Function[{r},Evaluate[If[r<=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]hor[r]],Evaluate[Exp[-I*\[Omega]*rtor[r]]\[Psi]in[r]]]],Listable];
		dRin=Function[{r},Evaluate[If[r<=rin,Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]hor[r]+\[Psi]hor'[r])],Evaluate[Exp[-I*\[Omega]*rtor[r]]((-((I*\[Omega]*r)/(r-rp)))\[Psi]in[r]+d\[Psi]in[r])]]],Listable];
	];

	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Rin,dRin,{rin,rout}}
]


(* ::Subsubsection::Closed:: *)
(*Up solution*)


ZsolverUp[l_,\[Omega]_,r2g_]:=Module[{workODE,precBC,precGoal,accGoal,\[Lambda],rin,rout,rinplus,rinmin,rp,p,qZ,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,\[Psi]up,d\[Psi]up,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,\[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,Rup,dRup,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;
	If[(Precision[\[Omega]]==MachinePrecision)||(Precision[r2g]==MachinePrecision),
		workODE=MachinePrecision;
		precGoal=13;
		accGoal=13;
		precBC=13;
		,
		workODE=Min[{Precision[\[Omega]]-5,Precision[r2g]-5}];
		precGoal=workODE-5;
		accGoal=workODE-5;
		precBC=Min[{Precision[\[Omega]],Precision[r2g]}];
	];
	\[Lambda]=1/2(-1+l)(2+l);
	p[r_,H_]:=2(1+I*H*r^2*\[Omega])/((-2+r)r);
	qZ[r_]:=-((2(9+r*\[Lambda](9+r*\[Lambda](3+r+r*\[Lambda]))))/((-2+r)r^2(3+r*\[Lambda])^2));

	{cOutinf,rout}=bcinfplusZ[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusZ[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminZ[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	eqinf={
			X'[r]== Y[r],
			Y'[r]==-p[r,1]Y[r]-qZ[r] X[r],
			X[rout]==\[Psi]inf[rout],Y[rout]==\[Psi]inf'[rout]
		}; 
		 
	If[r2g>=rout,
		Rup=Function[{r},Evaluate[Exp[I \[Omega] rtor[r]]\[Psi]inf[r]],Listable];
		dRup=Function[{r},Evaluate[Exp[I \[Omega] rtor[r]](((I \[Omega] r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])],Listable];
		,
		If[r2g<=rin,
			{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,rin,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];
			{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2 Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];
		
			Rup=Function[{r},Evaluate[If[r>=rin,
										Evaluate[If[r<=rout,Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]up[r]],Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]inf[r]]]]
										,
										Evaluate[{Cinc,Cref} . {Rhorplus[r],Rhormin[r]}]
									]],Listable];
			dRup=Function[{r},Evaluate[If[r>=rin,
										Evaluate[If[r<=rout,Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]up[r]+d\[Psi]up[r])],Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])]]]
										,
										Evaluate[{Cinc,Cref} . {dRhorplus[r],dRhormin[r]}]
									]],Listable];
			,
			{\[Psi]up,d\[Psi]up}={X,Y}/.First@NDSolve[eqinf,{X,Y},{r,r2g,rout},Method->"StiffnessSwitching",WorkingPrecision->workODE,PrecisionGoal->precGoal,AccuracyGoal->accGoal,InterpolationOrder->All];

			Rup=Function[{r},Evaluate[If[r>rout,Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]inf[r]],Evaluate[Exp[I*\[Omega]*rtor[r]]\[Psi]up[r]]]],Listable];
			dRup=Function[{r},Evaluate[If[r>rout,Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]inf[r]+\[Psi]inf'[r])],Evaluate[Exp[I*\[Omega]*rtor[r]](((I*\[Omega]*r)/(r-rp))\[Psi]up[r]+d\[Psi]up[r])]]],Listable];
		];
	]

	(* Remove local variables not garbage collected*)
	ClearSystemCache[];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`X$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[Evaluate[ToExpression[Pick[Names["HSCSolverHomogeneousRWZ`Private`*"],StringMatchQ[#,"HSCSolverHomogeneousRWZ`Private`Y$"~~__]&/@Names["HSCSolverHomogeneousRWZ`Private`*"],True]]]];
	Remove[r];

	{Rup,dRup,{rin,rout}}
]


(* ::Subsubsection::Closed:: *)
(*Up solution near the horizon*)


ZsolverUpNearHorizon[l_,\[Omega]_,{\[Psi]up_,d\[Psi]up_}]:=Module[{precBC,rin,rout,rinplus,rinmin,rp,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,\[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;

	If[(Precision[\[Omega]]==MachinePrecision),
		precBC=13;
		,
		precBC=Precision[\[Omega]];
	];

	{cOutinf,rout}=bcinfplusZ[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusZ[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminZ[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2*Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];

	{
		{Function[{r},Evaluate[Cinc*Rhorplus[r]],Listable],Function[{r},Evaluate[Cinc*dRhorplus[r]],Listable]},
		{Function[{r},Evaluate[Cref*Rhormin[r]],Listable],Function[{r},Evaluate[Cref*dRhormin[r]],Listable]},
		{rin,rout},{Cinc,Cref}
	}
]


(* ::Subsubsection::Closed:: *)
(*Up solution near the horizon rescaled*)


ZsolverUpNearHorizonres[l_,\[Omega]_,{\[Psi]up_,d\[Psi]up_}]:=Module[
	{precBC,rin,rout,rinplus,rinmin,rp,\[Psi]inf,cOutinf,nmaxinf,eqinf,r,rtor,X,Y,cInHplus,cInHmin,nmaxhorplus,nmaxhormin,
	 \[Psi]horplus,\[Psi]hormin,Cref,Cinc,Rhorplus,dRhorplus,Rhormin,dRhormin,Rhorplusres,dRhorplusres,Rhorminres,dRhorminres,C1,C2},
	rp=2;
	rtor[r_]:=2Log[(r-rp)/2]+r;

	If[(Precision[\[Omega]]==MachinePrecision),
		precBC=13;
		,
		precBC=Precision[\[Omega]];
	];

	{cOutinf,rout}=bcinfplusZ[precBC,l,\[Omega]];
	nmaxinf=Length[cOutinf];
	\[Psi]inf=Evaluate[Sum[cOutinf[[i]]#^(-i+1),{i,nmaxinf}]]&;

	{cInHplus,rinplus}=bchorplusZ[precBC,l,\[Omega]];
	{cInHmin,rinmin}=bchorminZ[precBC,l,\[Omega]];
	nmaxhorplus=Length[cInHplus];
	nmaxhormin=Length[cInHmin];
	rin=Max[{rinplus,rinmin}];

	\[Psi]horplus=Evaluate[Sum[cInHplus[[i]](#-rp)^(i-1),{i,nmaxhorplus}]]&;
	\[Psi]hormin=Evaluate[Sum[cInHmin[[i]](#-rp)^(i-1),{i,nmaxhormin}]]&;

	Rhorplus=Evaluate[Exp[I*\[Omega]*rtor[#]]\[Psi]horplus[#]]&;
	dRhorplus=Evaluate[Rhorplus[#]((I*\[Omega]*#)/(#-rp))+Exp[I*\[Omega]*rtor[#]]\[Psi]horplus'[#]]&;
	Rhormin=Evaluate[Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin[#]]&;
	dRhormin=Evaluate[Rhormin[#](-((I*\[Omega]*#)/(#-rp)))+Exp[-I*\[Omega]*rtor[#]]\[Psi]hormin'[#]]&;

	{Cinc,Cref}={C1,C2}/.NSolve[C1*Rhorplus[rin]+C2*Rhormin[rin]==(Exp[I*\[Omega]*rtor[rin]]\[Psi]up[rin])&&C1*dRhorplus[rin]+C2*dRhormin[rin]==(Exp[I*\[Omega]*rtor[rin]](((I*\[Omega]*rin)/(rin-rp))\[Psi]up[rin]+d\[Psi]up[rin])),{C1,C2}][[1]];
	
	Rhorplusres=Evaluate[\[Psi]horplus[#]]&;
	dRhorplusres=Evaluate[\[Psi]horplus'[#]]&;
	Rhorminres=Evaluate[\[Psi]hormin[#]]&;
	dRhorminres=Evaluate[\[Psi]hormin'[#]]&;
	
	{
		{Function[{r},Evaluate[Cinc*Rhorplusres[r]],Listable],Function[{r},Evaluate[Cinc*dRhorplusres[r]],Listable]},
		{Function[{r},Evaluate[Cref*Rhorminres[r]],Listable],Function[{r},Evaluate[Cref*dRhorminres[r]],Listable]},
		{rin,rout},{Cinc,Cref}
	}
]


(* ::Section::Closed:: *)
(*End package*)


End[];


EndPackage[];
