#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
(*Assignment-06*)
(*Name: Santanu Das*)
(*Roll No.-PH22C040*)
(**)
(*Q1.*)


Clear["Global`*"]
solve=DSolve[{0.5*y''[x]+4*y[x]==x^2,y[0]+5*y'[1]==0,y[1]==0},y[x],x];(*Robin Boundary condition*)
plot1=Plot[y[x]/.solve,{x,0,1},PlotStyle->Black]


eq[x_]:=0.5*y''[x]+4*y[x]== f[x];
f[x_]:= x^2;
y[0]+5*y'[1]==0; y[1]=0; 
n=51;(*no of points*)
h=1/(n-1);
b1=Table[f[x],{x,0,1-h,h}];
b1=Join[b1,{0}];
di1[i_]:= Expand[(0.5*(Subscript[y, (i/h)+1] -2Subscript[y, (i/h)]+ Subscript[y, (i/h)-1]))/h^2+4*Subscript[y, i/h]] == f[i] 


unkn= Table[di1[i],{i,h,1-h,h}];
bc1=Subscript[y, 0]+5*((Subscript[y, (i/h)]-Subscript[y, (i/h)-1])/h)==0/.i->1 ;

mat1 = Join[{bc1},unkn,{Subscript[y, 1/h]==0}];
fun1=mat1[[2]]
fun2=mat1[[1]];

lhs1=fun1[[1]]
lhs2=fun2[[1]];
coef1={Coefficient[lhs1,Subscript[y, 0]],Coefficient[lhs1,Subscript[y, 1]],Coefficient[lhs1,Subscript[y, 2]]};
coef2={Coefficient[lhs2,Subscript[y, 0]],Coefficient[lhs2,Subscript[y, (1/h)-1]],Coefficient[lhs2,Subscript[y, 1/h]]};
MatrixForm@mat1;


matspa=SparseArray[{Band[{1,1}]->coef1[[2]], Band[{1,2}]->coef1[[3]], Band[{2,1}]->coef1[[1]]},{(1/h)+1,(1/h)+1}];
matspa;
A1=ReplacePart[matspa,{{1,1}->coef2[[1]],{1,2}->0,{1,1/h}->coef2[[2]],{1,(1/h)+1}->coef2[[3]],{(1/h)+1,1/h}->0,{(1/h)+1,(1/h)+1}->1}];
a=ReplacePart[matspa,{{1,1}->coef2[[1]],{1,2}->0,{1,1/h}->coef2[[2]],{1,(1/h)+1}->coef2[[3]],{(1/h)+1,1/h}->0,{(1/h)+1,(1/h)+1}->1}]//MatrixForm;
b2=b1//MatrixForm;
yfitval=LinearSolve[A1,b1]//N;
xvals=Table[x,{x,0,1,h}];
solx=Thread[List[xvals,yfitval]];
plot2=ListPlot[solx]
Show[plot1,plot2]


(*Local trucation error*)
ybudin=Table[y[x]/.solve,{x,0,1,h}]//Flatten;
yfitval;
lotruerr=Abs[ybudin-yfitval]


(*Global trucation error*)
glotruerr=Max[lotruerr](*it is defined as maximum of local trucation error*)


(*order of the error*)
order=Log[glotruerr]/Log[h]//N


(* ::Subsection:: *)
(*Q2. *)


(* ::Subsection:: *)
(*Here we are Defining a 2nd order non-linear differential equation "van der pol oscillator"*)
(* and the equation is given below*)


Clear["Global`*"]
Nlieqn=x''[t]-\[Mu]*(1-x[t]^2)*x'[t]+x[t]==0;

sol=NDSolve[{Nlieqn, x[0] == 0, x'[0] ==9 }, x[t], {t, 0, 20}]/.\[Mu]->1/30;
Plot[x[t]/. sol[[1]], {t, 0, 20}, PlotRange -> All, PlotStyle -> Blue]
X[t_]=x[t]/. sol[[1]];
ParametricPlot[{X[t],X'[t]},{t,0,20},PlotRange->All]
Animate[ParametricPlot[{X[t],X'[t]},{t,0,T},PlotRange->All],{T,0,20}]




