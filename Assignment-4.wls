#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
(*Assignment-04*)
(*Name: Santanu Das*)
(*Roll No.-PH22C040*)


(*(1).we define a transfer function*)
Clear["Global`*"]
f[s_]:=1/(I*s+1)
(*(2)12 data points from the transfer function*)
datax=Table[s,{s,-5,6,1}]
datay=Table[Re[f[s]],{s,-5,6,1}]
datapoints=Table[{s,Re[f[s]]},{s,-5,6,1}]
scatterplot=ListPlot[datapoints,PlotStyle->Black]
realf=Plot[Re[f[s]],{s,-5,6}]


(*(3)here we are using a 11th order polynomial to fit the data points*)
G[x_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(l = 1\), \(12\)]\(datay[\([\)\(l\)\(]\)]*\(
\*UnderoverscriptBox[\(\[Product]\), \(j = 1\), \(12\)]\((If[l == j, 1, \((x - datax[\([\)\(j\)\(]\)])\)/\((datax[\([\)\(l\)\(]\)] - datax[\([\)\(j\)\(]\)])\)])\)\)\)\)
G[x]//Expand
p2=Plot[G[x],{x,-5,6},AxesLabel->{"x","g[x]"},PlotStyle->Red,PlotRange->All,PlotLegends->{"lagrange interpolation"}]
com1=Show[scatterplot,p2,PlotRange->All]


(*(4)interpolation using cubic splines*)
(*Real Part*)
Clear["Global`*"]
f[s_]:=1/(s*I+1)
datax=Table[s,{s,-5,6,1}]
datay=Table[Re[f[s]],{s,-5,6,1}]
g[A_,B_,C_,D_,i_,x_]:=A*(x-datax[[i]])^3+B*(x-datax[[i]])^2+C*(x-datax[[i]])+D
g1[A_,B_,C_,i_,x_]:=3*A*(x-datax[[i]])^2+2*B*(x-datax[[i]])+C
g2[A_,B_,i_,x_]:=6*A*(x-datax[[i]])+2*B
A1=Table[a[r],{r,0,10,1}]
B1=Table[b[r],{r,0,10,1}]
C1=Table[c[r],{r,0,10,1}]
D1=Table[d[r],{r,0,10,1}]

splines=Table[g[A1[[t]],B1[[t]],C1[[t]],D1[[t]],t,x],{t,1,11,1}]
constraints1={Table[datay[[p]]==g[A1[[1]],B1[[1]],C1[[1]],D1[[1]],1,datax[[p]]],{p,1,2}],
Table[datay[[p]]==g[A1[[2]],B1[[2]],C1[[2]],D1[[2]],2,datax[[p]]],{p,2,3}],Table[datay[[p]]==g[A1[[3]],B1[[3]],C1[[3]],D1[[3]],3,datax[[p]]],{p,3,4}],
Table[datay[[p]]==g[A1[[4]],B1[[4]],C1[[4]],D1[[4]],4,datax[[p]]],{p,4,5}],Table[datay[[p]]==g[A1[[5]],B1[[5]],C1[[5]],D1[[5]],5,datax[[p]]],{p,5,6}],
Table[datay[[p]]==g[A1[[6]],B1[[6]],C1[[6]],D1[[6]],6,datax[[p]]],{p,6,7}],Table[datay[[p]]==g[A1[[7]],B1[[7]],C1[[7]],D1[[7]],7,datax[[p]]],{p,7,8}],
Table[datay[[p]]==g[A1[[8]],B1[[8]],C1[[8]],D1[[8]],8,datax[[p]]],{p,8,9}],Table[datay[[p]]==g[A1[[9]],B1[[9]],C1[[9]],D1[[9]],9,datax[[p]]],{p,9,10}],
Table[datay[[p]]==g[A1[[10]],B1[[10]],C1[[10]],D1[[10]],10,datax[[p]]],{p,10,11}],Table[datay[[p]]==g[A1[[11]],B1[[11]],C1[[11]],D1[[11]],11,datax[[p]]],{p,11,12}]};

constraints1=Flatten[constraints1,1];
constraints2=Table[g1[A1[[t]],B1[[t]],C1[[t]],t,datax[[t+1]]]==g1[A1[[t+1]],B1[[t+1]],C1[[t+1]],t+1,datax[[t+1]]],{t,1,10,1}];
AppendTo[constraints2,g1[A1[[1]],B1[[1]],C1[[1]],1,datax[[1]]]==0];
AppendTo[constraints2,g1[A1[[11]],B1[[11]],C1[[11]],11,datax[[12]]]==0];
constraints3=Table[g2[A1[[t]],B1[[t]],t,datax[[t+1]]]==g2[A1[[t+1]],B1[[t+1]],t+1,datax[[t+1]]],{t,1,10,1}];

constraint=Join[constraints1,constraints2,constraints3]
matrix3 = MatrixForm[constraint]
Var={a[0],a[1],a[2],a[3],a[4],a[5],a[6],a[7],a[8],a[9],a[10],b[0],b[1],
b[2],b[3],b[4],b[5],b[6],b[7],b[8],b[9],b[10],c[0],c[1],c[2],c[3],c[4],
c[5],c[6],c[7],c[8],c[9],c[10],d[0],d[1],d[2],d[3],d[4],d[5],d[6],d[7],d[8],d[9],d[10]}
sol=Solve[constraint,Var]


(*defining splines*)
A2={a[0]=-(13813923/1358149195),a[1]=244828807/23088536315,a[2]=133094789/23088536315,a[3]=496470021/4617707263,a[4]=-(774038093/1776041255),a[5]=17620643/40435265,a[6]=-(495361965/4617707263),a[7]=-(154147853/23088536315),a[8]=-(166156831/23088536315),a[9]=-(58798149/23088536315),a[10]=126895993/23088536315};
B2={b[0]=1409930517/46177072630,b[1]=910371/46177072630,b[2]=1469883213/46177072630,b[3]=2268451947/46177072630,b[4]=17162552577/46177072630,b[5]=-(2541906981/2716298390),b[6]=17155904241/46177072630,b[7]=2295045291/46177072630,b[8]=1370158173/46177072630,b[9]=10086951/1248028990,b[10]=20428293/46177072630};
C2={c[0]=0,c[1]=705420444/23088536315,c[2]=1440817236/23088536315,c[3]=3309984816/23088536315,c[4]=13025487078/23088536315,c[5]=554028/23088536315,c[6]=-(2605540638/4617707263),c[7]=-(3302228424/23088536315),c[8]=-(1469626692/23088536315),c[9]=-(597939012/23088536315),c[10]=-(401116272/23088536315)};
D2={d[0]=1/26,d[1]=1/17,d[2]=1/10,d[3]=1/5,d[4]=1/2,d[5]=1,d[6]=1/2,d[7]=1/5,d[8]=1/10,d[9]=1/17,d[10]=1/26};

gs=Table[g[A2[[t]],B2[[t]],C2[[t]],D2[[t]],t,x],{t,1,11,1}](*the splines between two points. there are 11 cubic splines*)



cubicinterpol=Show[Plot[gs[[1]],{x,-5,-4},PlotRange->All, PlotStyle->Red,PlotLegends->{" 1st spline"}],
Plot[gs[[2]],{x,-4,-3},PlotRange->All,PlotStyle->{Yellow,Opacity[2]},PlotLegends->{" 2nd spline"}],
Plot[gs[[3]],{x,-3,-2},PlotRange->All,PlotStyle->Green,PlotLegends->{" 3rd spline"}],
Plot[gs[[4]],{x,-2,-1},PlotRange->All,PlotStyle->Blue,PlotLegends->{" 4th spline"}],
Plot[gs[[5]],{x,-1,0},PlotRange->All,PlotStyle->Magenta,PlotLegends->{" 5th spline"}],
Plot[gs[[6]],{x,0,1},PlotRange->All,PlotStyle->Black,PlotLegends->{" 6th spline"}],
Plot[gs[[7]],{x,1,2},PlotRange->All,PlotStyle->Pink,PlotLegends->{" 7th spline"}],
Plot[gs[[8]],{x,2,3},PlotRange->All,PlotStyle->Purple,PlotLegends->{" 8th spline"}],
Plot[gs[[9]],{x,3,4},PlotRange->All,PlotStyle->Orange,PlotLegends->{" 9th spline"}],
Plot[gs[[10]],{x,4,5},PlotRange->All,PlotStyle->Gray,PlotLegends->{" 10th spline"}],
Plot[gs[[11]],{x,5,6},PlotRange->All,PlotStyle->{Green,Opacity[0.5]},PlotLegends->{" 11th spline"}]]
datapoints=Table[{s,Re[f[s]]},{s,-5,6,1}];
scatterplot=ListPlot[datapoints,PlotStyle->Black,PlotLegends->{"datapoints"}]
com2=Show[scatterplot,cubicinterpol]


(*(5)We can use the Built in Interpolation function in Mathematica
 to interpolate the data points using cubic spline interpolation:*)
datapoints=Table[{s,Re[f[s]]},{s,-5,6,1}]
spline = Interpolation[datapoints, Method -> "Spline"]
splinefit=Plot[spline[s], {s, -5, 5},PlotStyle->Gray]
com3=Show[splinefit,scatterplot]
com4=Show[splinefit,cubicinterpol,scatterplot]


(* ::Subsection:: *)
(*(6).Conclusion about the results*)
(*1. If we want to compair between these three method of comparison then the cubic splines method is best for fitting the curves .It is better than linear and Lagrange interpolation method.*)
(*2. Result using built in spline function and defining a cubic spline are same as seen in early figure and it matches perfectly with the function itself.*)
(*3.Lagrange interpolation is  a good approximation method but in end point it fluctuates more.*)
