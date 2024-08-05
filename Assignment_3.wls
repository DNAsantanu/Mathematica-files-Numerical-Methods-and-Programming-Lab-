#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Chapter:: *)
(*Assignment -3*)
(*Name:-Santanu Das*)
(*Roll No.-PH22C040*)
(**)


(* ::Subsection:: *)
(*Q1.*)


(*Lagrangian Interpolation*)
Clear["Global`*"]
X={Subscript[x, 0],Subscript[x, 1],Subscript[x, 2],Subscript[x, 3]}
F={Subscript[f, 0],Subscript[f, 1],Subscript[f, 2],Subscript[f, 3]}
g[x_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(4\)]\(F[\([i]\)]*\(
\*UnderoverscriptBox[\(\[Product]\), \(j = 1\), \(4\)]\((If[i == j, 1, \((x - X[\([j]\)])\)/\((X[\([i]\)] - X[\([j]\)])\)])\)\)\)\)
g[x]



(* ::Subsection:: *)
(*Q2.*)


Clear["Global`*"]
(*Creating six pair of data points from a second order polynomial*)
f[x_]:=0.034+2.099x-45x^2+5x^3//N
dataset=Table[{i,f[i]},{i,1,12,2}](*6 pair of data point*)
datax=Table[i,{i,1,12,2}](*x values*)
datay=Table[f[i],{i,1,12,2}](*y values*)
p1=Plot[f[x],{x,-1,10},AxesLabel->{"x","f[x]"}]


(*Defining the lagrangian Basis function*)
lag[x_,i_]:=\!\(
\*UnderoverscriptBox[\(\[Product]\), \(j = 1\), \(6\)]\((If[i == j, 1, \((x - datax[\([j]\)])\)/\((datax[\([i]\)] - datax[\([j]\)])\)])\)\)
lag[x,1]//Expand
lag[x,2]//Expand
lag[x,3]//Expand
lag[x,4]//Expand
lag[x,5]//Expand
lag[x,6]//Expand
Plot[{lag[x,1],lag[x,2],lag[x,3],lag[x,4],lag[x,5],lag[x,6]},{x,0,10},PlotLegends->"Expressions"]


(*The fitted function via lagrange interpolation*)
g[x_]:=\!\(
\*UnderoverscriptBox[\(\[Sum]\), \(i = 1\), \(6\)]\(datay[\([i]\)]*\(
\*UnderoverscriptBox[\(\[Product]\), \(j = 1\), \(6\)]\((If[i == j, 1, \((x - datax[\([j]\)])\)/\((datax[\([i]\)] - datax[\([j]\)])\)])\)\)\)\)
g[x]//Expand
p2=Plot[g[x],{x,-1,10},PlotLegends->"Expressions",AxesLabel->{"x","g[x]"},PlotStyle->Red]


(*Now i choose 3 pair of data points in between the given range which is given below*)
expdata=Table[{k,f[k]},{k,4.3,9,2}]
(*Now pair of points found from the fitted function*)
fitteddata=Table[{k,g[k]},{k,4.3,9,2}](*interpolated data*)



(*error calculation*)
e[x_]:=f[x]-g[x]
error=SetPrecision[Table[e[i],{i,4.3,9,2}],9](*error is very low order*)


(* ::Subsection:: *)
(**)
