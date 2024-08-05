#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
(*Assignment 2:*)
(*Name: Santanu Das*)
(*Roll No.: PH22C040*)
(*Q1.*)


(* the data from the CSV file *)

data = Import["D:\\books\\computational physics\\Book1.csv"];
ListPlot[data, PlotStyle -> Red, PlotMarkers -> Automatic, PlotLegends -> {"Data"},AxesLabel -> {"Temperature (K)", "Resistivity (\[CapitalOmega]*m)"}]
Length[data]



(* Step 2: Define the Bloch-Gr\[UDoubleDot]neisen function *)
\[Rho][T_] := r+\[Alpha]*(T/\[Theta])^5*NIntegrate[(t^5 / ((1-Exp[-t])*(Exp[t]-1))),{t,0,\[Theta]/T}]//FullSimplify




(* Step 3: Use the FindFit function to fit the data to the Bloch-Gr\[UDoubleDot]neisen function *)
NonlinearFit=NonlinearModelFit[data,\[Rho][T],{{r,0.5},{\[Alpha],1},{\[Theta],201}},T]



(*parameterTable*)
NonlinearFit["ParameterTable"]


(* Step 4: Plot the data along with the fitted curve *)

Show[ListPlot[data, PlotRange -> All], Plot[NonlinearFit[T], {T, 0, 210}, PlotRange -> All],AxesLabel -> {"Temperature (K)", "Resistivity (\[CapitalOmega]*m)"}]
 Plot[NonlinearFit[T], {T, 0, 1000},PlotStyle->Gray,AxesLabel -> {"Temperature (K)", "Resistivity (\[CapitalOmega]*m)"}](*at high temperature due to electron-phonon interaction resistivity varies linear with temperature*)
 Plot[NonlinearFit[T], {T, 0, 100},PlotStyle->Red,AxesLabel -> {"Temperature (K)", "Resistivity (\[CapitalOmega]*m)"}](*at low temperature due to electron-eledtron interaction resistivity varies with temperature as T^5*)


fittedata=Table[{2*i,NonlinearFit[2*i]},{i,1,Length[data]}]


(*Export the fitted data and original data*)

Export["D:\\books\\computational physics\\fitted data and original data.csv",fittedata]



(* ::Subsection:: *)
(*Q2.*)


(* Consider an infinite conducting plate situatated at z=0 and two charges q1=-2 and q2 =3 are positioned at z=5 and z=10, find the equipotential surfaces and electric field in the region of interest*)
Clear["Global`*"]
e=8.85*10^-12(*the value of permitivity in free space*)
k=1/(4*\[Pi]*e);
q1=-2;
q2=3;
(*the expression of the potential for the charge q1 and q2 and the contribution for the image charge is written below*)
V[x_,y_,z_]:= (q1/(x^2+y^2+(z-5)^2))+(q2/(x^2+y^2+(z-10)^2))+(-q1/(x^2+y^2+(z+5)^2))+(-q2/((x^2+y^2+(z+10)^2)));
(*the electric field is defined below*)
field[x_,y_,z_]=-Grad[V[x,y,z],{x,y,z},"Cartesian"]
field1=-Grad[V[0,y,z],{y,z},"Cartesian"]
field2=-Grad[V[x,0,z],{x,z},"Cartesian"]
ContourPlot[V[x,0,z],{x,-10,10},{z,0,20}](*equipotential curves at x-z plane*)
ContourPlot[V[0,y,z],{y,-10,10},{z,0,20}](*equipotential curves at x-z plane*)
ContourPlot[V[x,y,0],{x,-10,10},{y,-10,10}](*zero potential at z=0 plane*)
VectorPlot3D[field[0,y,z],{x,-10,10},{y,-10,10},{z,0,20}](*In 3D *)
StreamPlot[field1,{y,-10,10},{z,0,20}]
StreamPlot[field2,{x,-10,10},{z,0,20}]
VectorPlot[field1,{y,-10,10},{z,0,20}](* in y-z plane*)
VectorPlot[field2,{x,-10,10},{z,0,20}](* in x-z plane*)










