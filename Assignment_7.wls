#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Subsection:: *)
(*Q1. Here i choose  heat equation and solve it in the particular region.*)
(*i choose Neumann Boundary condition to solve it. *)


coordinates={{179.2`,300.0196972126255`},{174.39999999999998`,318.6196972126255`},{191.2`,304.2196972126254`}
,{185.8`,288.6196972126255`},{226.`,314.41969721262546`},{231.39999999999998`,364.2196972126254`},{213.39999999999998`,441.61969721262545`}
,{244.`,366.6196972126255`},{241.59999999999997`,304.81969721262544`},{191.2`,279.6196972126255`},{275.2`,279.6196972126255`},
{268.59999999999997`,360.0196972126255`},{289.`,315.6196972126255`},{286.59999999999997`,266.4196972126255`},{190.`,266.4196972126255`},
{319.`,236.41969721262552`},{294.4`,123.01969721262552`},{217.59999999999997`,12.61969721262551`},{274.59999999999997`,121.81969721262553`},
{289.59999999999997`,225.01969721262552`},{187.`,256.81969721262556`},{268.59999999999997`,210.01969721262552`},{259.`,138.61969721262554`},
{223.`,56.41969721262558`},{245.2`,127.81969721262553`},{252.39999999999998`,207.61969721262554`},{181.60000000000002`,251.41969721262552`},
{216.39999999999998`,222.61969721262554`},{213.39999999999998`,173.41969721262555`},{166.60000000000002`,126.61969721262554`},
{116.8`,172.81969721262553`},{115.60000000000001`,222.01969721262552`},{145.`,248.41969721262552`},{82.00000000000001`,208.2196972126255`},
{83.8`,148.81969721262553`},{107.8`,55.819697212625556`},{74.2`,130.81969721262553`},{62.2`,213.01969721262552`},{143.2`,258.0196972126255`},
{43.000000000000014`,225.01969721262552`},{60.400000000000006`,114.01969721262552`},{106.00000000000001`,24.019697212625488`},
{36.400000000000006`,129.61969721262554`},{16.000000000000014`,237.0196972126255`},{137.8`,265.81969721262556`},
{42.400000000000006`,266.4196972126256`},{42.400000000000006`,313.21969721262553`},{61.60000000000001`,358.81969721262556`},
{58.000000000000014`,282.0196972126255`},{140.2`,277.21969721262553`},{89.2`,304.81969721262556`},{87.40000000000002`,364.81969721262556`},
{113.8`,432.0196972126255`},{101.8`,358.2196972126255`},{107.8`,314.41969721262546`},{146.2`,288.0196972126255`},
{138.39999999999998`,301.81969721262544`},{155.8`,316.2196972126255`},{151.60000000000002`,297.0196972126255`},
{165.39999999999998`,295.2196972126254`}};
spidy=Polygon[coordinates];(*use polygon function to connect the coordinates*)
Graphics[spidy]


T=10;(*the time need to express the solution*)
eqn = D[u[x, y, t], t] == Laplacian[u[x, y, t], {x, y}]+NeumannValue[12,True];
(* Neumann Boundary condition
 i have taken and the value of this is 12 (it specifies a flux across the edges if the boundary*)
ic = u[x, y, 0] == 5;(*initial condition*)
sol = NDSolve[{eqn,  ic}, u, {x,y}\[Element]spidy, {t, 0, T}];

Plot3D[Evaluate[u[x, y, 2] /. sol],{x, y}\[Element]spidy ,
  PlotRange -> All, PlotLegends->"Gold spidy",AxesLabel -> {"x", "y","u"}]
ContourPlot[Evaluate[u[x, y, 1] /. sol], {x, y}\[Element]spidy,
            PlotRange -> All, ColorFunction->"Hue"]                   
           


(* ::Subsection:: *)
(*Q2. SVD computation:*)


Clear["Global`*"]
A={{0,1,0},{0,1,0},{1,1,0},{1,0,1},{0,1,1}};
a=MatrixForm[A];
TransA=Transpose[A];
a1=MatrixForm[TransA];

(*computation of V*)
B=TransA . A;
{eigenvals,eigenvecs} = Eigensystem[B];
eigenvecs(*eigenvectors*);
V=Orthogonalize[eigenvecs];(*normalize the eigen vector*)
V[[1]]//MatrixForm;
TransV=Transpose[V]//MatrixForm

(*computation of sigma matrix*)
eigenvals(*eigenvalues*);
singularvals={};
For[k=1,k<=5,k++,For[i=1,i<=3,i++,If[i==k,a=Sqrt[eigenvals[[i]]],a=0];AppendTo[singularvals,a]]]
sigmamatrix=Partition[singularvals,3]//MatrixForm

(*computation of U*)
B=A . TransA;
b=A . TransA//MatrixForm;
{eigenvals,eigenvecs}=Eigensystem[B];
eigenvals;
eigenvecs;(*eigenvectors*)
eigenvecs1=Delete[eigenvecs,4];
orthonormality={};
orthonormality1={};
For[i=1,i<=4,i++,For[j=1,j<=4,j++,If[i!=j,
p=eigenvecs[[i]] . Transpose[eigenvecs[[j]]];AppendTo[orthonormality,p]]]](*the set {x1,x2,x3,x4} is a orthonormal set*)
For[i=1,i<=4,i++,For[j=1,j<=4,j++,If[i!=j
,p=eigenvecs1[[i]] . Transpose[eigenvecs1[[j]]];AppendTo[orthonormality1,p]]]](*the set {x1,x2,x3,x5} is a orthonormal set*)
orthonormality
orthonormality1
eigenvecs[[4]] . Transpose[eigenvecs[[5]]](*the last two eigen vectors are orthogonal to each other*)
U1=Orthogonalize[eigenvecs];(*Orthogonalize the eigen vector*)
U1[[4]] . Transpose[U1[[5]]](*Now the last two eigen vectors are orthogonal to each other*)
U=Transpose[U1]//MatrixForm(*U matrix*)
u=Transpose[U1];
sigmamatrix=Partition[singularvals,3]//MatrixForm(*S matrix*)
s=Partition[singularvals,3];
TransV=Transpose[V]//MatrixForm(*V matrix*)
transv=V;(*V^T matrix*)
(*Comprison*)
u . s . transv//MatrixForm
A//MatrixForm


Clear["Global`*"]
(*using bulid in function*)
A={{0,1,0},{0,1,0},{1,1,0},{1,0,1},{0,1,1}};
a=MatrixForm[A]
{U,S,TransV}=SingularValueDecomposition[A];
U//MatrixForm
S//MatrixForm
TransV//MatrixForm
U . S . Transpose[TransV]//MatrixForm


(* ::Subsubsection:: *)
(*Application of SVD: (Image Processing)*)


Clear["Global`*"]
img = Import["D:\\books\\computational physics\\spidy.jpg"]


grayImg = ColorConvert[img, "Grayscale"]
mat = ImageData[grayImg];

(*decompose the imagedata*)
{u, s, v} = SingularValueDecomposition[mat];

k = 20;
dias=Diagonal[s];
news = Take[dias,k];
newu = Take[u,All,k];
newtransv = Take[v,All, k];
compressedMat=newu . DiagonalMatrix[news] . Transpose[newtransv];
compressedImg=Image[compressedMat,ImageSize->Medium]


(* ::Subsection:: *)
(*1.*)
(*Suppose we have an m * n matrix A, SVD decomposes this matrix as A = U\[CapitalSigma]V^T, where U is an m * m orthogonal matrix, \[CapitalSigma] is an m * n diagonal matrix with non-negative entries known as the singular values, and V is an n * n orthogonal matrix.*)
(**)
(*SVD is considered one of the powerful tools in data analysis for several reasons:*)
(**)
(*Dimensionality Reduction: SVD can be used to reduce the dimensionality of a dataset by retaining only the top k no of singular values and their corresponding singular vectors. This reduces the amount of noise in the data and can improve the efficiency of algorithms that work with high-dimensional data.*)


(* ::Subsection:: *)
(*2.*)
(*PCA (Principal Component Analysis) and SVD (Singular Value Decomposition) are both linear algebra techniques used for dimensionality reduction, data compression, noise reduction etc.*)
(**)
(**)
(*1)PCA is a statistical technique that focuses on finding a new set of uncorrelated variables called principal components that capture the maximum variance in the original data. PCA is based on the covariance matrix of the data and uses eigen decomposition to obtain the principal components.*)
(*SVD, on the other hand, is a matrix factorization technique that decomposes a matrix into three matrices: U, \[CapitalSigma], and V. Here, U and V are orthogonal matrices, and \[CapitalSigma] is a diagonal matrix containing the singular values of the original matrix. SVD can be used for data compression, low-rank approximation, and solving linear systems of equations.*)
(**)
(**)
(*2)Both PCA and SVD can be used for dimensionality reduction and feature extraction, PCA is specifically designed for covariance matrices, while SVD decomposes a matrix into three matrices and is more general.*)
