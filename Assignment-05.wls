#!/usr/bin/env wolframscript
(* ::Package:: *)

(* ::Section:: *)
(*Assignment-05*)
(*Name: Santanu Das*)
(*Roll No.- PH22C040*)
(**)


(* ::Subsection:: *)
(*Q1.*)
(* Here we define an linear congruential generator (LCG)  by the equation*)
(*  Xn \[Congruent] (a \:22c5 Xn-1 + c) mod m, where the values of m (the modulus, a positive integer),*)
(*   a (the multiplier, a positive integer less than m) ,c (the increment, either 0 or a positive integer less than m)*)


Clear["Global`*"]
Subscript[x,0]=x;(*Subscript[x,0] is called seed, it is used as starting point) *)
m=m;
n=5;(* n (the numbers of pseudo random numbers)*) 
pseudorannum={};
For[i=0,i<n,i++,Subscript[x,i+1]=Mod[a*Subscript[x,i]+c,m];
AppendTo[pseudorannum,Subscript[x,i+1]/m]]
pseudorannum



Clear["Global`*"](*We Generate pseudo-random numbers 50000 *)
Subscript[x,0]=0;
m=2^30-1;
a=1254;
c=1234;
n=50000;
pseudorannum={};
For[i=0,i<n,i++,Subscript[x,i+1]=Mod[a*Subscript[x,i]+c,m]//N;
AppendTo[pseudorannum,Subscript[x,i+1]/m]]
pseudorannum;

(* Histogram of the pseudo random numbers *)
Histogram[pseudorannum,Automatic, "PDF"]
ListPlot[pseudorannum]
(*various moments for these pseudo numbers*)
Mean[pseudorannum]//N
Variance[pseudorannum]//N
StandardDeviation[pseudorannum]//N
Skewness[pseudorannum]//N
Kurtosis[pseudorannum]//N


Clear["Global`*"](*We Generate pseudo-random numbers 50000 *)
Subscript[x,0]=0;
m=2^30;
a=1254;
c=1234;
n=50000;
pseudorannum={};
For[i=0,i<n,i++,Subscript[x,i+1]=Mod[a*Subscript[x,i]+c,m]//N;
AppendTo[pseudorannum,Subscript[x,i+1]/m]]
pseudorannum;

(* Histogram of the pseudo random numbers *)
Histogram[pseudorannum,Automatic, "PDF"]
ListPlot[pseudorannum]
(*various moments for these pseudo numbers*)
Mean[pseudorannum]//N
Variance[pseudorannum]//N
StandardDeviation[pseudorannum]//N
Skewness[pseudorannum]//N
Kurtosis[pseudorannum]//N


(* ::Section:: *)
(*Comment : *)


(* ::Subsection:: *)
(*Linear congruential operator produces pseudo random numbers *)
(*but the produced number may or may not be random ,it all depends on how wisely we choose the value of a,c and m in the above expression*)
(**)
(**)


(* ::Section:: *)
(*Q2.Central limit theorem:*)


(* ::Subsection:: *)
(*Let's we have a random variable X (Population) which has a mean \[Mu] and a standard deviation \[Sigma] and if it  follows a certain kind of distribution (it could be anything uniform, Gaussian,Poisson,exponential or any random distribution etc.).Now if we make a sampling of sample n (preferable if n>= 30)*)
(*and do the average of each sampling group and try to find out it's distribution, it would be a Gaussian with same mean \[Mu] and standard deviation \[Sigma]/Sqrt[n]*)
(*if we increase the sampling group size n then the fit would be eventually better. *)


Clear["Global`*"]

popudata1=RandomReal[{1,10},30000];(*here, i am taking a large no. of randomreal and 
                                     found out what kind of distribution it follows using histogram plot*)
plot1=Histogram[popudata1]
Mean[popudata1]       (* here, the mean and standard deviation*)
StandardDeviation[popudata1]

popudata2=RandomVariate[PoissonDistribution[25],30000];(*Here, i am taking poisson distribution and taking some large random no. 
                                                      i also plot the histogram to show the behaviour of the numbers which i take.*)
plot1=Histogram[popudata2]
Mean[PoissonDistribution[25]](* here, the mean and standard deviation*)
StandardDeviation[PoissonDistribution[25]]//N


Clear["Global`*"](*here i have done group sampling the two random variables for n=10,n=30,n=50 and averaging the sample *)
(*1*)
sampling01={};
Group={};
For[j=1,j<101,j++,For[n=1,n<11,n++,s=RandomReal[{1,10}];AppendTo[Group,s]];AppendTo[sampling01,Group]] 
Samplemeans01=Mean/@sampling01;

sampling02={};
Group1={};
For[j=1,j<101,j++,For[n=1,n<31,n++,\[NonBreakingSpace]s=RandomReal[{1,10}];AppendTo[Group1,s]];AppendTo[sampling02,Group1]] 
Samplemeans02=Mean/@sampling02;

sampling03={};
Group2={};
For[j=1,j<101,j++,For[n=1,n<51,n++,s=RandomReal[{1,10}];AppendTo[Group2,s]];AppendTo[sampling03,Group2]]
Group2;
Samplemeans03=Mean/@sampling03;


(*2*)
sampling11={};
For[i=1,i<101,i++,k=RandomVariate[PoissonDistribution[25],10];AppendTo[sampling11,k]]

Samplemeans11=Mean/@sampling11//N

sampling12={};
For[i=1,i<101,i++,k=RandomVariate[PoissonDistribution[25],30];AppendTo[sampling12,k]]

Samplemeans12=Mean/@sampling12//N

sampling13={};
For[i=1,i<101,i++,k=RandomVariate[PoissonDistribution[25],50];AppendTo[sampling13,k]]

Samplemeans13=Mean/@sampling13//N


(*1*)Histogram[Samplemeans01,20, "PDF"](*in this histogram plot their is a left skewness 
                                          that means the won't fit perfectly as the sample size is 10*)
Mean[Samplemeans01]
StandardDeviation[Samplemeans01]

Histogram[Samplemeans02,20, "PDF"]

Mean[Samplemeans02]
StandardDeviation[Samplemeans02]

Histogram[Samplemeans03,20, "PDF"]

Mean[Samplemeans03]
StandardDeviation[Samplemeans03]

(*2*)Histogram[Samplemeans11,20, "PDF"]

Mean[Samplemeans11]
StandardDeviation[Samplemeans11]

Histogram[Samplemeans12,20, "PDF"]

Mean[Samplemeans12]
StandardDeviation[Samplemeans12]

Histogram[Samplemeans13,20, "PDF"]

Mean[Samplemeans13]
StandardDeviation[Samplemeans13]


(* ::Subsection:: *)
(*Q3.Buffon's Needle Problem:*)


Clear["Global`*"]
Clear["Global`*"]
g1=Graphics[{{Red,Line[{{0, 0}, {0, 8}}],Line[{{8, 0}, {8, 8}}],Line[{{0, 0}, {8, 0}}],
Line[{{0, 8}, {8, 8}}]}, {Green,Line[{{2, 0}, {2, 8}}],Line[{{4, 0}, {4, 8}}],Line[{{6, 0}, {6, 8}}],{Blue,needles}}
}];(*red line definig for the boundary *)
d=2;(*the difference between two lines*)
L=1;(*length of the needle*)
n=200;(*number of needles*)
xvals=RandomReal[{1,7},n];
yvals=RandomReal[{1,7},n];
\[Theta]vals=RandomReal[{0,180\[Degree]},n];
coordinates=Table[{xvals[[i]],yvals[[i]],\[Theta]vals[[i]]},{i,1,n}];
needle[ x_, y_, \[Theta]_] := Arrow[{{x-(L/2)*Cos[\[Theta]], y-(L/2)*Sin[\[Theta]]},
{x + (L/2)*Cos[\[Theta]], y +(L/2)*Sin[\[Theta]]}}](*a single needle is defined*)
nee1=Dynamic[Graphics[g1,needles]]
 needles={};
 For[h=1,h<n+1,h++,s=needle[xvals[[h]],yvals[[h]],\[Theta]vals[[h]]];AppendTo[needles,s];Pause[0.1]]
 needles;
nee=Graphics[{Blue,needles}];(*here we get needles*)

(*intersection condition of needles and the lines*)
intersect={};
sp={};
For[i=1,i<n+1,i++,r=If[Abs[xvals[[i]]-2]<=(L/2)*Sin[\[Theta]vals[[i]]]||
Abs[xvals[[i]]-4]<=(L/2)*Sin[\[Theta]vals[[i]]]||Abs[xvals[[i]]-6]<=(L/2)*Sin[\[Theta]vals[[i]]],1,0];
AppendTo[intersect,r]];
intersect;
n1=Count[intersect,1];(*intersecting-point=n1*)
valpi=(n/n1)//N(*value of \[Pi]*)
actuval=N[Pi]
percenerror=(Abs[actuval-valpi]/actuval)*100//N




(* ::Subsection:: *)
(*Q4.*)
(* Monte Carlo integration method.*)


f[x_]:=x^2; (* here we are defining  the function *)
n=100000; (* number of random points within the define x and y intervals*)
n1=0; (* n1=points that fall under the function's area *)
a=0;
b=1;
c=0;
d=1;
h=d-c;
For[i=1,i<=n,i++,x=RandomReal[{a,b}]; y=RandomReal[{c,d}];If[y <= f[x],n1++]]
integral1 =N[h*(b-a)] * N[n1/n](* estimate the integral  using monte carlo method *)


(* ::Subsection:: *)
(*i)Integration using Trapezoidal rule*)


f[x_]:= x^2; (* here we are defining  the function *)
a1=0; (* lower limit of integration *)
b1=1; (* upper limit of integration *)
n2=100; 
h2=(b1-a1)/n2;(* number of intervals *)
sum=(f[a1]+f[b1])/2;
For[i=1,i<n2,i++,sum=sum + f[a1 + i*h2]]
integral2 = N[h2*sum](* estimate the integral using the trapezoidal rule *)


(* ::Subsection:: *)
(*ii)Integration using Simpson 1/3 rd rule*)


f[x_]:=x^2; (* here we are defining  the function *)
a2=0; (* lower limit of integration *)
b2=1; (* upper limit of integration *)
n3=100; 
h3=(b2 - a2)/n3;(* number of intervals *)
sum1=(h3*(f[a2]+f[b2]))/3;
sum02=0;
sum03=0;
For[i=1,i<=(n3-2)/2,i++,sum02=sum02+f[a2+2*i*h3]]
sum2=(h3*2*sum02)/3;
For[i=1,i<=(n3-2)/2,i++,sum03=sum03+f[a2+(2*i+1)*h3]]
sum3=(h3*4*sum03)/3;
integral3=sum1+sum2+sum3//N


(* using build in function *)
actualval=NIntegrate[x^2, {x, 0, 1}]


(*percentage error*)
percentageerror1=(Abs[actualval-integral1]/actualval)*100
percentageerror2=(Abs[actualval-integral2]/actualval)*100
percentageerror3=(Abs[actualval-integral3]/actualval)*100


(* ::Subsubsection:: *)
(*Comment: *)
(*1). In Monte Carlo method is the error is fluctuating because the points are random*)
(*       but if we increase the no. points error is reduced.*)
(*2).  In case of trapezoidal and Simpson 1/3 rd method ,If we increase the number of interval*)
(*       the error is reduced gradually.*)
(*      The error is less Simpson method compare to trapezoidal method.*)
(* 3).But both in trapezoidal and Simpson the error is less than Monte Carlo method.     *)


(* ::Section:: *)
(**)
(*Q5.  Metropolis-Hastings Algorithm*)


(* ::Subsection:: *)
(*The Metropolis-Hastings algorithm is a method used to generate random samples from a probability distribution(called target distribution) that is difficult to directly sample from. It is used in statistical modeling.*)
(* Metropolis-Hastings algorithm is to create a Chain where each step depends on the previous step and the proposed new step.*)
(* At each step, the algorithm creates a new value from the target distribution, based on a "proposal distribution"(that we have give to find out the new values).*)
(*  Then, it calculates the ratio of the target distribution at the new proposed value to the target distribution at the current value. If this ratio is greater than or equal to 1, the new proposed value is accepted. Otherwise, the this value is accepted with a probability equal to the ratio we calculate before. *)
(*In simple terms, the Metropolis-Hastings algorithm repeatedly generates new samples from a proposal distribution, and then decides whether to accept or reject them based on a comparison of the proposed value to the current value in the chain.And eventually this chain moves towards the desired target distribution.*)
