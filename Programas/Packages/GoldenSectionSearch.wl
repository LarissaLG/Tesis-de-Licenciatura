(* ::Package:: *)

BeginPackage["GoldenSectionSearch`"];


GoldenSearchLowUppIterMin::usage = "GoldenSearchLowUppIterMin[lower, upper, iterations, tolRes] computes the local minimum of a globally defined function f[x] 
							within the interval [lower, upper] using the Golden Section Search algorithm. It iteratively narrows the search interval until the 
							width is less than tolRes or the maximum number of iterations is reached. 
							It prints the iteration history and returns the x-value of the estimated minimum.";
						
GoldenSearchMin::usage = "GoldenSearchMin[f, lower, upper,tolRes] computes the local minimum of the function f within the interval [lower, upper] 
						using the Golden Section Search algorithm. It relies on a resolution 'tolRes' for its stopping criterion and prints the estimated x-value.";

GoldenSearchMax::usage = "GoldenSearchMax[f, lower, upper,tolRes] computes the local maximum of the function f within the interval [lower, upper] 
						using the Golden Section Search algorithm. It relies on aresolution 'tolRes' for its stopping criterion and returns the estimated x-value.";


Begin["Private`"];


GoldenSearchLowUppIterMin[f_,lower_, upper_, iterations_, tolRes_] := 
 Module[{a = N[lower], b = N[upper], c, d, k = 1},
  
  Print["Iteration Results are:"];
  
  (* El While revisa TANTO la tolerancia como el n\[UAcute]mero de iteraciones *)
  While[
   Abs[(b - (b - a)/GoldenRatio) - (a + (b - a)/GoldenRatio)] > tolRes && k <= iterations,
   
   c = b - (b - a)/GoldenRatio;
   d = a + (b - a)/GoldenRatio;
   
   If[f[c] < f[d],
    (* Si es menor, ajusta el l\[IAcute]mite superior *)
    b = d,
    (* Si es mayor o igual, ajusta el l\[IAcute]mite inferior *)
    a = c
   ];
   
   (* Ahora s\[IAcute] imprime en TODAS las iteraciones *)
   Print[{"Lowerbound at iteration...", k, PaddedForm[a, {7, 6}]}];
   Print[{"Upperbound at iteration...", k, PaddedForm[b, {7, 6}]}];
   Print["Midpoint: ", (b + a)/2];
   
   k++;
  ];
  
  (* Devuelve el valor \[OAcute]ptimo (punto medio del intervalo final) *)
  Return[(b + a)/2];
 ]
 
 
 
 
 GoldenSearchMin[f_, lower_, upper_, tolRes_] := Module[{a = N[lower], b = N[upper], c, d},
  While[
   Abs[(b - (b - a)/GoldenRatio) - (a + (b - a)/GoldenRatio)] > tolRes,
   
   c = b - (b - a)/GoldenRatio;
   d = a + (b - a)/GoldenRatio;
   
   If[f[c] < f[d],
    b = d,
    a = c
   ];
  ];
  Return[(b + a)/2];
]

GoldenSearchMax[f_, lower_, upper_, tolRes_] := Module[{a = N[lower], b = N[upper], c, d},
  While[
   Abs[(b - (b - a)/GoldenRatio) - (a + (b - a)/GoldenRatio)] > tolRes,
   
   c = b - (b - a)/GoldenRatio;
   d = a + (b - a)/GoldenRatio;
   
   (* La \[UAcute]nica diferencia matem\[AAcute]tica es el signo > en lugar de < *)
   If[f[c] > f[d],
    b = d,
    a = c
   ];
  ];
  Return[(b + a)/2];
]
 


End[];


EndPackage[];
