(* ::Package:: *)

BeginPackage["LorentzFitSSKK`"];


ajusteSinPico::usage =
"ajusteSinPico[epsIm, intervalos, tolRes, xmin, xmax] fits the imaginary part of the dielectric function with five Lorentz oscillators, excluding the small Soret peak.";

ajustePico::usage =
"ajustePico[epsIm] fits the small Soret peak using a single Lorentz oscillator.";

ajusteTotal::usage =
"ajusteTotal[...] performs the complete Lorentzian fit using the results of ajusteSinPico and ajustePico as initial parameters.";


Begin["Private`"];


(*Directorio donde est\[AAcute]n las paqueter\[IAcute]as*)
SetDirectory["C:\\Users\\danal\\Documents\\Github\\Tesis-de-Licenciatura\\Programas\\Packages"];

Needs["KramersKronig`"];
Needs["GoldenSectionSearchMaxMin`"];


(*Funciones auxiliares*)

(*Funci\[OAcute]n para encontrar los m\[AAcute]ximos alrededor de intervalos dados*)
omega0i[epsFunc_, intervalos_, tolRes_] := GoldenSectionSearchMax[epsFunc, #1, #2, tolRes] & @@@ intervalos;

(*Parte imaginaria de lorentziana*)
lorentzIm[omega_,omega0_,gamma_]:=(gamma omega)/((omega0^2-omega^2)^2+(gamma omega)^2)
(*Suma lorentzianas*)
lorentzSumIm[omega_, A_List, omega0_List, gamma_List] :=  Sum[A[[i]] * lorentzIm[omega, omega0[[i]], gamma[[i]]], {i, Length[A]}]


(*Ajuste lorentzianas*)

(*Ajuste ignorando el pico peque\[NTilde]o*)

(*Intervalos de b\[UAcute]squeda de los m\[AAcute]ximos*)

ajusteSinPico[epsIm_,intervalos_,tolRes_,xmin_,xmax_]:=Module[{epsData0,epsData1,omega0i0,epsData2,fitSinPico,fitSinPicoValues},														
															epsData0=Table[{x,epsIm[x]},{x,xmin,2.65,0.01}];
															epsData1=Table[{x,epsIm[x]},{x,2.76,xmax,0.01}];
															omega0i0=omega0i[epsIm,intervalos,tolRes];
															epsData2=Join[epsData0,epsData1];
															fitSinPico = NonlinearModelFit[epsData2,
												                           lorentzSumIm[omega, {A1, A2, A3, A4, A5}, {omega01,omega02,omega03,omega04,omega05}, 
																	                           {gamma1, gamma2, gamma3, gamma4, gamma5} ], 
												                          {{A1, 0.01}, {omega01, omega0i0[[1]]}, {gamma1, 0.05}, 
												                           {A2, 0.01}, {omega02, omega0i0[[2]]}, {gamma2, 0.1}, 
												                           {A3, 0.01}, {omega03, omega0i0[[3]]}, {gamma3, 0.1}, 
												                           {A4, 0.01}, {omega04, omega0i0[[4]]}, {gamma4, 0.1}, 
												                            {A5, 0.01}, {omega05, omega0i0[[5]]}, {gamma5, 0.1}}, 
												                            omega];
												             fitSinPicoValues = Values[fitSinPico["BestFitParameters"]];

															{epsData2,fitSinPico, fitSinPicoValues}]
												                           
ajustePico[epsIm_]:=Module[{epsData3,fitPico,fitPicoValues},
							epsData3=Table[{x,epsIm[x]},{x,2.72,2.74,0.0007}]; 												
							fitPico = NonlinearModelFit[epsData3,A1*lorentzIm[omega,omega01,gamma1],{{A1,0.0001},{omega01,2.75},{gamma1,0.01}},omega];
							fitPicoValues = Values[fitPico["BestFitParameters"]];{epsData3,fitPico,fitPicoValues}]
									


ajusteTotal[epsIm_,intervalos_,tolRes_,xmin_,xmax_]:=Module[{epsData,sinPico,pico,fitTotal},
															epsData=Table[{x,epsIm[x]},{x,xmin,xmax,0.01}];
															sinPico=ajusteSinPico[epsIm,intervalos,tolRes,xmin,xmax][[3]];
															
															pico=ajustePico[epsIm][[3]];
																		
															
															fitTotal=NonlinearModelFit[epsData,
															{lorentzSumIm[omega, {A1, A2, A3, A4, A5,A6}, {omega01,omega02,omega03,omega04,omega05,omega06}, 
																										{gamma1, gamma2, gamma3, gamma4, gamma5,gamma6} ], 
											                         A2>0 && 0.002>A6>0 && 0.2>gamma6>0 && 2.77>omega06>2.66}, 
											                         {{A1,sinPico[[1]]},{omega01,sinPico[[2]]},{gamma1,sinPico[[3]]},
											                         {A2,sinPico[[4]]},{omega02,sinPico[[5]]},{gamma2,sinPico[[6]]},
											                         {A3,sinPico[[7]]},{omega03,sinPico[[8]]},{gamma3,sinPico[[9]]},
											                         {A4,sinPico[[10]]},{omega04,sinPico[[11]]},{gamma4,sinPico[[12]]},
											                         {A5,sinPico[[13]]},{omega05,sinPico[[14]]},{gamma5,sinPico[[15]]},
											                         {A6,pico[[1]]},{omega06,pico[[2]]},{gamma6,pico[[3]]}},
											                         omega, Weights->(If[2.8<#<3.1,100,2]&/@(epsData[[All,1]]))];
											                  {epsData,fitTotal, Values[fitTotal["BestFitParameters"]]}]


End[];


EndPackage[];
