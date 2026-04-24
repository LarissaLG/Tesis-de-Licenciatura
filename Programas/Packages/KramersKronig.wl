(* ::Package:: *)

BeginPackage["KramersKronig`"];


kkrebook::usage = "kkrebook[omega, parteImag] computes the real part of the dielectric function from its imaginary part using the Kramers-Kronig relations. 
					The input lists for frequencies (omega) and the imaginary part (parteImag) must be equispaced.";

kkimbook::usage = "kkimbook[omega, parteRe] computes the imaginary part of the dielectric function from its real part using the Kramers-Kronig relations. 
					The input lists for frequencies (omega) and the real part (parteRe) must be equispaced.";

selfconsbook::usage = "selfconsbook[omega, realpart, imagpart, N, mu] iteratively computes a self-consistent dielectric function. 
						It requires an equispaced frequency list (omega), initial experimental  arrays for the real and imaginary parts, 
						the number of iterations (N), and a weighting factor (mu) between 0 and 1 to balance the experimental data with the calculated KK values. 
						Returns a list containing the refined {real part, imaginary part}.";

sskkrebook::usage = "sskkrebook[omega, parteImag, omega1, parteReal1] applies the Singly Subtractive Kramers-Kronig (SSKK) relation to compute the real part of 
							the dielectric function. It requires equispaced input arrays for frequencies and the imaginary part, 
								plus a single known reference point {omega1, parteReal1}.";
					
sskkimbook::usage = "sskkimbook[omega, parteRe, omega1, parteImag1] applies the Singly Subtractive Kramers-Kronig (SSKK) relation to compute the imaginary part of
					 the dielectric function. It requires equispaced input arrays for frequencies and the real part, 
						plus a single known reference point {omega1, parteImag1}.";


Begin["Private`"];


kkrebook[omega_List,parteImag_List]:=Module[{g,parteRe,a,b,deltaomega,j,k},
									(*Asegura que omega y parteImag sean listas horizontales*)
									If[MatrixQ[{omega}],omega=Flatten[omega]];
									If[MatrixQ[{parteImag}],parteImag=Flatten[parteImag]];
									g=Length[omega];
									parteRe=ConstantArray[0.,g];
									a=ConstantArray[0.,g];
									b=ConstantArray[0.,g];
									deltaomega=omega[[2]]-omega[[1]];
									
									(*Primer punto (excluye omega[[1]])*)

									b[[1]]=Sum[parteImag[[k]]*omega[[k]]/(omega[[k]]^2-omega[[1]]^2),{k,2,g}];
									parteRe[[1]]=(2/Pi)*deltaomega*b[[1]]+1;
									
									(*\[CapitalUAcute]ltimo punto (excluye omega[[g]])*)
									
									a[[g]]=Sum[parteImag[[k]]*omega[[k]]/(omega[[k]]^2-omega[[g]]^2),{k,1,g-1}];
									parteRe[[g]]=(2/Pi)*deltaomega*a[[g]]+1;
									
									(*Puntos intermedios*)
									For[j=2,j<=g-1,j++,a[[j]]=Sum[parteImag[[k]]*omega[[k]]/(omega[[k]]^2-omega[[j]]^2),{k,1,j-1}];
									b[[j]]=Sum[parteImag[[k]]*omega[[k]]/(omega[[k]]^2-omega[[j]]^2),{k,j+1,g}];
									
									parteRe[[j]]=(2/Pi)*deltaomega*(a[[j]]+b[[j]])+1;];
						parteRe]


kkimbook[omega_List,parteRe_List]:=Module[{g,parteImag,a,b,deltaomega,j,k},(*Asegura que sean listas horizontales*)
									If[MatrixQ[{omega}],omega=Flatten[omega]];
									If[MatrixQ[{parteRe}],parteRe=Flatten[parteRe]];
									g=Length[omega];
									parteImag=ConstantArray[0.,g];
									a=ConstantArray[0.,g];
									b=ConstantArray[0.,g];
									deltaomega=omega[[2]]-omega[[1]];
									
									(*Primer punto (excluye omega[[1]])*)
									
									b[[1]]=Sum[(parteRe[[k]]-1)/(omega[[k]]^2-omega[[1]]^2),{k,2,g}];
									parteImag[[1]]=(-2*deltaomega*b[[1]]*omega[[1]])/Pi;
									
									(*\[CapitalUAcute]ltimo punto (excluye omega[[g]])*)
									
									a[[g]]=Sum[(parteRe[[k]]-1)/(omega[[k]]^2-omega[[g]]^2),{k,1,g-1}];
									parteImag[[g]]=(-2*deltaomega*a[[g]]*omega[[g]])/Pi;
									
									(*Puntos intermedios*)
									
									For[j=2,j<=g-1,j++,a[[j]]=Sum[(parteRe[[k]]-1)/(omega[[k]]^2-omega[[j]]^2),{k,1,j-1}];
									b[[j]]=Sum[(parteRe[[k]]-1)/(omega[[k]]^2-omega[[j]]^2),{k,j+1,g}];
									parteImag[[j]]=(-2*deltaomega*(a[[j]]+b[[j]])*omega[[j]])/Pi;];
						parteImag]

selfconsbook[omega_List,nreal_List,kimag_List,N_Integer,mu_?NumericQ]:=Module[{comodo1,comodo2,refin,imfin,j},
									(*Asegura vectores planos*)
									omega=Flatten[omega];
									nreal=Flatten[nreal];
									kimag=Flatten[kimag];
									
									comodo1=nreal;(*Parte real inicial*)
									comodo2=kimag;(*Parte imaginaria inicial*)
									
									For[j=1,j<=N,j++,(*Estima nueva parte real desde la parte imaginaria actual*)
									comodo1=kkrebook[omega,comodo2];
									comodo1=mu*nreal+(1-mu)*comodo1;
									
									(*Estima nueva parte imaginaria desde la parte real actualizada*)
									comodo2=kkimbook[omega,comodo1];
									comodo2=mu*kimag+(1-mu)*comodo2;];
									
									{refin,imfin}={comodo1,comodo2};
						{refin,imfin}]


sskkrebook[omega_List,kimag_List,omega1_,nreal1_]:=Module[{nreal,k},
									(*Encontrar \[IAcute]ndice m\[AAcute]s cercano a omega1*)
									k=First@FirstPosition[omega,_?(Abs[#-omega1]<10^-6&),Missing["NotFound"]];
									If[MissingQ[k],Return[$Failed,Module]];
									
									(*Calcular la parte real con KK normal*)
									nreal=kkrebook[omega,kimag];
									(*Valor en omega1*)
									
									(*Correcci\[OAcute]n substractiva*)
									nreal=nreal+(nreal1-nreal[[k]]);
							nreal]

sskkimbook[omega_List,nreal_List,omega1_,kimag1_]:=Module[{kimag,k},
									(*Encontrar \[IAcute]ndice m\[AAcute]s cercano a omega1*)
									k=First@FirstPosition[omega,_?(Abs[#-omega1]<10^-6&),Missing["NotFound"]];
									
									(*Calcular parte imaginaria inicial con Kramers-Kronig*)
									kimag=kkimbook[omega,nreal];
									
									(*Aplicar correcci\[OAcute]n sustractiva*)
									kimag=kimag+(omega/omega1)*(kimag1-kimag[[k]]);
							kimag]


End[];


EndPackage[];
