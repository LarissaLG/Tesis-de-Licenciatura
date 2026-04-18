(* ::Package:: *)

BeginPackage["MieTheory`"];


mieQsca::usage = "mieQsca[np,nm, a, lambda] recibes as first entry the particle refractive index, then the medium refractive index, the sphere radius and the wavelength. 
					As a result it gives for results: the an coefficient list for first index, the bn coefficient for second index, scattering cross section for third index and
					the scattering efficiency for the fourth index. ";

mieQext::usage = "mieQsca[np,nm, a, lambda] recibes as first entry the particle refractive index, then the medium refractive index, the sphere radius and the wavelength. 
					As a result it gives for results: the an coefficient list for first index, the bn coefficient for second index, extinction cross section for third index and
					the extinction efficiency for the fourth index. ";

mieQabs::usage = "mieQsca[np,nm, a, lambda] recibes as first entry the particle refractive index, then the medium refractive index, the sphere radius and the wavelength. 
					As a result it gives for results: the an coefficient list for first index, the bn coefficient for second index, absorption cross section for third index and
					the absorption efficiency for the fourth index. ";


Begin["Private`"];


\[Psi]n[n_,x_]:= x*SphericalBesselJ[n,x]
\[Xi]n[n_,x_]:=x*SphericalHankelH1[n,x]
d\[Psi]n[n_,x_]:=-n *SphericalBesselJ[n,x]+ x * SphericalBesselJ[n-1,x] 
d\[Xi]n[n_,x_]:=-n * SphericalHankelH1[n,x] + x* SphericalHankelH1[n-1,x] 


criterioWiscombe[nm_,a_,lambda_]:=Module[{x,Wiscombe},
						x=(2 Pi nm a)/lambda;							
						Wiscombe=Round[x+ 4*x^(1/3)+2];
						Wiscombe]

mieCoefAn[np_,nm_,a_,lambda_,n_]:=Module[{m,x,an},
						m=np/nm;
						x=(2 Pi nm a)/lambda;							
						an=((m \[Psi]n[n,m x] d\[Psi]n[n,x])-(\[Psi]n[n,x] d\[Psi]n[n,m x]))/((m \[Psi]n[n,m x] d\[Xi]n[n,x])-(\[Xi]n[n,x] d\[Psi]n[n,m x])); 
						an]
mieCoefBn[np_,nm_,a_,lambda_,n_]:=Module[{m,x,bn},
						m=np/nm;
						x=(2 Pi nm a)/lambda;							
						bn=((\[Psi]n[n,m x] d\[Psi]n[n,x])-(m \[Psi]n[n,x] d\[Psi]n[n,m x]))/((\[Psi]n[n,m x] d\[Xi]n[n,x])-(m \[Xi]n[n,x] d\[Psi]n[n,m x]));
						bn]
						
						
						
mieQsca[np_,nm_,a_,lambda_]:=Module[{nmax,an,bn,k,Csca,Qsca},
							nmax=criterioWiscombe[nm,a,lambda];
							an=mieCoefAn[np,nm,a,lambda,#]&/@Range[1,nmax];
							bn=mieCoefBn[np,nm,a,lambda,#]&/@Range[1,nmax];
							Csca= (lambda^2 )/(2*\[Pi]* Norm[nm]^2)*(Sum[(2 n +1)(Norm[an[[n]] ]^2+Norm[bn[[n]]] ^2),{n,nmax}]);
							Qsca=Csca/(Pi* a^2) ;
							{an,bn,Csca,Qsca}]
mieQext[np_,nm_,a_,lambda_]:=Module[{nmax,an,bn,Cext,Qext},
							nmax=criterioWiscombe[nm,a,lambda];
							an=mieCoefAn[np,nm,a,lambda,#]&/@Range[1,nmax];
							bn=mieCoefBn[np,nm,a,lambda,#]&/@Range[1,nmax];
							Cext=(lambda^2)/(2 Pi Norm[nm]^2)*Sum[(2 n+1) (Re[an[[n]]]+Re[bn[[n]]]),{n,1,nmax}];
							Qext=Cext/(Pi a^2);
							{an,bn,Cext,Qext}]
mieQabs[np_,nm_,a_,lambda_]:=mieQext[np,nm,a,lambda]-mieQsca[np,nm,a,lambda]



End[];


EndPackage[];
