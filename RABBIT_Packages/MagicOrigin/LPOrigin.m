(* Mathematica Package *)

BeginPackage["MagicOrigin`LPOrigin`"]

lpnonIBDProb2::usage = "lpnonIBDProb2  "

lpnonIBDProb3::usage = "lpnonIBDProb3  "

lpmapExpansion::usage = "lpmapExpansion  "

lpjunc1122::usage = "lpjunc1122  "

lpjunc1232::usage = "lpjunc1232  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

lpnonIBDProb2[inialpha12_, coalProb_, t_] :=
    Module[ {lam1},
        lam1 = 1-coalProb/2;
        inialpha12 lam1^t
    ]

lpnonIBDProb3[inialpha123_, coalProb_, t_] :=
    Module[ {lam3},
        lam3 = 1-coalProb 3/2;
        inialpha123 lam3^t
    ]
    
lpmapExpansion[inialpha12_, iniR_, coalProb_, t_] :=
    Module[ {lam1},
        lam1 = 1-coalProb/2;
        iniR+inialpha12 (1-lam1^t)/(1-lam1)
    ]     
    
lpjunc1122[inialpha12_, iniR_, iniJ1122_, coalProb_, t_] :=
    Module[ {lam1},
        lam1 = 1-coalProb/2;
        iniJ1122 lam1^t+iniR (1-lam1^t) +inialpha12/(1-lam1) (1-(1+(1/lam1-1) t) lam1^t)
    ]    

lpjunc1232[inialpha123_, iniJ1232_, coalProb_, t_]:=
    Module[ {lam1,lam3},
    	lam1 = 1-coalProb/2;
    	lam3 = 1-coalProb 3/2;
    	iniJ1232 lam1^t+inialpha123 (2 lam1-1)/(lam1-lam3) (lam1^t-lam3^t)
    ]
        
End[] (* End Private Context *)

EndPackage[]