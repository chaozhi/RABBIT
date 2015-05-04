(* Mathematica Package *)

BeginPackage["MagicOrigin`AILOrigin`",{"MagicOrigin`LPOrigin`"}]

inbredProbAIL::usage = "inbredProbAIL  "

mapExpansionAIL::usage = "mapExpansionAIL  "

juncDensityAIL::usage = "juncDensityAIL  "

sumJuncDensityAIL::usage = "sumJuncDensityAIL  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*initial conditions refer to the F2 populations of AIL*)
(*inialpha12 =alpha_0(12) ~ beta_0(12)*)
(*inialpha123 =alpha_0(123)~ beta_0(123)}*)    
getinialpha12[nFounder_] :=
    1 - 1/ nFounder

getinialpha123[nFounder_] :=
    (1 - 1/ nFounder) (1 - 2/ nFounder)

getiniR[] :=
    1      
    
getiniJ1122[] :=
    0

getiniJ1232[nFounder_] :=
    1-2/nFounder
 
inbredProbAIL[nFounder_,coalProb_,gCross_] :=
    Module[ {inialpha12,a12,t},
        inialpha12 = getinialpha12[nFounder];
        a12 = Table[Evaluate[lpnonIBDProb2[inialpha12,coalProb,t]],{t,gCross-1}];
        Join[{1,0,1-inialpha12},1 - a12]
    ]    

mapExpansionAIL[nFounder_,coalProb_,gCross_] :=
    Module[ {inialpha12,iniR,Ra,t},
        inialpha12 = getinialpha12[nFounder];
        iniR = getiniR[];
        Ra = Table[Evaluate[lpmapExpansion[inialpha12, iniR, coalProb,t]],{t,gCross-1}];
        Join[{0,0,iniR},Ra]
    ]       
    
junc1122AIL[nFounder_, coalProb_, gCross_] :=
    Module[ {inialpha12,iniR, iniJ1122,j1122,t},
        inialpha12 = getinialpha12[nFounder];
        iniR = getiniR[];
        iniJ1122 = getiniJ1122[];
        j1122 = Table[Evaluate[lpjunc1122[inialpha12, iniR, iniJ1122, coalProb, t]],{t,gCross-1}];
        Join[{0,0,iniJ1122}, j1122]
    ]   
    
junc1232AIL[nFounder_, coalProb_, gCross_] :=
    Module[ {inialpha12, inialpha123,iniJ1232,j1232,t},
        inialpha12 = getinialpha12[nFounder];
        inialpha123 = getinialpha123[nFounder];
        iniJ1232 = getiniJ1232[nFounder];
        j1232 = Table[Evaluate[lpjunc1232[inialpha123, iniJ1232,coalProb,t]],{t,gCross-1}];
        Join[{0,0,iniJ1232},j1232]
    ]    

juncDensityAIL[nFounder_, coalProb_, gCross_] :=
    Module[ {Ra,j1122,j1222,j1232},
        Ra = mapExpansionAIL[nFounder,coalProb,gCross];
        j1122 = junc1122AIL[nFounder,coalProb,gCross];
        j1232 = junc1232AIL[nFounder, coalProb, gCross];
        j1222 = (Ra - j1232 - j1122)/2;
        {j1122,j1222,j1232}
    ]

sumJuncDensityAIL[nFounder_, coalProb_, gCross_] :=
    Module[ {Ra,j1122},
        Ra = mapExpansionAIL[nFounder,coalProb,gCross];
        j1122 = junc1122AIL[nFounder,coalProb,gCross];
        2 Ra - j1122
    ]
        
End[] (* End Private Context *)

EndPackage[]