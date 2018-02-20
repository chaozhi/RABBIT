(* Mathematica Package *)

BeginPackage["MaginOrigin`DOOrigin`", {"MagicOrigin`","MagicOrigin`LPOrigin`","MagicOrigin`SIBOrigin`"}]

inbredProbDO::usage = "inbredProbDO  "

mapExpansionDO::usage = "mapExpansionDO  "

juncDensityDO::usage = "juncDensityDO  "

sumJuncDensityDO::usage = "sumJuncDensityDO  "

origSummaryDO::usage = "origSummaryDO  "

magicOrigPriorDO::usage = "magicOrigPriorDO  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*subscript 0 refers to the founder population t=0;
  The pre-CC lines are assumed to be independent, so that 
  K_0(1122)=0; K_0(1232) =R_0;
 *)
 
(*initial conditions refer to the F1 populations of DO*)
(*inialpha12 ~ beta_0(12)*)
(*inialpha123  ~ beta_0(123)*)    
(*iniR =R_0+alpha_0(12)*)
(*iniJK1122~0*)
(*iniJK1232 ~ K_0(1232)+alpha_0(123); 
*)

getinialpha12[nStrain_] :=
    1 - 1/nStrain  
    
getinialpha123[nStrain_] :=
    (1 - 1/nStrain) (1 - 2/nStrain)

getiniR[founderAlpha12_,founderR_] :=
    founderR +founderAlpha12
    
getiniJK1122[] :=
    0
       
getiniJK1232[founderAlpha12_,founderR_,nStrain_] :=
    founderR (1-2/nStrain)+founderAlpha12 (1 - 2/nStrain)

inbredProbDO[founderAlpha12_,nStrain_,coalProb_,gCross_] :=
    Module[ {inialpha12,a12,t},
        inialpha12 = getinialpha12[nStrain];
        a12 = Table[Evaluate[lpnonIBDProb2[inialpha12,coalProb,t]],{t,gCross}];
        Join[{1-founderAlpha12,1-inialpha12},1 - a12]
    ]    

mapExpansionDO[founderAlpha12_,founderR_,nStrain_,coalProb_,gCross_] :=
    Module[ {inialpha12,iniR,Ra,t},
        inialpha12 = getinialpha12[nStrain];
        iniR = getiniR[founderAlpha12,founderR];
        Ra = Table[Evaluate[lpmapExpansion[inialpha12, iniR, coalProb,t]],{t,gCross}];
        Join[{founderR,iniR},Ra]
    ]       
    
junc1122DO[founderAlpha12_,founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {founderJ1122, inialpha12,iniR, iniJK1122,j1122,t},
        founderJ1122 = 0;
        inialpha12 = getinialpha12[nStrain];
        iniR = getiniR[founderAlpha12,founderR];
        iniJK1122 = getiniJK1122[];
        j1122 = Table[Evaluate[lpjunc1122[inialpha12, iniR, iniJK1122, coalProb, t]],{t,gCross}];
        Join[{founderJ1122,iniJK1122}, j1122]
    ]   
    
junc1232DO[founderAlpha12_,founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {founderJ1232, inialpha12, inialpha123,iniJK1232,j1232,t},
        founderJ1232 = founderR (1-2/nStrain);
        inialpha12 = getinialpha12[nStrain];
        inialpha123 = getinialpha123[nStrain];
        iniJK1232 = getiniJK1232[founderAlpha12,founderR,nStrain];
        j1232 = Table[Evaluate[lpjunc1232[inialpha123, iniJK1232,coalProb,t]],{t,gCross}];
        Join[{founderJ1232,iniJK1232},j1232]
    ]    

juncDensityDO[founderAlpha12_,founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {Ra,j1122,j1222,j1232},
        Ra = mapExpansionDO[founderAlpha12,founderR,nStrain,coalProb,gCross];
        j1122 = junc1122DO[founderAlpha12,founderR, nStrain,coalProb, gCross];
        j1232 = junc1232DO[founderAlpha12,founderR,nStrain,coalProb, gCross];
        j1222 = (Ra - j1232 - j1122)/2;
        {j1122,j1222,j1232}
    ]

sumJuncDensityDO[founderAlpha12_,founderR_, nStrain_,coalProb_, gCross_] :=
    Module[ {Ra,j1122},
        Ra = mapExpansionDO[founderAlpha12,founderR,nStrain,coalProb,gCross];
        j1122 = junc1122DO[founderAlpha12,founderR,nStrain,coalProb, gCross];
        2 Ra - j1122
    ]

origSummaryDO[nPower_, preCCfreq_, crossPopSize_, gCross_, crossScheme_] :=
    Module[ {inigamma12, iniR, freq = preCCfreq, nStrain, coalProb, 
      founderalpha12, founderR, finb, mapR,sumRho, juncRho},
        inigamma12 = {1, 1};
        iniR = nPower - 2;
        {inigamma12,iniR}=N[{inigamma12,iniR}];
        freq[[All, 2]] = N[Normalize[freq[[All, 2]], Total]];
        founderalpha12 = Total[#[[2]] sibIBDProb2[inigamma12, #[[1]]][[1]] & /@ freq];
        founderR = Total[#[[2]] sibmapExpansion[inigamma12, iniR, #[[1]]] & /@ freq];
        coalProb = N[getCoalProb2[crossScheme, crossPopSize]];
        nStrain = 2^nPower;
        finb = inbredProbDO[founderalpha12, nStrain, coalProb, gCross];
        mapR = mapExpansionDO[founderalpha12, founderR, nStrain, coalProb, gCross];
        sumRho = sumJuncDensityDO[founderalpha12, founderR, nStrain, coalProb,gCross];
        juncRho = juncDensityDO[founderalpha12, founderR, nStrain, coalProb, gCross];
        {finb, mapR, sumRho, juncRho}
    ]
    
magicOrigPriorDO[nPower_, preCCfreq_, crossPopSize_, gCross_, crossScheme_]:=
	Module[{temp,inbred,j1122,j1211,j1213,j1222,j1232},    
		temp = origSummaryDO[nPower, preCCfreq, crossPopSize, gCross, crossScheme];
		{inbred,j1122,j1222,j1232}=N[Prepend[temp[[4, All, -1]], temp[[1, -1]]]];
        (*{inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp} *)
        {j1211,j1213}={j1222,j1232};
        {inbred,j1122,j1211,j1213,j1222,j1232}
	]    
  
End[] (* End Private Context *)

EndPackage[]