(* Mathematica Package *)

BeginPackage["MagicOriginXY`AILOriginXY`",{"MagicOriginXY`LPOriginXY`"}]

inbredProbAILXY::usage = "inbredProbAILXY  "

mapExpansionAILXY::usage = "mapExpansionAILXY  "

mapExpansionAILXY2::usage = "mapExpansionAILXY2  "

juncDensityAILXY::usage = "juncDensityAILXY  "

sumJuncDensityAILXY::usage = "sumJuncDensityAILXY  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*There are nFounder/2 inbred females and nFounder/2 males;
Initial population refers to the F1 of AIL*)

(*inibeta12 ={beta^mm_0(12),beta^mp_0(12),beta^pp_0(12)}*)
getinibeta12[nFounder_] := {1 - 2/ nFounder, 1, 1 - 2/nFounder}
    

(*inibeta123 ={beta^mmm_0(123),beta^mmp_0(123),beta^mpp_0(123),beta^ppp_0(123)}*)    
getinibeta123[nFounder_] := {(1 - 2/nFounder)(1- 4/nFounder), 1 - 2/nFounder, 1 - 2/nFounder, (1 - 2/nFounder)(1- 4/nFounder)}

getiniR[]:={0,0}

(*iniJ1122 ={J1122^mm_0,J1122^mp_0,J1122^pp_0}*)
getiniJ1122[]:={0,0,0}

(*iniJ1232 ={J1232^mm_0,J1232^avg_0,J1232^pp_0}*)
getiniJ1232[]:={0,0,0}

getiniJ1232diff[]:=0
 
inbredProbAILXY[nFounder_,coalProb_,gCross_] :=
    Module[ {inibeta12,a12,t},        
        inibeta12 = getinibeta12[nFounder];
        a12 = Table[Evaluate[lpnonIBDProbXY2[inibeta12,coalProb,t][[2]]],{t,gCross}];
        Join[{1,1-inibeta12[[2]]},1 - a12]
    ]    

mapExpansionAILXY[nFounder_,coalProb_,gCross_] :=
    Module[ {inibeta12,iniR,Rm,Rp,t},
        inibeta12 = getinibeta12[nFounder];
        iniR = getiniR[];
        {Rm, Rp} = Transpose[Table[Evaluate[lpmapExpansionXY[inibeta12, iniR, coalProb,t]],{t,gCross}]];
        {Join[{0,0},Rm],Join[{0,0},Rp]}
    ]   
    
mapExpansionAILXY2[nFounder_, coalProb_,gCross_] :=
    Module[ {inibeta12,iniR,Rx,Rd,t},
        inibeta12 = getinibeta12[nFounder];
        iniR = getiniR[];  
        {Rx, Rd} = Transpose[Table[Evaluate[lpmapExpansionXY2[inibeta12, iniR, coalProb,t]],{t,gCross}]];
        {Join[{0,0},Rx],Join[{0,0},Rd]}
    ]           
    
junc1122AILXY[nFounder_, coalProb_, gCross_] :=
    Module[ {inibeta12,iniR,iniJ1122,j1122mp,t},
        inibeta12 = getinibeta12[nFounder];
        iniR = getiniR[];
        iniJ1122 = getiniJ1122[];
        j1122mp = Table[Evaluate[lpjunc1122XY[inibeta12, iniR, iniJ1122, coalProb, t][[2]]],{t,gCross}];
        Join[{0,0}, j1122mp]
    ]   
    
junc1232AILXY[nFounder_, coalProb_, gCross_] :=
    Module[ {inibeta12, inibeta123, iniR, iniJ1232,iniJ1232diff,j1232avg,j1232diff,j1232mp,j1232pm,t},        
        inibeta12 = getinibeta12[nFounder];
        inibeta123 = getinibeta123[nFounder];
        iniR = getiniR[];
        iniJ1232=getiniJ1232[];
        iniJ1232diff= getiniJ1232diff[];    
        j1232avg = Table[Evaluate[lpjunc1232XYsym[inibeta123, iniJ1232,coalProb,t][[2]]],{t,gCross}];
        j1232diff=Table[Evaluate[lpjunc1232XYdiff[inibeta123, iniJ1232diff, coalProb, t]],{t,gCross}];        
        j1232mp=j1232avg+j1232diff;
        j1232pm=j1232avg-j1232diff;
        {Join[{0,0},j1232mp],Join[{0,0},j1232pm]}
    ]         

juncDensityAILXY[nFounder_, coalProb_, gCross_] :=
    Module[ {Rm,Rp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {Rm,Rp} = mapExpansionAILXY[nFounder,coalProb,gCross];
        j1122mp = junc1122AILXY[nFounder,coalProb,gCross];
        {j1232mp,j1213mp} = junc1232AILXY[nFounder, coalProb, gCross];        
        j1222mp = (Rm-j1122mp-j1232mp)/2;
        j1211mp = (Rp-j1122mp-j1213mp)/2;
        {j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}
    ]

sumJuncDensityAILXY[nFounder_, coalProb_, gCross_] :=
    Module[ {Rm,Rp,j1122mp},
        {Rm,Rp} = mapExpansionAILXY[nFounder,coalProb,gCross];
        j1122mp = junc1122AILXY[nFounder,coalProb,gCross];
        Rm+Rp - j1122mp
    ]
        

End[] (* End Private Context *)

EndPackage[]