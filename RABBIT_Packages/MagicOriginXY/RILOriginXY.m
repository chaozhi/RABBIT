(* Mathematica Package *)

BeginPackage["MagicOriginXY`RILOriginXY`",{"MagicOriginXY`SIBOriginXY`"}]

inbredProbRILXY::usage = "inbredProbRILXY  "

mapExpansionRILXY::usage = "mapExpansionRILXY  "

mapExpansionRILXY2::usage = "mapExpansionRILXY2  "

juncDensityRILXY::usage = "juncDensityRILXY  "

sumJuncDensityRILXY::usage = "sumJuncDensityRILXY  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

getgCross[nPower_] :=
    Max[nPower - 2, 0]   
  
(*inibeta12 ={beta^mm_0(12),beta^mp_0(12)}*)
getinibeta12[nPower_] :=
    {Boole[nPower > 1], 1}

(*inialpha123 =alpha^mmp_0(123)*)    
getinialpha123[nPower_] := Boole[nPower > 1]
      
inbredProbRILXY[nPower_, gInbred_] :=
    Module[ {gCross, inibeta12,a12,t},
        gCross = getgCross[nPower];
        inibeta12 = getinibeta12[nPower];
        a12 = Table[Evaluate[sibIBDProbXY2[inibeta12, t][[2]]],{t,0,gInbred}];
        Join[{1}, Table[0, {gCross}], 1 - a12]
    ]    
    
getRCross[nPower_] :=
    Module[ {gCross,ls,RmCross,RpCross},
        gCross = getgCross[nPower];
        ls = Range[0, gCross];
        RmCross = 2 (1 - (-2)^(-ls) + 3 ls)/9;
        RpCross = 4 (-1 + (-2)^(-ls) + 3/2 ls)/9;
        {RmCross,RpCross}
    ]

mapExpansionRILXY[nPower_, gInbred_] :=
    Module[ {gCross, RmCross,RpCross,inibeta12,iniR,Rm,Rp,t},
        gCross = getgCross[nPower];
        inibeta12 = getinibeta12[nPower];
        {RmCross,RpCross} = getRCross[nPower];
        iniR = Last[#]&/@{RmCross,RpCross};
        {Rm, Rp} = Transpose[Table[Evaluate[sibmapExpansionXY[inibeta12, iniR, t]],{t,1,gInbred}]];
        {Join[{0}, RmCross, Rm],Join[{0}, RpCross, Rp]}
    ]   
    
mapExpansionRILXY2[nPower_, gInbred_] :=
    Module[ {gCross, RmCross,RpCross,RxCross,RdCross,inibeta12,iniR,Rx,Rd,t},
        gCross = getgCross[nPower];
        inibeta12 = getinibeta12[nPower];
        {RmCross,RpCross} = getRCross[nPower];
        iniR = Last[#]&/@{RmCross,RpCross};     
        {Rx, Rd} = Transpose[Table[Evaluate[sibmapExpansionXY2[inibeta12, iniR, t]],{t,1,gInbred}]];
        RxCross =2/3 RmCross+1/3 RpCross;
        RdCross =1/2 RmCross -1/2 RpCross;
        {Join[{0}, RxCross, Rx],Join[{0}, RdCross, Rd]}
    ]           
    
junc1122RILXY[nPower_, gInbred_] :=
    Module[ {gCross, RmCross,RpCross, inibeta12,iniR,iniJ1122,j1122mp,t},
        gCross = getgCross[nPower];
        inibeta12 = getinibeta12[nPower];
        {RmCross,RpCross} = getRCross[nPower];
        iniR = Last[#]&/@{RmCross,RpCross};
        (*iniJ1122 ={J1122^mm_0,J1122^mp_0}*)
        iniJ1122 = {0,0};
        j1122mp = Table[Evaluate[sibjunc1122XY[inibeta12, iniR, iniJ1122, t][[2]]],{t,1,gInbred}];
        Join[{0},Table[0, {gCross + 1}], j1122mp]
    ]   
    
junc1232RILXY[nPower_, gInbred_] :=
    Module[ {gCross, RmCross,RpCross,inibeta12, inialpha123, iniR, iniJ1232,iniJ1232diff,j1232avg,j1232diff,j1232mp,j1232pm,t},
        gCross = getgCross[nPower];
        {RmCross,RpCross} = getRCross[nPower];
        gCross = getgCross[nPower];
        inibeta12 = getinibeta12[nPower];
        inialpha123 = getinialpha123[nPower];
        {RmCross,RpCross} = getRCross[nPower];
        iniR = Last[#]&/@{RmCross,RpCross};
        (*iniJ1232 ={J1232^mm_0,J1232^avg_0}*)   
        iniJ1232={iniR[[1]],(iniR[[1]]+iniR[[2]])/2};
        iniJ1232diff= (iniR[[1]]-iniR[[2]])/2;          
        j1232avg = Table[Evaluate[sibjunc1232XYsym[inialpha123,iniJ1232, t][[2]]],{t,1,gInbred}];
        j1232diff=Table[Evaluate[sibjunc1232XYdiff[inialpha123, iniJ1232diff, t]],{t,1,gInbred}];        
        j1232mp=j1232avg+j1232diff;
        j1232pm=j1232avg-j1232diff;
        {Join[{0},RmCross, j1232mp],Join[{0},RpCross, j1232pm]}
    ]         

juncDensityRILXY[nPower_, gInbred_] :=
    Module[ {Rm,Rp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {Rm,Rp} = mapExpansionRILXY[nPower,gInbred];
        j1122mp = junc1122RILXY[nPower,gInbred];
        {j1232mp,j1213mp} = junc1232RILXY[nPower, gInbred];        
        j1222mp = (Rm-j1122mp-j1232mp)/2;
        j1211mp = (Rp-j1122mp-j1213mp)/2;
        {j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}
    ]

sumJuncDensityRILXY[nPower_, gInbred_] :=
    Module[ {Rm,Rp,j1122mp},
        {Rm,Rp} = mapExpansionRILXY[nPower,gInbred];
        j1122mp = junc1122RILXY[nPower,gInbred];
        Rm+Rp - j1122mp
    ]
        
End[] (* End Private Context *)

EndPackage[]