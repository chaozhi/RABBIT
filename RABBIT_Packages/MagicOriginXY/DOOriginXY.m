(* Mathematica Package *)

BeginPackage["MagicOriginXY`DOOriginXY`", {"MagicOriginXY`","MagicOriginXY`PairingOriginXY`","MagicOriginXY`SIBOriginXY`","MagicOriginXY`LPOriginXY`"}]

inbredProbDOXY::usage = "inbredProbDOXY  "

mapExpansionDOXY::usage = "mapExpansionDOXY  "

mapExpansionDOXY2::usage = "mapExpansionDOXY2  "

juncDensityDOXY::usage = "juncDensityDOXY  "

sumJuncDensityDOXY::usage = "sumJuncDensityDOXY  "

origSummaryDOXY::usage = "origSummaryDOXY  "

magicOrigPriorDOXY::usage = "magicOrigPriorDOXY  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*DO founder population consists of N/2 females and N/2 males, sampled from independent funnels of CC populations*)

(*Subscript 0 refer to the founder population t=0;
  The pre-CC lines are assumed to be independent, so that 
  K^x_0(1122)=0 for x=mm,mp,pp;
  K^{mp}_0(1232) =R^m_0; K^{pm}_0(1232) =R^p_0; *) 
(*Initial population refers to the F1 of DO, t=1*)
(*inibeta12 ={beta^mm_0(12),beta^mp_0(12),beta^pp_0(12)}*) 

getinibeta12[nStrain_] :=
    { 1 - 1/nStrain, 1 - 1/nStrain, 1 - 1/nStrain}

(*inibeta123 ={beta^mmm_0(123),beta^mmp_0(123),beta^mpp_0(123),beta^ppp_0(123)}*)  
getinibeta123[nStrain_] :=
    {1, 1, 1, 1} (1 - 1/nStrain) (1 - 2/nStrain)

(*founderAlpha12mp: alpha^mp_0(12) in founder population at generation 0*) 
getiniR[founderAlpha12mp_,founderR_] :=
    Module[ {Rm0,Rp0},
        {Rm0,Rp0} = founderR;
        {(Rm0 + Rp0)/2 + founderAlpha12mp, Rm0}
    ]



    
(*iniK1122 ={K1122^mm_0,K1122^mp_0,K1122^pp_0}*)
getiniJ1122[] :=
    {0,0,0}

(*founder K1232 at generation 0: {K1232mm0,K1232mp0,K1232pm0,K1232pp0}={Rm0,Rm0,Rp0,Rp0}*)
(*initial K1232 at generation 1*)   
(*iniJ1232 ={{J1232^mm_0,J1232^avg_0,J1232^pp_0},J1232^avg_diff}*)
getiniJ1232[founderAlpha12mp_,founderR_,nStrain_] :=
    Module[ {Rm0,Rp0,jj},
        {Rm0,Rp0} = founderR;
        jj = {((Rm0 + Rp0)/2 + founderAlpha12mp), 
            ((Rm0 + Rp0)/2 + founderAlpha12mp), Rm0, Rm0}  (1 - 2/nStrain);
        {{jj[[1]],(jj[[2]]+jj[[3]])/2,jj[[4]]},(jj[[2]]-jj[[3]])/2}
    ] 

inbredProbDOXY[founderAlpha12mp_,nStrain_,coalProb_,gCross_] :=
    Module[ {inibeta12,a12,t},
        inibeta12 = getinibeta12[nStrain];
        a12 = Table[Evaluate[lpnonIBDProbXY2[inibeta12,coalProb,t][[2]]],{t,gCross}];
        Join[{1-founderAlpha12mp,1-inibeta12[[2]]},1 - a12]
    ]    

mapExpansionDOXY[founderAlpha12mp_,founderR_,nStrain_,coalProb_,gCross_] :=
    Module[ {inibeta12,iniR,Rm,Rp,t},
        inibeta12 = getinibeta12[nStrain];
        iniR = getiniR[founderAlpha12mp,founderR];
        {Rm, Rp} = Transpose[Table[Evaluate[lpmapExpansionXY[inibeta12, iniR, coalProb,t]],{t,gCross}]];
        {Join[{founderR[[1]],iniR[[1]]},Rm],Join[{founderR[[2]],iniR[[2]]},Rp]}
    ]   
    
mapExpansionDOXY2[founderAlpha12mp_, founderR_,nStrain_,coalProb_,gCross_] :=
    Module[ {inibeta12,iniR,Rx,Rd,t},
        inibeta12 = getinibeta12[nStrain];
        iniR = getiniR[founderAlpha12mp,founderR];
        {Rx, Rd} = Transpose[Table[Evaluate[lpmapExpansionXY2[inibeta12, iniR, coalProb,t]],{t,gCross}]];
        {Join[{2/3,1/3}.#&/@{founderR,iniR},Rx],Join[{1/2,-1/2}.#&/@{founderR,iniR},Rd]}
    ]           
 
junc1122DOXY[founderAlpha12mp_, founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {founderJ1122mp,inibeta12,iniR,iniJ1122,j1122mp,t},
        founderJ1122mp = 0;
        inibeta12 = getinibeta12[nStrain];
        iniR = getiniR[founderAlpha12mp,founderR];
        iniJ1122 = getiniJ1122[];
        j1122mp = Table[Evaluate[lpjunc1122XY[inibeta12, iniR, iniJ1122, coalProb, t][[2]]],{t,gCross}];
        Join[{founderJ1122mp,iniJ1122[[2]]}, j1122mp]
    ]   
    
junc1232DOXY[founderAlpha12mp_, founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {founderJ1232mp,founderJ1232pm,inibeta123, iniJ1232,iniJ1232diff,j1232avg,j1232diff,j1232mp,j1232pm,t},
        {founderJ1232mp,founderJ1232pm} = founderR (1-2/nStrain);
        inibeta123 = getinibeta123[nStrain];
        {iniJ1232,iniJ1232diff} = getiniJ1232[founderAlpha12mp,founderR, nStrain];
        j1232avg = Table[Evaluate[lpjunc1232XYsym[inibeta123, iniJ1232,coalProb,t][[2]]],{t,gCross}];
        j1232diff = Table[Evaluate[lpjunc1232XYdiff[inibeta123, iniJ1232diff, coalProb, t]],{t,gCross}];
        j1232mp = j1232avg+j1232diff;
        j1232pm = j1232avg-j1232diff;
        {Join[{founderJ1232mp,iniJ1232[[2]]+iniJ1232diff},j1232mp],
         Join[{founderJ1232pm,iniJ1232[[2]]-iniJ1232diff},j1232pm]}
    ]         

juncDensityDOXY[founderAlpha12mp_,founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {Rm,Rp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {Rm,Rp} = mapExpansionDOXY[founderAlpha12mp,founderR,nStrain,coalProb,gCross];
        j1122mp = junc1122DOXY[founderAlpha12mp, founderR,nStrain,coalProb, gCross];
        {j1232mp,j1213mp} = junc1232DOXY[founderAlpha12mp, founderR,nStrain,coalProb, gCross];
        j1222mp = (Rm-j1122mp-j1232mp)/2;
        j1211mp = (Rp-j1122mp-j1213mp)/2;
        {j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}
    ]

sumJuncDensityDOXY[founderAlpha12mp_,founderR_,nStrain_,coalProb_, gCross_] :=
    Module[ {Rm,Rp,j1122mp},
        {Rm,Rp} = mapExpansionDOXY[founderAlpha12mp,founderR,nStrain,coalProb,gCross];
        j1122mp = junc1122DOXY[founderAlpha12mp, founderR,nStrain,coalProb, gCross];
        Rm+Rp - j1122mp
    ]       
    
origSummaryDOXY[nPower_, preCCfreq_, crossPopSize_, gCross_, crossScheme_] :=
    Module[ {freq = preCCfreq, inibeta12, iniR, founderAlpha12mp, 
      founderR, nStrain, coalProb, finb, mapR, sumRho, juncRho},
        freq[[All, 2]] = N[Normalize[freq[[All, 2]], Total]];
        inibeta12 = Take[pairIBDProbXY2[], 2];
        iniR = pairmapExpansionXY[nPower - 2];
        {inibeta12,iniR}=N[{inibeta12,iniR}];
        founderAlpha12mp = Total[#[[2]] sibIBDProbXY2[inibeta12, #[[1]]][[2]] & /@ freq];
        founderR = Total[#[[2]] sibmapExpansionXY[inibeta12, iniR, #[[1]]] & /@ freq];
        nStrain = 2^nPower;
        coalProb = N[First[getCoalProbXY2[crossScheme, crossPopSize, crossPopSize]]];
        finb = inbredProbDOXY[founderAlpha12mp, nStrain, coalProb, gCross];
        mapR = mapExpansionDOXY[founderAlpha12mp, founderR, nStrain, coalProb, gCross];
        sumRho = sumJuncDensityDOXY[founderAlpha12mp, founderR, nStrain, coalProb,gCross];
        juncRho = juncDensityDOXY[founderAlpha12mp, founderR, nStrain, coalProb,gCross];
        {finb, mapR, sumRho, juncRho}
    ]    

magicOrigPriorDOXY[nPower_, preCCfreq_, crossPopSize_, gCross_, crossScheme_]:=
	Module[{temp},    
		temp = origSummaryDOXY[nPower, preCCfreq, crossPopSize, gCross, crossScheme];
        (*juncdist={inbred,j1122,j1222,j1232}*)
        N[Prepend[temp[[4, All, -1]], temp[[1, -1]]]]    
	]
	
End[] (* End Private Context *)

EndPackage[]