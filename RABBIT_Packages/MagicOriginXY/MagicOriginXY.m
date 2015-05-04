(* Mathematica Package *)

(* Created by the Wolfram Workbench Nov 1, 2013 *)

BeginPackage["MagicOriginXY`"]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicOriginXY`*"];
ClearAll @@ Names["MagicOriginXY`*"];

$VersionMagicOriginXY = "MagicOriginXY Version 1.0!";

getPopSizeXY::usage = "getPopSizeXY  "

getCoalProbXY2::usage = "getCoalProbXY2  "

getCoalProbXY3::usage = "getCoalProbXY3  "

(*getfounderGamma12::usage = "getfounderGamma12  "

getfounderGamma123::usage = "getfounderGamma123  "*)

nonIBDProbXY2::usage = "nonIBDProbXY2  "

nonIBDProbXY3::usage = "nonIBDProbXY3  "

inbredProbXY::usage = "inbredProbXY  "

mapExpansionXY::usage = "mapExpansionXY  "

sumJuncDensityXY::usage = "sumJuncDensityXY  "

juncDensityXY::usage = "juncDensityXY  "

origSummaryXY::usage = "origSummaryXY  "

magicOrigPriorXY::usage = "magicOrigPriorXY  "

magicStationaryProbXY::usage = "magicStationaryProbXY  "

magicBasicRateXY::usage = "magicBasicRateXY  "

magicRateMatrixXY::usage = "magicRateMatrixXY  "


Begin["`Private`"]
(* Implementation of the package *)

getPopSizeXY[scheme_, preSize_] :=
    Which[
        StringMatchQ[scheme, "RM1-"|"RM2-"~~ "NE-" | "E-" ~~ DigitCharacter ..],
        ToExpression[Last[StringSplit[scheme, "-"]]],  
        True,
        Switch[scheme,
            "Pairing", preSize/2,
            "Sibling", 2,
            "FullDiallel", (preSize/2)^2,
            (*"HalfDiallel",  preSize*(preSize - 1)/2,*)
            "RM1-NE" | "RM1-E" |"RM2-NE" | "RM2-E", preSize,
            (*"CM-E" | "CPM-E" | "MAI-E", preSize,*)
            _, preSize
        ]
    ]
    
(*return {s^m,s^p};
s^m: the coalescence probability that two Materanlly derived genes 
     come from a single female of the previous generation;
s^p: the coalescence probability that two Pateranlly derived genes 
     come from a single male of the previous generation;
Assuming sex ratio=1, s^m=s^p;     
*)

getCoalProbXY2[scheme_, preSize_,nowSize_] :=
    Module[ {nf0,nm0},
        nf0 = nm0 = preSize/2;
        Which[
           StringMatchQ[scheme, "RM1-NE"~~ "" | "-" ~~DigitCharacter ...],
           {1/nf0,1/nm0},
           StringMatchQ[scheme, "RM2-NE"~~ "" | "-" ~~DigitCharacter ...],
           {1/(nowSize-1)+(nowSize-2)/(nowSize-1) 1/nf0,
               1/(nowSize-1)+(nowSize-2)/(nowSize-1) 1/nm0},
           StringMatchQ[scheme, "RM1-E"~~ "" | "-" ~~DigitCharacter ...],
           1/(nowSize-1) {1,1},      
           StringMatchQ[scheme, "RM2-E"~~ "" | "-" ~~DigitCharacter ...],
           1/(nowSize-1) {1,1},    
           True,
           Switch[scheme,
                "Pairing", {0,0},
                "Sibling", {1,1},
                "FullDiallel", {(nm0-1)/(nowSize-1),(nf0-1)/(nowSize/2-1)},
                _, Beep[];
                   Abort[];
                   {0,0}
           ]
        ]
    ]
    
(*return {q^m,q^p};
q^m: the coalescence probability that three Materanlly derived genes 
     come from a single female of the previous generation;
q^p: the coalescence probability that three Pateranlly derived genes 
     come from a single male of the previous generation;
Assuming sex ratio=1, s^m=s^p;     
*) 
qfun[nowSize_,ngender_] :=
    Module[ {temp},
        temp = 1-2/(nowSize-2)-(nowSize-4)/(nowSize-2) 1/ngender;
        temp*=(nowSize-2)/(nowSize-1) 1/ngender;
        temp+=1/(nowSize-1) (1-1/ngender);
        temp
    ];        
getCoalProbXY3[scheme_, preSize_,nowSize_] :=
    Module[ {nf0,nm0},
        nf0 = nm0 = preSize/2;
        Which[
            StringMatchQ[scheme, "RM1-NE"~~ "" | "-" ~~DigitCharacter ...],
            {1/nf0 (1 - 1/nf0),1/nm0 (1 - 1/nm0)},
            StringMatchQ[scheme, "RM2-NE"~~ "" | "-" ~~DigitCharacter ...],            
            {qfun[nowSize,nf0],qfun[nowSize,nm0]},
            StringMatchQ[scheme, "RM1-E"~~ "" | "-" ~~DigitCharacter ...],
            1/(nowSize-1) {1,1},  
            StringMatchQ[scheme, "RM2-E"~~ "" | "-" ~~DigitCharacter ...],
            1/(nowSize-1) {1,1}, 
            True,
            Switch[scheme,
             "Pairing", {0,0},
             "Sibling", {0,0},
             "FullDiallel",
                 {(nm0-1)/(nowSize-1) (1-(nm0-2)/(nowSize-2)),
                 (nf0-1)/(nowSize/2-1) (1-(nf0-2)/(nowSize/2-2))},     
             _, Beep[];
                Abort[];
             ]
        ]
    ]

(*
getfounderGamma12[isInbredFounder_] :=
    If[ isInbredFounder,
        {0,1, 1, 1},
        {1,1, 1, 1}
    ]
           
getfounderGamma123[isInbredFounder_] :=
    If[ isInbredFounder,
        {0,0,1,1,1, 1},
        {1,1,1,1,1, 1}
    ]    
*)
    
tranMatrix[{sm_,sp_}] :=
    {{0,1/2, 1/2, 0}, 
    {sm/2,(1 - sm)/4, (1 - sm)/2, (1 - sm)/4}, 
    {0,1/2, 1/2, 0}, 
    {0,(1 - sp), 0, 0}}    
            
nonIBDProbXY2[founderGamma12_, coalProb2_] :=
    Module[ {alpha,t},
         (*calculate non-IBD probability 
         alpha^{mp}(12) for two genes in a single individual;
         beta^X(12) for two genes in two distinct individuals, with the parental origins: X = mm, mp, pp*)
        alpha = Table[0, {t, Length[coalProb2]}];
        alpha[[1]] = founderGamma12;
        Do[alpha[[t]] = tranMatrix[coalProb2[[t]]].alpha[[t - 1]], {t, 2, Length[alpha]}];
        Transpose[alpha]
    ]    
        
nonIBDProbXY3[founderGamma123_, coalProb2_, coalProb3_] :=
    Module[ {alpha,matrix,t},
        (*calculate non-IBD probability
        alpha^X (123)  for three genes in two distinct individuals, with the parental origins X=mmp,mpp 
        beta^X(123) for three genes in three distinct individuals, with the parental origins X=mmm,mmp,mpp,ppp by recurrence relation*)
        alpha = Table[0, {t, Length[coalProb2]}];
        alpha[[1]] = founderGamma123;
        matrix[tt_] :=
            Module[ {sm,sp,qm,qp,rm,rp},
                {sm,sp} = coalProb2[[tt]];
                {qm,qp} = coalProb3[[tt]];
                rm = 1-sm-2 qm;
                rp = 1-sp-2 qp;
                {{sm/2,0, (1-sm)/4, (1-sm)/2,(1-sm)/4,0},
                 {0,0, (1-sp)/2,(1-sp)/2,0,0},
                 {qm 3/4, qm 3/4, rm/8,rm 3/8,rm 3/8,rm/8},
                 {sm/2,0, (1-sm)/4, (1-sm)/2,(1-sm)/4,0},
                 {0,0, (1-sp)/2,(1-sp)/2,0,0},
                 {0,0,rp,0,0,0}}
            ];
        Do[alpha[[t]] = matrix[t].alpha[[t - 1]], {t, 2, Length[alpha]}];
        Transpose[alpha]
    ]
    
junc1122XY[founderJK1122_,coalProb2_, Rm_,Rp_] :=
    Module[ {j1122,t},
        (*calculate expected junc density {J^mp(1122),K^mm(1122),K^mp(1122),K^pp(1122)}*)
        j1122 = ConstantArray[0, {Length[coalProb2]}];
        j1122[[1]] = founderJK1122;
        Do[j1122[[t]] = tranMatrix[coalProb2[[t]]].j1122[[t - 1]]
            +{0,coalProb2[[t, 1]] (Rm[[t-1]]+Rp[[t-1]])/4,0,coalProb2[[t, 2]] Rm[[t-1]]}, {t, 2, Length[j1122]}];
        Transpose[j1122]
    ]    

tranMatrix1232[{sm_,sp_}] :=
    {{0,0,1/2, 0, 1/2,0},
     {0,0,1/2, 1/2, 0,0},
     {0,0,(1-sm)/4,1/4, 1/4, (1-sm)/4}, 
     {0,0,1/2, 0, 1/2,0},
     {0,0,1/2, 1/2, 0,0},
     {0,0,(1 - sp),0, 0, 0}}
            
junc1232XY[founderGamma123_,founderJK1232_, coalProb2_,coalProb3_] :=
    Module[ {alpha,j1232, t},
         (*calculate expected junc density J^X(1232) 
         for two haplotypes with the parental origins: X = mm, mp,pm, pp*)
        alpha = nonIBDProbXY3[founderGamma123, coalProb2, coalProb3];
        j1232 = ConstantArray[0, {Length[coalProb2]}];
        j1232[[1]] = founderJK1232;
        Do[j1232[[t]] = tranMatrix1232[coalProb2[[t]]].j1232[[t - 1]]
            +{alpha[[1,t-1]],0,(1-coalProb2[[t,1]]) (alpha[[1,t-1]]+alpha[[2,t-1]])/2,alpha[[1,t-1]],0,0}, {t, 2, Length[j1232]}];
        Transpose[j1232]
    ]    
    
inbredProbXY[founderGamma12_, coalProb2_] :=
    1- First[nonIBDProbXY2[founderGamma12, coalProb2]];

mapExpansionXY[founderGamma12_, founderR_, coalProb2_] :=
    Module[ {amp,Rm,Rp, t},
        amp = First[nonIBDProbXY2[founderGamma12, coalProb2]];
        Rm = Rp = Table[0, {Length[amp]}];
        {Rm[[1]],Rp[[1]]} = founderR;
        Do[
         Rm[[t]] = Rm[[t - 1]]/2 + Rp[[t - 1]]/2 + amp[[t - 1]];
         Rp[[t]] = Rm[[t - 1]], {t, 2, Length[amp]}];
        {Rm,Rp}
    ]

          
sumJuncDensityXY[founderGamma12_, founderR_,founderJK1122_,coalProb2_] :=
    Module[ {Rm,Rp,j1122},
        {Rm,Rp} = mapExpansionXY[founderGamma12, founderR, coalProb2];
        j1122 = junc1122XY[founderJK1122,coalProb2, Rm,Rp];
        (Rm+Rp) - j1122[[1]]
    ]
          
juncDensityXY[founderGamma12_,founderGamma123_, founderR_,founderJK1122_,founderJK1232_,coalProb2_, coalProb3_] :=
    Module[ {Rm,Rp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp,j1232pm},
        {Rm,Rp} = mapExpansionXY[founderGamma12, founderR, coalProb2];
        j1122mp = junc1122XY[founderJK1122,coalProb2, Rm,Rp][[1]];
        {j1232mp,j1232pm} = junc1232XY[founderGamma123,founderJK1232, coalProb2, coalProb3][[;;2]];
        j1213mp = j1232pm;
        j1222mp = (Rm-j1122mp-j1232mp)/2;
        j1211mp = (Rp-j1122mp-j1213mp)/2;
        {j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}
    ]

origSummaryXY[nFounder_, mateScheme_, founderGamma12_, founderGamma123_, founderR_, founderJ1122_, founderJ1232_] :=
    Module[ {nGeneration, popSize, coalProb2, coalProb3, finb, mapR, 
      sumRho, juncRho,t},
        nGeneration = Length[mateScheme] + 1;
        popSize = FoldList[getPopSizeXY[#2, #1] &, nFounder, mateScheme];
        coalProb2 = ConstantArray[0, {nGeneration, 2}];
        coalProb3 = ConstantArray[0, {nGeneration, 2}];
        Do[coalProb2[[t]] = getCoalProbXY2[mateScheme[[t - 1]], popSize[[t - 1]], popSize[[t]]];
           coalProb3[[t]] = getCoalProbXY3[mateScheme[[t - 1]], popSize[[t - 1]], popSize[[t]]], {t, 2, Length[coalProb2]}];
        finb = inbredProbXY[founderGamma12, coalProb2];
        mapR = mapExpansionXY[founderGamma12, founderR, coalProb2];
        sumRho = sumJuncDensityXY[founderGamma12, founderR, founderJ1122, coalProb2];
        (*{j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}*)
        juncRho = juncDensityXY[founderGamma12, founderGamma123, founderR, founderJ1122, founderJ1232, coalProb2, coalProb3];
        {finb, mapR, sumRho, juncRho}
    ]
    
origSummaryXY[nFounder_,mateScheme_] :=
    Module[ {founderGamma12, founderGamma123, founderR, founderJK1122, 
      founderJK1232},
        founderGamma12 = {0, 1, 1, 1};
        founderGamma123 = Join[{0, 0}, {1, 1, 1, 1} If[ nFounder >= 3,
                                                        1,
                                                        0
                                                    ]];
        founderR = Table[0, {2}];
        founderJK1122 = Table[0, {4}];
        founderJK1232 = Table[0, {6}];
        (*{inbred,mapR, sumRho,junction} from origSummary*)
        origSummaryXY[nFounder, mateScheme, founderGamma12, founderGamma123,
          founderR, founderJK1122, founderJK1232]
    ]      

magicOrigPriorXY[nFounder_, mateScheme_]:=
	Module[{temp},    
		temp = origSummaryXY[nFounder, mateScheme];
        (*juncdist={finb_mp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}*)
        N[Prepend[temp[[4, All, -1]], temp[[1, -1]]]]    
	]
	        
magicStationaryProbXY[nFgl_,inbredf_] :=
    Module[ {prob},
        prob = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        prob /. {{i_, i_} :> inbredf/nFgl, {_, _} :> (1 - inbredf)/(nFgl (nFgl - 1))}
    ]

(*{"J1122","J1211","J1213","J1222","J1232"}*)
magicBasicRateXY[nFgl_, inbredf_, juncrho_] :=
    Module[ {r00m, r00p, r01m, r01p, r10m, r10p, r11},
        r11 = If[ inbredf == 0,
                  0,
                  juncrho[[1]]/(inbredf (nFgl - 1))
              ];
        {r01p, r01m} = juncrho[[{2, 4}]]/(1 - inbredf);
        {r00p, r00m} = 
         If[ nFgl > 2,
             juncrho[[{3, 5}]]/((1 - inbredf) (nFgl - 2)),
             {0, 0}
         ];
        {r10p, r10m} = 
         If[ inbredf == 0,
             {0, 0},
             {r01p, r01m} (1 - inbredf)/(inbredf (nFgl - 1))
         ];
        {r00m, r00p, r01m, r01p, r10m, r10p, r11}
    ]
      
magicRateMatrixXY[nFgl_Integer?(#>=2&), baseR:{_?AtomQ, _?AtomQ,_?AtomQ,_?AtomQ, _?AtomQ,_?AtomQ, _?AtomQ}] :=
    Module[ {lab, qq, q00, q11, r00m, r00p, r01m, r01p, r10m, r10p, r11},
        {r00m, r00p, r01m, r01p, r10m, r10p, r11} = baseR;
        lab = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        qq = Outer[List, lab, lab, 1];
        q11 = -(nFgl - 1) ( r11 + r10m+r10p);
        q00 = -(nFgl - 2) (r00m+r00p) - r01m - r01p;
        qq = qq /. {{{x1_, x1_}, {x1_, x1_}} -> q11, 
               {{x1_, x2_}, {x1_,x2_}} -> q00, 
               {{x1_, x1_}, {x2_,x2_}} -> r11, 
               {{x1_, x1_}, {x1_, x2_}} -> r10p, {{x1_, x1_}, {x2_, x1_}} -> r10m, 
               {{x1_, x2_}, {x1_, x1_}} -> r01p, {{x1_, x2_}, {x2_, x2_}} -> r01m,  
               {{x1_, x2_}, {x1_, x3_}} -> r00p, {{x1_, x2_}, {x3_, x2_}} -> r00m, {{_, _}, {_, _}} -> 0};
        qq
    ]
    
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicOriginXY`*"];

EndPackage[]

