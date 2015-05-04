(* Mathematica Package *)

(* Created by the Wolfram Workbench Oct 30, 2013 *)

BeginPackage["MagicOrigin`"]

Unprotect @@ Names["MagicOrigin`*"];
ClearAll @@ Names["MagicOrigin`*"];

getSelfing::usage = "getSelfing  "

getPopSize::usage = "getPopSize  "

getCoalProb2::usage = "getCoalProb2  "

getCoalProb3::usage = "getCoalProb3  "

(*getfounderGamma12::usage = "getfounderGamma12  "

getfounderGamma123::usage = "getfounderGamma123  "*)

nonIBDProb2::usage = "nonIBDProb2  "

nonIBDProb3::usage = "nonIBDProb3  "

inbredProb::usage = "inbredProb  "

mapExpansion::usage = "mapExpansion  "

juncDensity::usage = "juncDensity  "

sumJuncDensity::usage = "sumJuncDensity  "

origSummary::usage = "origSummary  "

magicOrigPrior::usage = "magicOrigPrior  "

magicStationaryProb::usage = "magicStationaryProb  "

magicRateMatrix::usage = "magicRateMatrix  "

magicBasicRate::usage = "magicBasicRate  "

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package *)

getSelfing::wrongscheme :=
    "Wrong mating scheme `1`!"

getSelfing[scheme_] :=
    Which[
        StringMatchQ[scheme, "Selfing" ~~ DigitCharacter .. ~~ "-" ~~ DigitCharacter ..],
        True,
        StringMatchQ[scheme, "WF-"~~ "NE" | "E" ~~ "" | "-" ~~DigitCharacter ...],
        True,
        StringMatchQ[scheme, "RM1-" | "RM11-" | "RM2-" ~~ "NE" | "E" ~~ "" | "-" ~~DigitCharacter ...],
        False,    
        True,    
        Switch[scheme,
            "FounderPopulation", False,
            "Selfing", True,
            "Pairing"|"Sibling"|"FullDiallel"|"HalfDiallel",False,
            "CM-E" | "CPM-E" | "MAI-E", False,
            _, Print["getSelfing::Wrong mating scheme ",scheme];
               Beep[];
               Abort[];
        ]
    ]
 
getPopSize[scheme_, preSize_] :=
    Which[
        StringMatchQ[scheme, "Selfing" ~~ DigitCharacter .. ~~ "-" ~~ DigitCharacter ..],
        1,
        StringMatchQ[scheme, "WF-"|"RM1-"|"RM11-" |"RM2-"~~ "NE-" | "E-" ~~ DigitCharacter ..],
        ToExpression[Last[StringSplit[scheme, "-"]]],  
        True,
        Switch[scheme,
            "Pairing", preSize/2,
            "Sibling", 2,
            "Selfing", 1,            
            "FullDiallel", preSize*(preSize - 1),
            "FullDiallel2", 2 (preSize/2)^2,            
            "HalfDiallel", preSize*(preSize - 1)/2,
            "HalfDiallel2", (preSize/2)^2,
            "WF-NE" | "WF-E" | "RM1-NE" | "RM1-E" | "RM11-NE" | "RM11-E" | 
            "RM2-NE" | "RM2-E" | "CM-E" | "CPM-E" | "MAI-E", preSize,
            _, Print["getPopSize::Wrong mating scheme ",scheme];
               Beep[];
               Abort[];
        ]
    ]   
                      
getCoalProb2[scheme_, preSize_] :=
	Module[{nowsize},
    Which[
       StringMatchQ[scheme, "Selfing" ~~ DigitCharacter .. ~~ "-" ~~ DigitCharacter ..],
       1,
       StringMatchQ[scheme, "WF-NE"|"RM1-NE"|"RM11-NE"~~ "" | "-" ~~DigitCharacter ...],
       1/preSize,
       StringMatchQ[scheme, "RM2-NE"~~ "" | "-" ~~DigitCharacter ...],
       nowsize=getPopSize[scheme, preSize];
       1/(2 (nowsize - 1)) + (1 - 1/(nowsize - 1)) 1/preSize,
       StringMatchQ[scheme, "WF-E"~~ "" | "-" ~~DigitCharacter ...],
       1/(2 preSize),  
       StringMatchQ[scheme, "RM1-E"|"RM11-E"|"RM2-E"~~ "" | "-" ~~DigitCharacter ...],
       1/(2 (preSize-1)),  
       True,
       Switch[scheme,
            "Pairing", 0,
            "Selfing", 1,
            "Sibling", 1/2,
            "HalfDiallel",1/(preSize + 1),
            "FullDiallel",(2 (preSize - 1) - 1)/(preSize (preSize - 1) - 1) 1/2,
            _, Print["getCoalProb2::Wrong mating scheme ",scheme];
               Beep[];
               Abort[];
       ]
    ]
    ]
    


getCoalProb3[scheme_, preSize_] :=
    Module[ {p, temp,nowsize},
        Which[
            StringMatchQ[scheme, "Selfing" ~~ DigitCharacter .. ~~ "-" ~~ DigitCharacter ..],
               0,
            StringMatchQ[scheme, "WF-NE"|"RM1-NE"|"RM11-NE"~~ "" | "-" ~~DigitCharacter ...],
            1/preSize (1 - 1/preSize),
            StringMatchQ[scheme, "RM2-NE"~~ "" | "-" ~~DigitCharacter ...],
            nowsize=getPopSize[scheme, preSize];
            temp = 1-1/(nowsize-2)-(1-2/(nowsize-2))/preSize;
            temp*=(1 - 1/(nowsize - 1)) 1/preSize;
            temp+=1/(2 (nowsize - 1)) (1-1/preSize);
            temp,
            StringMatchQ[scheme, "WF-E"~~ ""|"-" ~~DigitCharacter ...],
            1/(2 preSize),  
            StringMatchQ[scheme, "RM1-E"|"RM11-E"|"RM2-E"~~ "" | "-" ~~DigitCharacter ...],
            1/(2 (preSize-1)), 
            True,
            Switch[scheme,
             "Pairing", 0,
             "Selfing", 0,
             "Sibling", 1/2,
             "HalfDiallel",
                  (*((n-1)-1)/(n(n-1)/2-1) =2/n+1*)
                 p = 1/(preSize + 1);
                 temp = (preSize - 3)/(preSize (preSize - 1)/2 - 2);
                 p (1 - temp/2),
             "FullDiallel",
                 p = (2 (preSize - 1) - 1)/(preSize (preSize - 1) - 1) 1/2;
                 temp = (2 (preSize - 1) - 2)/(preSize (preSize - 1) - 2);
                 p (1 - temp/2),     
             _, Print["getCoalProb3::Wrong mating scheme ",scheme];
                Beep[];
                Abort[];
             ]
        ]
    ]

(*
getfounderGamma12[isInbredFounder_] :=
    {1 - Boole[isInbredFounder],1}    
    
getfounderGamma123[isInbredFounder_,nFounder_] :=
    {1 - Boole[isInbredFounder],Boole[nFounder>= 3]}    
*)
    
nonIBDProb2[isRandomSelfing_, founderGamma12_,coalProb2_,popSize_] :=
    Module[ {a, b,selfprob,t},
         (*calculate non-IBD probability Alpha(12) and Beta(12) by recurrent equation*)
        a = b = Table[0, {Length[coalProb2]}];
        {a[[1]],b[[1]]} = founderGamma12;
        Do[
         If[ popSize[[t]] >= 2,
             b[[t]] = coalProb2[[t]]/2 a[[t - 1]] + (1 - coalProb2[[t]]) b[[t - 1]];
         ];
         a[[t]] = If[ isRandomSelfing[[t-1]],
                       (*"WF-*" or "Selfing"*)
                      selfprob = 1/popSize[[t]];
                      selfprob/2 a[[t - 1]] + (1 - selfprob) b[[t - 1]],
                      b[[t - 1]]
                  ], {t, 2, Length[a]}];
        (*Alpha(11)=1-Alpha(12)*)
        {a, b}
    ]

nonIBDProb3[isRandomSelfing_, founderGamma123_,coalProb2_, coalProb3_, popSize_] :=
    Module[ {a, b, t},
        (*calculate non-IBD probability a=Alpha(123) and b= Beta(123) by recurrent equation*)
        a = b = Table[0, {Length[coalProb2]}];
        {a[[1]],b[[1]]} = founderGamma123;
        Do[
         If[ popSize[[t]] >= 3,
             b[[t]] = 3/2 coalProb3[[t]] a[[t - 1]] + (1 - coalProb2[[t]]- 2 coalProb3[[t]]) b[[t - 1]];
         ];
         a[[t]] = If[ isRandomSelfing[[t-1]],
                      b[[t]],
                      coalProb2[[t]] a[[t - 1]] + (1 - 2 coalProb2[[t]]) b[[t - 1]]
                  ], {t, 2, Length[a]}];
        {a, b}
    ]
    
junc1232[isRandomSelfing_, founderJK1232_,coalProb2_, alpha123_] :=
    Module[ {a, b, t},
        (*calculate non-IBD probability a=J(1213) and b= K(1213) by recurrent equation*)
        a = b = Table[0, {Length[coalProb2]}];
        {a[[1]],b[[1]]} = founderJK1232;
        Do[
         b[[t]] = coalProb2[[t]]/2 a[[ t - 1]] + (1 - coalProb2[[t]]) (b[[t - 1]] + alpha123[[t - 1]]);
         a[[t]] = If[ isRandomSelfing[[t-1]],
                      b[[t]],
                      b[[t - 1]] + alpha123[[t - 1]]
                  ], {t, 2, Length[a]}];
        {a, b}
    ]

(*junc1222[isRandomSelfing_, founderJK1222_,coalProb2_, alpha12_, alpha121_] :=
    Module[ {a, b,t},
        (*calculate non-IBD probability a=J(1211) and b= K(1211) by recurrent equation*)
        a = b = Table[0, {Length[coalProb2]}];
        {a[[1]],b[[1]]}=founderJK1222;
        Do[
         b[[t]] = coalProb2[[t]]/2 (a[[t - 1]] + alpha12[[t - 1]]) 
                   +(1 - coalProb2[[t]]) (b[[t - 1]] + alpha121[[t - 1]]);
         a[[t]] = If[ isRandomSelfing[[t-1]],
                      b[[t]],
                      b[[t - 1]] + alpha121[[t - 1]]
                  ], {t, 2, Length[a]}];
        {a, b}
    ]
*)

junc1122[isRandomSelfing_, founderJK1122_,coalProb2_, mapR_] :=
    Module[ {a, b,t},
        a = b = Table[0, {Length[coalProb2]}];
        {a[[1]],b[[1]]} = founderJK1122;
        Do[
         b[[t]] = coalProb2[[t]]/2 a[[t - 1]] + coalProb2[[t]]/2 mapR[[t - 1]] + (1 - coalProb2[[t]]) b[[t - 1]];
         a[[t]] = If[ isRandomSelfing[[t-1]],
                      b[[t]],
                      b[[t - 1]]
                  ], {t, 2, Length[a]}];
        {a, b}
    ]

inbredProb[isRandomSelfing_, founderGamma12_,coalProb2_,popSize_] :=
    1-First[nonIBDProb2[isRandomSelfing, founderGamma12, coalProb2,popSize]]

mapExpansion[isRandomSelfing_,founderGamma12_, founderR_,coalProb2_,popSize_] :=
    Module[ {alpha12},
        alpha12 = First[nonIBDProb2[isRandomSelfing, founderGamma12,coalProb2,popSize]];
        founderR+Prepend[Most[Accumulate[alpha12]], 0]
    ]
      
          
sumJuncDensity[isRandomSelfing_, founderGamma12_, founderR_,founderJK1122_,coalProb2_,popSize_] :=
    Module[ {mapR,j1122},
        mapR = mapExpansion[isRandomSelfing, founderGamma12,founderR,coalProb2,popSize];
        j1122 = First[junc1122[isRandomSelfing,founderJK1122, coalProb2, mapR]];
        2 mapR - j1122
    ]   
       
    
juncDensity[isRandomSelfing_, founderGamma12_,founderGamma123_, founderR_,founderJK1122_,founderJK1232_, coalProb2_,coalProb3_, popSize_] :=
    Module[ {mapR,alpha12, alpha123,j1232,j1222,j1122},
        mapR = mapExpansion[isRandomSelfing, founderGamma12,founderR,coalProb2,popSize];
        alpha12 = First[nonIBDProb2[isRandomSelfing, founderGamma12, coalProb2,popSize]];
        alpha123 = First[nonIBDProb3[isRandomSelfing, founderGamma123,coalProb2, coalProb3, popSize]];
        j1232 = First[junc1232[isRandomSelfing, founderJK1232,coalProb2, alpha123]];
        j1122 = First[junc1122[isRandomSelfing,founderJK1122, coalProb2, mapR]];
        (*alpha121 = (alpha12-alpha123)/2;
        j1222 = First[junc1222[isRandomSelfing, founderJK1222,coalProb2, alpha12, alpha121]];*)
        j1222 = (mapR - j1232 - j1122)/2;
        {j1122,j1222,j1232}
    ]

origSummary[nFounder_, mateScheme_, founderGamma12_, founderGamma123_, founderR_, founderJK1122_, founderJK1232_] :=
    Module[ {isRandomSelfing, nGeneration, popSize, coalProb2, coalProb3, 
      finb, mapR, sumRho, juncRho, t},
        nGeneration = Length[mateScheme] + 1;
        isRandomSelfing = getSelfing[#] & /@ mateScheme;
        nGeneration = Length[mateScheme] + 1;
        popSize = FoldList[getPopSize[#2, #1] &, nFounder, mateScheme];
        coalProb2 = ConstantArray[0, nGeneration];
        coalProb3 = ConstantArray[0, nGeneration];
        Do[coalProb2[[t]] = getCoalProb2[mateScheme[[t - 1]], popSize[[t - 1]]];
           coalProb3[[t]] = getCoalProb3[mateScheme[[t - 1]], popSize[[t - 1]]], {t, 2, Length[coalProb2]}];
        finb = inbredProb[isRandomSelfing, founderGamma12, coalProb2, popSize];
        mapR = mapExpansion[isRandomSelfing, founderGamma12, founderR, coalProb2, popSize];
        sumRho = sumJuncDensity[isRandomSelfing, founderGamma12,founderR,founderJK1122, coalProb2, popSize];
        (*{j1122,j1222,j1232}*)
        juncRho = juncDensity[isRandomSelfing, founderGamma12, founderGamma123, 
        	founderR, founderJK1122, founderJK1232, coalProb2, coalProb3, popSize];
        {finb, mapR, sumRho, juncRho}
    ]

origSummary[nFounder_, mateScheme_] :=
    Module[ {isInbredFounder = True, founderGamma12, founderGamma123, 
      founderR, founderJK1122, founderJK1232},
        founderGamma12 = {1-Boole[isInbredFounder], 1};
        founderGamma123 = {1-Boole[isInbredFounder], Boole[nFounder >= 3]};
        founderR = 0;
        founderJK1122 = founderJK1232 = {0, 0};        
        origSummary[nFounder, mateScheme, founderGamma12, founderGamma123, 
         founderR, founderJK1122, founderJK1232]
    ]

magicOrigPrior[nFounder_, mateScheme_]:=
	Module[{temp,inbred,j1122,j1222,j1232,j1211,j1213},    
		temp = origSummary[nFounder, mateScheme];
        {inbred,j1122,j1222,j1232}=N[Prepend[temp[[4, All, -1]], temp[[1, -1]]]];
        (*{inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp} *)
        {j1211,j1213}={j1222,j1232};
        {inbred,j1122,j1211,j1213,j1222,j1232}    
	]
      
magicStationaryProb[nFgl_,inbredf_] :=
    Module[ {prob},
        prob = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        prob /. {{i_, i_} :> inbredf/nFgl, {_, _} :> (1 - inbredf)/(nFgl (nFgl - 1))}
    ]
    
(*juncrho={j1122,j1222,j1232}; f=inbredf; L=nFgl;*)
(*f r11 (L-1)=J1122;
(1-f) r01 =J1222;
(1-f) r00 (L-2) =J1232;
f r10 (L-1) =J1121 (=J1222);*) 
magicBasicRate[nFgl_, inbredf_, juncrho_] :=
    Module[ {r00, r01, r10, r11},
        (*{r11, r01, r00} = (juncrho/{inbredf, 1 - inbredf, 1 - inbredf})/{nFgl -1, 1, nFgl - 2};*)
        r11 = If[ inbredf==0,
                  0,
                  juncrho[[1]]/(inbredf (nFgl-1))
              ];
        r01 = juncrho[[2]]/(1-inbredf);
        r00 = If[ nFgl>2,
                  juncrho[[3]]/((1-inbredf)(nFgl-2)),
                  0
              ];
        r10 = If[ inbredf==0,
                  0,
                  r01 (1 - inbredf)/(inbredf (nFgl - 1))
              ];
        {r00, r01, r10, r11}
    ]
    
magicRateMatrix[nFgl_Integer?(#>=2&), baseR:{_?AtomQ, _?AtomQ,_?AtomQ, _?AtomQ}] :=
    Module[ {lab, qq, q00, q11, r00, r01, r10, r11},
        {r00, r01, r10, r11} = baseR;
        lab = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        qq = Outer[List, lab, lab, 1];
        q11 = -(nFgl - 1) ( r11 + 2 r10);
        q00 = -2 (nFgl - 2)  r00 - 2 r01;
        qq = qq /. {{{x1_, x1_}, {x1_, x1_}} -> q11, 
               {{x1_, x2_}, {x1_,x2_}} -> q00, 
               {{x1_, x1_}, {x2_,x2_}} -> r11, 
               {{x1_, x1_}, {x2_, x1_}} | {{x1_,x1_}, {x1_, x2_}} -> r10, 
               {{x1_, x2_}, {x2_, x2_}} | {{x1_, x2_}, {x1_, x1_}} -> r01, 
               {{x1_, x2_}, {x1_, x3_}} | {{x1_, x2_}, {x3_, x2_}} -> r00, {{_, _}, {_, _}} -> 0};
        qq
    ]

magicRateMatrix[nFgl_,inbredf_,juncrho_] :=
    magicRateMatrix[nFgl, magicBasicRate[nFgl,inbredf,juncrho]]

End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicOrigin`*"];

EndPackage[]

