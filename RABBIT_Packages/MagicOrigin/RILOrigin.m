(* Mathematica Package *)

BeginPackage["MagicOrigin`RILOrigin`"]

inbredProbRIL::usage = "inbredProbRIL  "

mapExpansionRIL::usage = "mapExpansionRIL  "

sumJuncDensityRIL::usage = "sumJuncDensityRIL  "

juncDensityRIL::usage = "juncDensityRIL  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

getgCross[nPower_, inbredScheme_] :=
    Which[
        inbredScheme == "Selfing", 
        nPower - 1, 
        inbredScheme == "Sibling", 
        Max[nPower - 2, 0], 
        True, 
         Print["getgCross: inbredScheme must be either Selfing or Sibling!"];
         Abort[]
     ]
    
(*inbredProbRIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, inibeta12,a12,t},
        gCross = getgCross[nPower, inbredScheme];
        inibeta12 = getinibeta12[nPower];
        a12 = Table[Evaluate[sibIBDProbXY2[inibeta12, t][[2]]],{t,0,gInbred}];
        Join[{1}, Table[0, {gCross}], 1 - a12]
    ]
*)             

inbredProbRIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, gls, a12, lam1 = (1 + Sqrt[5])/4, lam2 = (1 - Sqrt[5])/4},
        gCross = getgCross[nPower, inbredScheme];
        gls = Range[0, gInbred];
        Which[
             inbredScheme == "Selfing",
             a12 = (1/2)^gls,
             inbredScheme == "Sibling" && nPower == 1,
             a12 = (5 + Sqrt[5])/10 lam1^gls + (5 - Sqrt[5])/10  lam2^gls,
             inbredScheme == "Sibling" && nPower >= 2,
             a12 = (5 + 3 Sqrt[5])/10 lam1^gls + (5 - 3 Sqrt[5])/10  lam2^gls,
             True,
             Print["inbredProbRIL: inbredScheme must be either Selfing or Sibling!"];
             Abort[]
        ];
        Join[{1}, Table[0., {gCross}], 1 - a12]
    ]

mapExpansionRIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, gls, mp, lam1 = (1 + Sqrt[5])/4, lam2 = (1 - Sqrt[5])/4},
        gCross = getgCross[nPower, inbredScheme];
        gls = Range[0, gInbred];
        Which[
            inbredScheme == "Selfing",
             mp = gCross + 2 (1 - (1/2)^gls),
             inbredScheme == "Sibling" && nPower == 1,
             mp = 4 - (10 + 4 Sqrt[5])/5 lam1^gls - (10 - 4 Sqrt[5])/5  lam2^gls,
             inbredScheme == "Sibling" && nPower >= 2,
             mp = gCross + 6 - (15 + 7 Sqrt[5])/5 lam1^gls - (15 - 7 Sqrt[5])/5  lam2^gls,
             True,
             Print["mapExpansionRIL: inbredScheme must be either Selfing or Sibling!"];
             Abort[]
         ];
        Join[{0}, Range[0, gCross - 1], mp]
    ]

junc1232RIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, gls, j1232, lam1 = (1 + Sqrt[5])/4, lam2 = (1 - Sqrt[5])/4},
        gCross = getgCross[nPower, inbredScheme];
        gls = Range[0, gInbred];
        Which[
            inbredScheme == "Selfing",
            j1232 = gCross (1/2)^gls,
            inbredScheme == "Sibling" && nPower == 1,
            j1232 = Table[0, {gInbred + 1}],
            inbredScheme == "Sibling" && nPower >= 2,
            j1232 = -2 (1/2)^gls + (5 + 3 Sqrt[5])/10 (gCross + 2) lam1^gls + (5 - 3 Sqrt[5])/10 (gCross + 2)  lam2^gls,
            True,
            Print["junc1232RIL: inbredScheme must be either Selfing or Sibling!"];
            Abort[]
         ];
        Join[{0}, Range[0, gCross - 1], j1232]
    ]

junc1122RIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, gls, j1122, lam1 = (1 + Sqrt[5])/4, lam2 = (1 - Sqrt[5])/4},
        gCross = getgCross[nPower, inbredScheme];
        gls = Range[0, gInbred];
        Which[
             inbredScheme == "Selfing",
             j1122 = gCross (1 - (1/2)^gls) + 2 (1 - (gls + 1) (1/2)^gls),
             inbredScheme == "Sibling" && nPower == 1,
             j1122 = 4 - ((50 + 22 Sqrt[5])/25 + (3 + Sqrt[5])/5 gls) lam1^gls 
                        - ((50 - 22 Sqrt[5])/25 + (3 - Sqrt[5])/5 gls) lam2^gls,
             inbredScheme == "Sibling" && nPower >= 2,
             j1122 = gCross + 6 - ((75 + 31 Sqrt[5])/25 + (5 + 3 Sqrt[5])/10 gCross + (4 + 2 Sqrt[5])/5 gls) lam1^gls 
                                 - ((75 - 31 Sqrt[5])/25 + (5 - 3 Sqrt[5])/10 gCross + (4 - 2 Sqrt[5])/5 gls) lam2^gls,
             True,
             Print["junc1122RIL: inbredScheme must be either Selfing or Sibling!"];
             Abort[]
         ];
        Join[Table[0, {gCross + 1}], j1122]
    ]

junc1222RIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {gCross, gls, j1222, lam1 = (1 + Sqrt[5])/4, lam2 = (1 - Sqrt[5])/4},
        gCross = getgCross[nPower, inbredScheme];
        gls = Range[0, gInbred];
        Which[
             inbredScheme == "Selfing",
             j1222 = gls (1/2)^gls,
             inbredScheme == "Sibling" && nPower == 1,
             j1222 = (Sqrt[5]/25 + (3 + Sqrt[5])/10 gls) lam1^gls 
                     + (-Sqrt[5]/25 + (3 - Sqrt[5])/10 gls) lam2^gls,
             inbredScheme == "Sibling" && nPower >= 2,
             j1222 = (1/2)^gls - ((25 + 19 Sqrt[5])/50 - (2 + Sqrt[5])/5 gls) lam1^gls 
                                - ((25 - 19 Sqrt[5])/50 - (2 - Sqrt[5])/5 gls) lam2^gls,
             True,
             Print["junc1222RIL: inbredScheme must be either Selfing or Sibling!"];
             Abort[]
         ];
        Join[Table[0, {gCross + 1}], j1222]
    ]

sumJuncDensityRIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {R,j1122},
        R=mapExpansionRIL[nPower, gInbred, inbredScheme];
        j1122 = junc1122RIL[nPower, gInbred, inbredScheme];
        2 R - j1122
    ]    
    
juncDensityRIL[nPower_, gInbred_, inbredScheme_] :=
    Module[ {j1232,j1222,j1122},
        j1232 = junc1232RIL[nPower, gInbred, inbredScheme];
        j1222 = junc1222RIL[nPower, gInbred, inbredScheme];
        j1122 = junc1122RIL[nPower, gInbred, inbredScheme];
        {j1122,j1222,j1232}
    ]    
    
End[] (* End Private Context *)

EndPackage[]