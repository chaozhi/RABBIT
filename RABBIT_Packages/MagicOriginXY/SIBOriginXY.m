(* Mathematica Package *)

BeginPackage["MagicOriginXY`SIBOriginXY`"]

sibIBDProbXY2::usage = "sibIBDProbXY2  "

sibIBDProbXY3::usage = "sibIBDProbXY3  "

sibmapExpansionXY::usage = "sibmapExpansionXY  "

sibmapExpansionXY2::usage = "sibmapExpansionXY2  "

sibjunc1122XY::usage = "sibjunc1122XY  "

sibjunc1232XYsym::usage = "sibjunc1232XYsym  "

sibjunc1232XYdiff::usage = "sibjunc1232XYdiff  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

sibIBDProbXY2[inibeta12_, t_] :=
    Module[ {lam, pp, cc},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        pp = {2 lam - 1, {1, 1}};
        cc = Inverse[pp].inibeta12;
        Simplify[pp.(lam^t cc)]
    ]

sibIBDProbXY3[inialpha123_, t_] :=
    Module[ {lam4 = 1/2},        
        inialpha123 (lam4)^t
    ]

getCoef8to11[inibeta12_, iniR_] :=
    Module[ {lam, pp, c1, c2, c8, c9, c10, c11},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        pp = {2 lam - 1, {1, 1}};
        {c1, c2} = Inverse[pp].inibeta12;
        c8 = {2/3, 1/3}.iniR + 2/3 {c1, c2}.(1/(1 - lam));
        {c10, c11} = -(2 lam {c1, c2})/((1 - lam) (2 lam + 1));
        c9 = First[iniR] - c8 - c10 - c11;
        Simplify[{c8, c9, c10, c11}]
    ]

sibmapExpansionXY[inibeta12_, iniR_, t_] :=
    Module[ {lam, c8, c9, c10, c11, Rm, Rp},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        {c8, c9, c10, c11} = getCoef8to11[inibeta12, iniR];
        Rm = c8 + c9 (-1/2)^t + {c10, c11}.lam^t;
        Rp = c8 - 2 c9 (-1/2)^t + {c10, c11}.lam^(t - 1);
        {Rm, Rp}
    ]

sibmapExpansionXY2[inibeta12_, iniR_, t_] :=
    Module[ {lam, c8, c9, c10, c11, cc, Rx, Rd},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        {c8, c9, c10, c11} = getCoef8to11[inibeta12, iniR];
        cc = Simplify[2/3 ({c10, c11} (2 lam + 1)/(2 lam))];
        Rx = c8 + cc.lam^t;
        cc = Simplify[{c10, c11} (1 - 1/lam)/2];
        Rd = 3/2 c9  (-1/2)^t + cc.lam^t;
        {Rx, Rd}
    ]

sibjunc1122XY[inibeta12_, iniR_, iniJ1122_, t_] :=
    Module[ {lam, c8, c9, c10, c11, c18to19, ww, pp, c15to16, j1122},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        {c8, c9, c10, c11} = getCoef8to11[inibeta12, iniR];
        c18to19 = {c10, c11}/(4 lam);
        ww = {2 lam, {0, 0}};
        pp = {2 lam - 1, {1, 1}};
        c15to16 = Simplify[Inverse[pp].(iniJ1122 - c8 - {1, -1/2} c9 - ww.c18to19)];
        j1122 = c8 + {1, -1/2} c9 (-1/2)^t + ww. (lam^t c18to19) + pp. (lam^t (c15to16 + c18to19 t));
        j1122
    ]

sibjunc1232XYsym[inialpha123_,iniJ1232sym_, t_] :=
    Module[ {lam4 = 1/2, lam, phi4, c24, ww, pp,c21to22, j1232sym},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        (*f12(x)=-1/4 for x=1/2*)
        phi4 = 1/2 (lam4)/(-1/4);        
        c24 = phi4 inialpha123;
        ww = {{-1/phi4 + (2 lam4 - 1)}, {1}};
        pp = {2 lam - 1, {1, 1}};
        c21to22 = Simplify[Inverse[pp].(iniJ1232sym - ww.{c24})];
        j1232sym = pp. (lam^t c21to22) + ww.{c24 lam4^t};
        j1232sym
    ]
    
sibjunc1232XYdiff[inialpha123_, iniJ1232diff_, t_] :=
    Module[ {j1232diff, lam4 = 1/2},
        j1232diff = (iniJ1232diff - inialpha123/2) (-1/2)^t +inialpha123/2 lam4^t;
        j1232diff
    ]    

End[] (* End Private Context *)

EndPackage[]