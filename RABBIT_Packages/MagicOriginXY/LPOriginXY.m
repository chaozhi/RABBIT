(* Mathematica Package *)

BeginPackage["MagicOriginXY`LPOriginXY`"]

lpnonIBDProbXY2::usage = "lpnonIBDProbXY2  "

lpnonIBDProbXY3::usage = "lpnonIBDProbXY3  "

lpmapExpansionXY::usage = "lpmapExpansionXY  "

lpjunc1122XY::usage = "lpjunc1122XY  "

lpjunc1232XYsym::usage = "lpjunc1232XYsym  "

lpjunc1232XYdiff::usage = "lpjunc1232XYdiff  "

lpmapExpansionXY2::usage = "lpmapExpansionXY2  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

lpnonIBDProbXY2[inibeta12_, coalProb_, t_] :=
    Module[ {lam, pp, cc},
        lam = {1 - coalProb/3, -1/2, 1/4};
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        cc = Inverse[pp].inibeta12;
        pp.(lam^t cc)
    ] 
       
getCoef4to7[inibeta123_, coalProb_] :=
    Module[ {qq, s = coalProb},
        qq={{1,-3/s, -1,-1/2}, {1,1,1,1},{1,3/s,0,-2},{1,6/s, -4,4}};   
        (Inverse[qq].inibeta123)
    ]  
      
lpnonIBDProbXY3[inibeta123_, coalProb_, t_] :=
    Module[ {s=coalProb,lam, qq, c4to7},
        lam = {1 - s, -1/2, 1/4,-1/8};
        qq={{1,-3/s, -1,-1/2}, {1,1,1,1},{1,3/s,0,-2},{1,6/s, -4,4}};        
        c4to7 = getCoef4to7[inibeta123, s];
        qq.(lam^t c4to7)
    ]    

getCoef8to12[inibeta12_, iniR_,coalProb_] :=
    Module[ {s=coalProb, lam, pp, c1, c2, c3, c8, c9, c10, c11, c12},    	
        lam = {1 - s/3, -1/2, 1/4};
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        {c1, c2, c3} = Inverse[pp].inibeta12;
        c10 = 2 lam[[1]]/(2 lam[[1]] + 1) c1;
        c11 = -2/3 c2;
        c12 = -4/9 c3;
        c8 = {2/3, 1/3}.iniR + 2/(3 (2  lam[[1]] + 1)) c1 + 4/9 c2 + 8/9 c3;
        c9 = iniR[[1]] - c8 - c12;
        {c8, c9, c10, c11, c12}/.s->0
    ]

lpmapExpansionXY[inibeta12_, iniR_, coalProb_, t_] :=
    Module[ {lam, c8, c9, c10, c11, c12, Rm, Rp},
        lam = {1 - coalProb/3, -1/2, 1/4};
        {c8, c9, c10, c11, c12} = getCoef8to12[inibeta12, iniR,coalProb];
        Rm = c8 + c10 (1 - lam[[1]]^t)/(1 - lam[[1]]) 
        	 + (c9 + c11 t) lam[[2]]^t + c12 lam[[3]]^t;
        Rp = c8 + c10 (1 - lam[[1]]^(t - 1))/(1 - lam[[1]]) 
        	 + (c9 +c11 (t - 1)) lam[[2]]^(t - 1) + c12 lam[[3]]^(t - 1);
        {Rm, Rp}
    ]

lpmapExpansionXY2[inibeta12_, iniR_, coalProb_, t_] :=
    Module[ {lam,pp,c1,c2,c3,cc, Rx, Rd},
        lam = {1 - coalProb/3, -1/2, 1/4};
        lam = {1 - coalProb/3, -1/2, 1/4};
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        {c1,c2,c3} = Inverse[pp].inibeta12;        
        Rx = {2/3, 1/3}.iniR + 2/3 c1 (1 - lam[[1]]^t)/(1 - lam[[1]]) 
        	 + 4/9 c2 (1 - lam[[2]]^t) + 8/9 c3 (1 - lam[[3]]^t);
        cc={-1/2, -1/2}.iniR- c1/(2 lam[[1]]+1)-2/3 lam[[3]];      
        Rd = c1/(2 lam[[1]]+1) lam[[1]]^t + (cc -c2 t) lam[[2]]^t+ 2/3 c3 lam[[3]]^t;        
        {Rx, Rd}
    ]
    
lpjunc1122XY[inibeta12_, iniR_, iniJ1122_, coalProb_, t_] :=
    Module[ {s = coalProb, lam, pp, c8, c9, c10, c11, c12, c15to17, j1122},
        lam = {1 - s/3, -1/2, 1/4};
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        {c8, c9, c10, c11, c12} = getCoef8to12[inibeta12, iniR,coalProb];
        c15to17 = Inverse[pp].(iniJ1122 - c8);
        j1122 = c8 + c10 (1 - lam[[1]]^t)/(1 - lam[[1]]) 
        		- c10 t lam[[1]]^t (1 + {2/9, 8/9, -4/9} s) + pp. (lam^t c15to17);
        j1122
    ]

getCoef21to23[inibeta123_, iniJ1232_, coalProb_] :=
    Module[ {s = coalProb, pp, c4to7, ww},
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        c4to7 = getCoef4to7[inibeta123, s];
        ww = {{1/3, -2/(3 s), -7/9, -4/9}, {0, 0, 0, -4/9}, {-1/3, -4/(3 s), -20/9, 32/9}};
        Inverse[pp].(iniJ1232 - ww.c4to7)
    ]

lpjunc1232XYsym[inibeta123_, iniJ1232_, coalProb_, t_] :=
    Module[ {c4, c5, c6, c7, c21, c22, c23, s = coalProb, lam, pp, ww, 
      ww2, j1232},
        {c4, c5, c6, c7} = getCoef4to7[inibeta123, s];
        {c21, c22, c23} = getCoef21to23[inibeta123, iniJ1232, s];
        lam = {1 - s/3, -1/2, 1/4, 1 - s, -1/2, 1/4, -1/8};
        pp = {{1, -2, -1/2}, {1, 1, 1}, {1, 4, -2}};
        ww = {{1/3, -2/(3 s), -7/9, -4/9}, {0, 0, 0, -4/9}, {-1/3, -4/(3 s), -20/9, 32/9}};
        ww2 = {{-4/( 3 s), -2/9}, {2/(3 s), 4/9}, {8/(3 s), -8/9}};
        j1232 = (lam[[1]]^t - lam[[4]]^t) c4/s + pp.(lam[[;; 3]]^t {c21, c22, c23}) 
        		+ ww.(lam[[{1, 2, 3, 7}]]^t {c4, c5, c6, c7}) + ww2.(t lam[[2 ;; 3]]^t {c5, c6});
        j1232
    ]

lpjunc1232XYdiff[inibeta123_, iniJ1232diff_, coalProb_, t_] :=
    Module[ {c4, c5, c6, c7, s = coalProb, lam, a28, j1232diff},
        lam = {1 - s, -1/2, 1/4, -1/8};
        {c4, c5, c6, c7} = getCoef4to7[inibeta123, s];
        a28 = iniJ1232diff - 1/3 c4 - 2/3 c6 - 4/3 c7;
        j1232diff = {1/3 c4 , a28, 2/3 c6, 4/3 c7}.lam^t;
        j1232diff
    ]

End[] (* End Private Context *)

EndPackage[]