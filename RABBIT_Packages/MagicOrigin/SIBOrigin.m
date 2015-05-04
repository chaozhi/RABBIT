(* Mathematica Package *)

BeginPackage["MagicOrigin`SIBOrigin`"]

sibIBDProb2::usage = "sibIBDProb2  "

sibIBDProb3::usage = "sibIBDProb3  "

sibmapExpansion::usage = "sibmapExpansion  "

sibjunc1122::usage = "sibjunc1122  "

sibjunc1232::usage = "sibjunc1232  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

sibIBDProb2[inigamma12_, t_] :=
    Module[ {lam,cc,alpha12},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};      
        cc=#.inigamma12&/@{{(5-Sqrt[5])/10,2 Sqrt[5]/5},{(5+Sqrt[5])/10,-2 Sqrt[5]/5}};     
        alpha12[x_]:=cc.lam^x;   
        {alpha12[t],alpha12[t+1]}
    ]

sibIBDProb3[inialpha123_, t_] := inialpha123 (1/2)^t

sibmapExpansion[inigamma12_, iniR_, t_] :=
    Module[ {lam,cc,Rinf},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        cc=#.inigamma12&/@{{(5+Sqrt[5])/5,(10+6 Sqrt[5])/5},{(5-Sqrt[5])/5,(10-6 Sqrt[5])/5}};
        Rinf=iniR+{2,4}.inigamma12;
        Rinf-cc.lam^t
    ]

sibjunc1122[inigamma12_, iniR_, iniJK1122_, t_] :=
    Module[ {lam,Rinf,cc,cc2,j1122},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        Rinf=iniR+{2,4}.inigamma12;
        cc=#.iniJK1122&/@{{(5-Sqrt[5])/10,2 Sqrt[5]/5},{(5+Sqrt[5])/10,-2 Sqrt[5]/5}};
        cc-={(5+3 Sqrt[5])/10,(5-3 Sqrt[5])/10} iniR;
        cc-=(#.inigamma12&/@{{(25+13 Sqrt[5])/25,(50+18 Sqrt[5])/25},{(25-13 Sqrt[5])/25,(50-18 Sqrt[5])/25}});
        cc2 =-#.inigamma12&/@ {{2/5, (2+2 Sqrt[5])/5},{2/5, (2-2 Sqrt[5])/5}};
        j1122[x_]:=Rinf+cc.lam^x+x (cc2.lam^x);
        {j1122[t],j1122[t+1]}        
    ]

sibjunc1232[inialpha123_,iniJK1232_, t_] :=
    Module[ {lam,cc,j1232},
        lam = {(1 + Sqrt[5])/4, (1 - Sqrt[5])/4};
        cc=#.iniJK1232&/@{{(5-Sqrt[5])/10,2 Sqrt[5]/5},{(5+Sqrt[5])/10,-2 Sqrt[5]/5}}; 
        cc+={(5+3 Sqrt[5])/5,(5-3 Sqrt[5])/5} inialpha123;
        j1232[x_]:=cc.lam^x-2 inialpha123 (1/2)^x;
        {j1232[t],j1232[t+1]-sibIBDProb3[inialpha123,t]}
    ]    
    

End[] (* End Private Context *)

EndPackage[]