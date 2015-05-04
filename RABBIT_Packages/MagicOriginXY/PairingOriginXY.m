(* Mathematica Package *)

BeginPackage["MagicOriginXY`PairingOriginXY`"]

pairIBDProbXY2::usage = "pairIBDProbXY2  "

pairIBDProbXY3::usage = "pairIBDProbXY3  "

pairmapExpansionXY::usage = "pairmapExpansionXY  "

pairmapExpansionXY2::usage = "pairmapExpansionXY2  "

pairjunc1122XY::usage = "pairjunc1122XY  "

pairjunc1232XYsym::usage = "pairjunc1232XYsym  "

pairjunc1232XYdiff::usage = "pairjunc1232XYdiff  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

(*2^npower (npower>=2) fully inbred founders*)

pairIBDProbXY2[]:={1,1,1}

pairIBDProbXY3[]:={1,1,1,1}

pairmapExpansionXY[t_] :=
    Module[ {Rm, Rp},
        Rm = 2 (1 - (-1/2)^t + 3 t)/9;
        Rp = 4 (-1 + (-1/2)^t + 3/2 t)/9;
        {Rm, Rp}
    ]
        
pairmapExpansionXY2[t_] :=
    Module[ {Rm, Rp},
        Rm = 2 (1 - (-1/2)^t + 3 t)/9;
        Rp = 4 (-1 + (-1/2)^t + 3/2 t)/9;
        {{2/3,1/3}.{Rm, Rp},(Rm-Rp)/2}
    ]

pairjunc1122XY[]:={0,0,0}

pairjunc1232XYsym[t_]:=Module[ {Rm, Rp},
	{Rm, Rp}=pairmapExpansionXY[t];
	{Rm,(Rm+Rp)/2,Rp}	
]

pairjunc1232XYdiff[t_]:=Module[ {Rm, Rp},
	{Rm, Rp}=pairmapExpansionXY[t];
	(Rm-Rp)/2
]

End[] (* End Private Context *)
            
EndPackage[]