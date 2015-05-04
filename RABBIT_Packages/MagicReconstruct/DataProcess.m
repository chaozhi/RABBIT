(* Mathematica Package *)

BeginPackage["MagicReconstruct`DataProcess`"]

Unprotect @@ Names["MagicReconstruct`DataProcess`*"];
ClearAll @@ Names["MagicReconstruct`DataProcess`*"];

SNPValidation::usage = "SNPValidation  "

transformMagicSNP::usage = "transformMagicSNP  "

(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

checkgeneticmap[magicSNP_] :=
    Module[ {ls, ls2,chr,res=True},
        ls = Transpose[magicSNP[[3 ;; 4, 2 ;;]]];
        ls = SplitBy[ls, First];
        chr = ls[[All, 1, 1]];
        If[ Length[chr] != Length[Union[chr]],
            Print["Markers on the same chromosome must be in neighbor columns"];
            res = False;
        ];
        ls2 = ls[[All, All, 2]];
        ls2 = OrderedQ[#] & /@ ls2;
        ls2 = Pick[ls[[All, 1, 1]], ls2, False];
        If[ ls2 =!= {},
            Print["The genetic distances must be increasing on chromosomes ", ls2];
            res=False
        ];
        res
    ]

checkfounderSNP[magicSNP_] :=
    Module[ {nfounder, ls},
        nfounder = magicSNP[[1, 2]];
        If[ ! (IntegerQ[nfounder] && nfounder >= 2),
            Print["the number of inbred founders must be an integer greater than or equal to 2!"]
        ];
        ls = Transpose[magicSNP[[5 ;; 4 + nfounder, 2 ;;]]];
        ls = ToString[#] & /@ Union[Flatten[ls]];
        ls = Complement[ls, {"1", "2", "N"}];
        If[ ls == {},
            True,
            Print["The alleles of founder haplotypes must be {1,2,N},where N denotes missing alelle.\nThe alleles ", ls, 
             " are not allowed! Make sure that the first founder haplotype is in row 5 and the number of founders is correct!"];
            False
        ]
    ]

checksampleSNP[magicSNP_] :=
    Module[ {nfounder, chrs, chr, posA, posX, geno, ls, pos, ls1, ls2, 
      res = True},
        nfounder = magicSNP[[1, 2]];
        chrs = Split[magicSNP[[3, 2 ;;]]][[All, 1]];
        chr = SplitBy[magicSNP[[3, 2 ;;]]][[All, 1]];
        If[ ! (Length[chr] == Length[Union[chr]]),
            Print["The SNP markers on the same chromsomes must be in order!"];
            res = False
        ];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        geno = SplitBy[Transpose[Join[magicSNP[[{3}, 2 ;;]], magicSNP[[5 + nfounder ;;, 2 ;;]]]], First];
        If[ posA=!={},
            ls = Transpose[Transpose[#] & /@ geno[[posA, All, 2 ;;]]];
            ls = Map[ToString, Union[Flatten[#]] & /@ ls, {2}];
            ls = Complement[#, {"NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"}] & /@ ls;
            pos = Flatten[Position[ls, _?(# =!= {} &), {1}, Heads -> False]];
            If[ pos =!= {},
                Print["The genotypes for sampled individuals on autosomals must be ", {"NN", "N1", "1N", "N2", "2N", "11", "12", "21", "22"}, ".\n", 
                 "The following sampled individuals have illegal genotypes: \n", 
                 Join[{{"Sampled individual", "Illegal genotypes"}}, 
                   Transpose[{magicSNP[[5 + nfounder ;;, 1]][[pos]], ls[[pos]]}]] //TableForm];
                res = False;
            ];
        ];
        If[ posX=!={},
            ls = Transpose[Transpose[#] & /@ geno[[posX, All, 2 ;;]]];
            ls = Map[ToString, Union[Flatten[#]] & /@ ls, {2}];
            ls1 = Complement[#, {"NN", "N1", "1N", "N2", "2N", "11", "12", "21",
                   "22"}] === {} & /@ ls;
            ls2 = Complement[#, {"N", "1", "2"}] === {} & /@ ls;
            pos = Flatten[Position[MapThread[Or, {ls1, ls2}], False]];
            If[ pos =!= {},
                Print["The genotypes on female XX chromsomes must be ", {"NN", 
                  "N1", "1N", "N2", "2N", "11", "12", "21", "22"}, ".\n", 
                 "The haplotypes on male X chromsomes must be ", {"N", "1", "2"}, 
                 ".\n", 
                 "The following sampled individuals have illegal genotypes on X chromosomes: \n", 
                 Join[{{"Sampled individual", "Illegal genotypes"}}, 
                   Transpose[{magicSNP[[5 + nfounder ;;, 1]][[pos]], ls[[pos]]}]] //TableForm];
                res = False;
            ];
        ];
        res
    ]

SNPValidation[magicSNP_] :=
    If[ ! (checkgeneticmap[magicSNP] && checkfounderSNP[magicSNP] && 
        checksampleSNP[magicSNP]),
        Abort[]
    ]



transformGeno[geno_] :=
    Module[ {strgeno = geno,res},
        strgeno[[3 ;;, 2 ;;]] = Map[ToString, strgeno[[3 ;;, 2 ;;]], {2}];
        res = SplitBy[Transpose[strgeno[[2 ;;, 2 ;;]]], First];
        res = Transpose[Transpose[#] & /@ res[[All, All, 2 ;;]]];
        res
    ]
    
jittersnp[snploc_] :=
    Module[ {factor = 5, gsnp, len, pos, j, delt},
        gsnp = SplitBy[snploc];
        len = Length[#] & /@ gsnp;
        pos = Flatten[Position[len - 1, _?Positive]];
        Do[
         Which[
          1 < j < Length[gsnp],
          delt = Differences[gsnp[[{j - 1, j, j + 1}, -1]]]/(len[[j]] factor);
          gsnp[[j]] +=Join[delt[[1]] Range[-Floor[len[[j]]/2], 0], 
            delt[[2]] Range[1, Ceiling[len[[j]]/2] - 1]],
          j == 1,
          delt = (gsnp[[j + 1, -1]] - gsnp[[j, -1]])/(len[[j]] factor);
          gsnp[[j]] += Range[0, len[[j]] - 1] delt,
          j == Length[gsnp],
          delt = (gsnp[[j, -1]] - gsnp[[j - 1, -1]])/(len[[j]] factor);
          gsnp[[j]] += Range[-len[[j]] + 1, 0] delt
          ], {j, pos}];
        gsnp = Flatten[gsnp];
        If[ Length[Union[gsnp]] < Length[gsnp],
            Print["Wrong in jittering overlapped SNPs!"];
            Abort[],
            gsnp
        ]
    ]
  
jittermap[snpMap_] :=
    Module[ {map = snpMap, ls, pos},
        ls = SplitBy[snpMap[[2 ;;, 2 ;;]], First];
        pos = Flatten[Position[Length[Union[#]] < Length[#] & /@ ls[[All, All, 2]],True]];
        If[ pos =!= {},
            Print["Warning: jittering the overlapped SNPs on chromosomes ", ls[[pos, 1, 1]],"!"];
            ls[[pos, All, 2]] = jittersnp[#] & /@ ls[[pos, All, 2]];
            map[[2 ;;, 2 ;;]] = Flatten[ls, 1];
        ];
        map
    ]    
    
transformMagicSNP[magicSNP_List] :=
    Module[ {snpMap, founderData, obsData,founderid,sampleid,deltd,founderHaplo,obsGeno,nFounder,chrs,posX,posA,genders},
        nFounder = magicSNP[[1,2]];
        snpMap = Transpose[magicSNP[[2 ;; 4]]];
        snpMap = jittermap[snpMap]; (*jitter the genetic map of magicSNP[[4]]*)
        founderData = Join[magicSNP[[2;;3]],magicSNP[[5 ;; 4 + nFounder]]];
        obsData = Join[magicSNP[[2;;3]],magicSNP[[5 + nFounder ;;]]];
        founderid = ToString[#]&/@founderData[[3 ;;, 1]];
        sampleid = ToString[#]&/@obsData[[3 ;;, 1]];
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = Differences[#] & /@ (SplitBy[snpMap[[2 ;;, 2 ;;]], First][[All, All, 2]]);
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = 0.01 deltd;
        founderHaplo = transformGeno[founderData];
        obsGeno = transformGeno[obsData];
        (*extract posA and posX*)
        chrs = Split[founderData[[2, 2 ;;]]][[All, 1]];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        (*extract gender from genotype on sex chromsomes*)
        If[ posX === {},
            genders = Table["Monoecious", {Length[obsData]-2}],
            genders = Union[Flatten[#]] & /@ obsGeno[[All, posX]];
            genders = (Complement[#, {"1", "2", "N"}] === {} & /@ genders);
            genders = genders /. {True -> "Male", False -> "Female"};
        ];
        (*return*)
        {deltd, founderHaplo, obsGeno,genders,posA,posX,nFounder,founderid,sampleid,Transpose[snpMap]}
    ]
             

End[] (* End Private Context *)

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicReconstruct`DataProcess`*"];

EndPackage[]