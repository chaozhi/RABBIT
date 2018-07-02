(* Mathematica Package *)

(* Created by the Wolfram Workbench Oct 11, 2016 *)

BeginPackage["MagicDataPreprocess`",{"MagicDefinition`"}]
(* Exported symbols added here with SymbolName::usage *) 
(* Mathematica Package *)

Unprotect @@ Names["MagicDataPreprocess`*"];
ClearAll @@ Names["MagicDataPreprocess`*"];

SNPValidation::usage = "SNPValidation  "

transformMagicSNP::usage = "transformMagicSNP  "

jitterGeneticMap::usage = "jitterGeneticMap  "

jitterMagicSNP::usage = "jitterMagicSNP  "

snpCluster::usage = "snpCluster  "

getClusterhaploMap::usage = "getClusterhaploMap  "

getDerivedGenotype::usage = "getDerivedGenotype  "

checkOptionValue::usage = "checkOptionValue  "

getsubMagicSNP::usage = "getsubMagicSNP  "

(*GBS data*)

calMagicReadDepth::usage = "calMagicReadDepth  "

randomMagicRead::usage = "randomMagicRead  "

rawGenotypeCall::usage = "rawGenotypeCall  "

magicRawGenotypeCall::usage = "magicRawGenotypeCall  "

magicInbredGenotypeCall::usage = "magicInbredGenotypeCall  "

magicSNPtoVCF::usage = "magicSNPtoVCF  "

vcftoMagicSNP::usage = "vcftoMagicSNP  "

magicSNPtoLBImpute::usage = "magicSNPtoLBImpute  "

magicSNPtompMap::usage = "magicSNPtompMap  "

mpMaptoMagicSNP::usage = "mpMaptoMagicSNP  "

(*magicSNP and R/HAPPY*)

happytoMagicSNP::usage = "happytoMagicSNP  "

getJoinMapCPInput::usage = "getJoinMapCPInput  "

calMagicMissingFraction::usage = "calMagicMissingFraction  "

recombinationRate::usage = "recombinationRate is an option"

tranformTarget::usage = "tranformTarget is an option"

isPhased::usage = "isPhased is an option"

getvcfsampleName::usage = "getvcfsampleName  "

genoofVCF::usage = "genoofVCF  "

genotoMagicSNP::usage = "genotoMagicSNP  "


(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *) 

getsubMagicSNP[magicsnp_, chrsubset_,snpsubset_] :=
    Module[ {ls},
        ls = Tally[magicsnp[[3, 2 ;;]]][[All, 2]];
        ls = Accumulate[Prepend[ls, 0]];
        ls = Transpose[{Most[ls] + 1, Rest[ls]}];
        ls = Range @@ # & /@ ls[[chrsubset]];
        ls = Flatten[#[[snpsubset]]&/@ls];
        Join[magicsnp[[{1}]], magicsnp[[2 ;;, Join[{1}, 1 + ls]]]]
    ]

checkgeneticmap[magicSNP_] :=
    Module[ {ls, ls2,chr,res = True},
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
            res = False
        ];
        res
    ]

checkfounderSNP[magicSNP_,isfounderinbred_] :=
    Module[ {nfounder,chrs,colX,colA,ls,set,res = True},
        nfounder = magicSNP[[1, 2]];
        If[ ! (IntegerQ[nfounder] && nfounder >= 1),
            Print["the number of inbred founders must be an integer greater than or equal to 1!"]
        ];
        chrs = Split[magicSNP[[3, 2 ;;]]][[All, 1]];
        If[ ! (Length[chrs] == Length[Union[chrs]]),
            Print["The SNP markers on the same chromsomes must be in order!"];
            res = False
        ];
        colX = Flatten[Position[magicSNP[[3, 2 ;;]], "X" | "x"]] + 1;
        colA = Complement[Range[2, Length[magicSNP[[3]]]], colX];
        ls = Map[ToString, Union[Flatten[magicSNP[[5 ;; nfounder + 4, #]]]] & /@ {colA, colX}, {2}];
        set = {If[ isfounderinbred,
                   {"1", "2", "N"},
                   {"11", "12", "1N", "21", "22", "2N", "N1", "N2", "NN"}
               ], 
        If[ isfounderinbred,
            {"1", "2", "N"},
            {"1", "2", "N", "11", "12", "1N", "21", "22", "2N", "N1", "N2", "NN"}
        ]};
        MapThread[
        If[ !SubsetQ[#1,#2],
            Print["The possible alleles of founder genotypes must be ",#[[1]],"; the alleles ", Complement[#2,#1], " are not allowed!"];
            res = False
        ]&,{set,ls}];
        res
    ]

checkoffspringSNP[magicSNP_] :=
    Module[ {nfounder, chrs, posA, posX, geno, ls, pos, ls1, ls2, 
      res = True},
        nfounder = magicSNP[[1, 2]];
        chrs = Split[magicSNP[[3, 2 ;;]]][[All, 1]];
        If[ ! (Length[chrs] == Length[Union[chrs]]),
            Print["The SNP markers on the same chromsomes must be in order!"];
            res = False
        ];
        posX = Flatten[Position[chrs, "X"|"x"]];
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

checkAllelicDepth[allelicdepth_] :=
    Module[ {geno, count, res},
        geno = StringSplit[#, "|"] & /@ allelicdepth;
        count = Union[Flatten[geno]];
        count = ToExpression[count];
        res = VectorQ[count, (NonNegative[#] && IntegerQ[#] &)];
        If[ ! res,
            Print["checkAllelicDepth: allelicDepth must be in form of x|y, where x and y are the read counts for alleles 1 and 2 repsecively."];
        ];
        res
    ]
  
SNPValidation[magicSNP_,isfounderinbred_,isfounderdepth_,isoffspringdepth_] :=
    Module[ {nfounder,bool1,bool2,bool3},
        nfounder = magicSNP[[1, 2]];
        bool1 = checkgeneticmap[magicSNP];
        bool2 = If[ isfounderdepth,
                    checkAllelicDepth[magicSNP[[5;;4+nfounder, 2 ;;]]],
                    checkfounderSNP[magicSNP,isfounderinbred]
                ];
        bool3 = If[ isoffspringdepth,
                    checkAllelicDepth[magicSNP[[5 + nfounder ;;, 2 ;;]]],
                    checkoffspringSNP[magicSNP]
                ];
        If[ (And @@ {bool1,bool2,bool3}),
            True,
            Print[{bool1,bool2,bool3}];
            Abort[]
        ]
    ]

transformGeno[geno_] :=
    Module[ {strgeno = geno,res},
        strgeno[[3 ;;, 2 ;;]] = Map[ToString, strgeno[[3 ;;, 2 ;;]], {2}];
        res = SplitBy[Transpose[strgeno[[2 ;;, 2 ;;]]], First];
        res = Transpose[Transpose[#] & /@ res[[All, All, 2 ;;]]];
        res
    ]
    
jittersnp[snploc_] :=
    Module[ {factor = 10., gsnp, len, pos, j, delt,maxdelt = 1},
        (*delt, maxdelt are in unit of cM*)
        gsnp = Split[N[snploc]];
        len = Length[#] & /@ gsnp;
        pos = Flatten[Position[len - 1, _?Positive]];
        Do[
         Which[
          1 < j < Length[gsnp],
          delt = Differences[gsnp[[{j - 1, j, j + 1}, -1]]];
          delt = Min[#,maxdelt]&/@delt;
          delt = delt/(len[[j]] factor);
          gsnp[[j]] +=Join[delt[[1]] Range[-Floor[len[[j]]/2], 0], 
            delt[[2]] Range[1, Ceiling[len[[j]]/2] - 1]],
          j == 1,
          delt = (gsnp[[j + 1, -1]] - gsnp[[j, -1]]);
          delt = Min[delt,maxdelt];
          delt = delt/(len[[j]] factor);
          gsnp[[j]] += Range[0, len[[j]] - 1] delt,
          j == Length[gsnp],
          delt = (gsnp[[j, -1]] - gsnp[[j - 1, -1]]);
          delt = Min[delt,maxdelt];
          delt = delt/(len[[j]] factor);
          gsnp[[j]] += Range[-len[[j]] + 1, 0] delt
          ], {j, pos}];
        gsnp = Flatten[gsnp];
        If[ Length[Union[gsnp]] < Length[gsnp],
            Print["jittersnp: there exist remaining overlapped SNPs! gsnp = ",gsnp];
            Put[snploc,gsnp,"temp.txt"];
            Abort[];
            (*Abort[]*)
        ];
        gsnp
    ]
  
jitterMagicSNP[magicSNP_List,isprint_:True] :=
    Module[ {map},
        map = Transpose[magicSNP[[2 ;; 4]]];
        Join[magicSNP[[{1}]],Transpose[jitterGeneticMap[map,isprint]],magicSNP[[5;;]]]
    ]  
    
jitterGeneticMap[snpMap_,isprint_:True] :=
    Module[ {map = snpMap, ls, pos},
        ls = SplitBy[snpMap[[2 ;;, 2 ;;]], First];
        pos = Flatten[Position[Length[Union[#]] < Length[#] & /@ ls[[All, All, 2]],True]];
        If[ pos === {},
            If[ isprint,
                Print["No overlapped SNPs!"]
            ],
            If[ isprint,
                Print["Jittering the overlapped SNPs on chromosomes ", ls[[pos, 1, 1]],"!"];
            ];
            ls[[pos, All, 2]] = jittersnp[#] & /@ ls[[pos, All, 2]];
            map[[2 ;;, 2 ;;]] = Flatten[ls, 1];
        ];
        map
    ]    
    
transformMagicSNP[magicSNP_List,isfounderinbred_,isfounderdepth_,isoffspringdepth_] :=
    Module[ {snpMap, founderData, obsData,founderid,sampleid,deltd,posfoundermale,ls,
        founderHaplo,obsGeno, nFounder,chrs,posX,posA,foundergender,offspringgender,haploid,haploMap},
        nFounder = magicSNP[[1,2]];
        foundergender = Table["NotApplicable",{nFounder}];
        snpMap = Transpose[magicSNP[[2 ;; 4]]];
        snpMap = jitterGeneticMap[snpMap,False]; (*jitter the genetic map of magicSNP[[4]]*)
        founderData = Join[magicSNP[[2;;3]],magicSNP[[5 ;; 4 + nFounder]]];
        obsData = Join[magicSNP[[2;;3]],magicSNP[[5 + nFounder ;;]]];
        founderid = ToString[#]&/@founderData[[3 ;;, 1]];
        sampleid = ToString[#]&/@obsData[[3 ;;, 1]];
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = Differences[#] & /@ (SplitBy[snpMap[[2 ;;, 2 ;;]], First][[All, All, 2]]);
        (*change the genetic distance from centiMorgan into Morgan*)
        deltd = 0.01 deltd;
        obsGeno = transformGeno[obsData];
        chrs = Split[founderData[[2, 2 ;;]]][[All, 1]];
        posX = Flatten[Position[chrs, "X"]];
        posA = Complement[Range[Length[chrs]], posX];
        If[ isoffspringdepth,
            obsGeno = ToExpression[Map[StringSplit[#, "|"] &, obsGeno, {2}]];
            If[ posX==={},
                offspringgender = Table["NotApplicable", {Length[sampleid]}],
                ls = StringSplit[sampleid, "_"];
                If[ Min[Length[#] & /@ ls] < 2,
                    Print["transformMagicSNP: offspring gender must be encoded in their id in form of offspringid_Male or offspringid_Female."];
                    Abort[]
                ];
                offspringgender = ls[[All, -1]] /. {"male" -> "Male", "female" -> "Female"};
                If[ !SubsetQ[{"female", "male"}, Union[ToLowerCase[offspringgender]]],
                    Print["transformMagicSNP: offspring gender must be encoded in their id in form of offspringid_Male or offspringid_Female."];
                    Abort[]
                ];
            ];
            sampleid = StringDelete[sampleid, "_Female" | "_Male" | "_female" | "_male"],
            offspringgender = getgenders[obsGeno[[All, posX]], True]
        ];
        founderHaplo = transformGeno[founderData];
        If[ isfounderdepth,
            founderHaplo = ToExpression[Map[StringSplit[#, "|"] &, founderHaplo, {2}]];
        ];
        If[ (!isfounderinbred)&&(isfounderdepth),
            If[ posX=!={},
                foundergender = StringSplit[founderid, "_"][[All, -1]];
            ];
            founderid = StringDelete[founderid, "_Female" | "_Male" | "_female" | "_male"];
        ];
        If[ (!isfounderinbred)&&(!isfounderdepth),
            foundergender = getgenders[founderHaplo[[All, posX]], True];
            posfoundermale = Flatten[Position[getgenders[founderHaplo[[All, posX]], True], "Male"]];
            founderHaplo[[posfoundermale, posX]] = founderHaplo[[posfoundermale, posX]] /. {"1" -> "1N", "2" -> "2N","N" -> "NN"};
            founderHaplo = Transpose[#] & /@ Transpose[founderHaplo];
            founderHaplo = Map[Flatten[Characters[#]] &, founderHaplo, {2}];
            founderHaplo = Transpose[Transpose[#] & /@ founderHaplo];
        ];
        snpMap = Transpose[snpMap];
        haploid = If[ isfounderinbred,
                      founderid,
                      Flatten[Outer[StringJoin, founderid, {"_Maternal", "_Paternal"}]]
                  ];
        haploMap = Join[snpMap, Join[Transpose[{haploid}], Flatten[#] & /@ founderHaplo, 2]];
        (*return*)
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,founderid,sampleid}
    ]

checkOptionValue[magicSNP_, isfounderdepth_, isoffspringdepth_,isfounderinbred_,isprint_] :=
    Module[ {isfounderGBS, isoffspringGBS, nfounder,fgenoset},
        {isfounderGBS, isoffspringGBS} = {isfounderdepth, isoffspringdepth};
        nfounder = magicSNP[[1, 2]];
        fgenoset = ToString[#]&/@Union[Flatten[Flatten[magicSNP[[5 ;; 4 + nfounder, 2 ;;]]]]];
        If[ isfounderinbred=!=Automatic,
            If[ (Total[StringLength[#] & /@ fgenoset] == Length[fgenoset]) && (!isfounderinbred),
                Print["checkOptionValue: the option value for isFounderInbred is not compatible the founder marker data!"];
                Abort[];
            ];
        ];
        If[ isfounderdepth == Automatic,
            isfounderGBS = And @@ StringContainsQ[fgenoset, "|"];
            If[ isprint,
                If[ isfounderGBS,
                    PrintTemporary["checkOptionValue: input marker data of founders are detected to be allelic depths."],
                    PrintTemporary["checkOptionValue: input marker data of founders are detected to be called genotypes."];
                ];
            ];
        ];
        If[ isoffspringdepth == Automatic,
            isoffspringGBS = And @@ StringContainsQ[ToString[#]&/@Union[Flatten[Flatten[magicSNP[[5 + nfounder ;;, 2 ;;]]]]], "|"];
            If[ isprint,
                If[ isoffspringGBS,
                    PrintTemporary["checkOptionValue: input marker data of offspring are detected to be allelic depths."],
                    PrintTemporary["checkOptionValue: input marker data of offspring are detected to be called genotypes."];
                ];
            ];
        ];
        {isfounderGBS, isoffspringGBS}
    ]
  
getgenders[genos_,isdiploid_] :=
    Module[ {genders},
        If[ (Dimensions[genos][[2]] == 0) || (!isdiploid),
            genders = Table["NotApplicable", {Length[genos]}],
            genders = Union[Flatten[#]] & /@ genos;
            genders = (Complement[#, {"1", "2", "N"}] === {} & /@ genders);
            genders /. {True -> "Male", False -> "Female"}
        ]
    ]
(*clustering nearby 1D gridpoints  so that the length of each \
cluster <threshhold (in unit consistent with gridpoint) except \
singleton cluster, return the a list of cluster, 
each clustering represented by the indices of gridpoints*)
snpCluster[points_?OrderedQ, threshhold_] :=
    Module[ {dis, res, cluster, len,index,weight,i},
        dis = Transpose[{Range[Length[points]], 
           Join[Differences[points], {0}]}];
        res = cluster = {};
        len = 0;
        res = Reap[
           Do[
            len += dis[[i, 2]];
            If[ len < threshhold||i==1,
                AppendTo[cluster, dis[[i]]],
                Sow[cluster];
                len = dis[[i, 2]];
                cluster = {dis[[i]]};
            ], {i, Length[dis]}];
           Sow[cluster];
           ][[2, 1]];
        index = res[[All, All, 1]];
        weight = Normalize[#, Total] & /@ res[[All, All, 2]];
        {index,weight}
    ]            
    
getClusterhaploMap[haploMap_,outputsmooth_] :=
    Module[ {temp,weights,clusters,clusterhaploMap},
        temp = SplitBy[Transpose[haploMap[[All, 2 ;;]]], #[[2]] &];
        {clusters,weights} = Transpose[Map[snpCluster[#, outputsmooth] &, temp[[All, All, 3]]]];
        clusterhaploMap = Transpose[Flatten[MapThread[#1[[#2]] &, {temp, clusters[[All, All, 1]]}],1]];
        clusterhaploMap = Join[haploMap[[All,{1}]],clusterhaploMap,2];
        {clusters,weights,clusterhaploMap}
    ]    

getDerivedGenotype[founderHaplo_,isdepModel_,isphased_] :=
    Module[ {nFgl, nchr, genotypes, derivedRule,i,j},
        {nFgl, nchr} = Take[Dimensions[founderHaplo], 2];
        If[ isdepModel,
            genotypes = Transpose[Table[Range[nFgl], {2}]],
            If[ isphased,
                (*origGenotype[nFounder][[1, 2]]*)
                genotypes = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1],
                (*origGenotype[nFounder][[1, 1]]*)
                genotypes = Flatten[Table[{i, j}, {i, nFgl}, {j, i}], 1];
            ];
        ];
        derivedRule = {{"N", "N"} -> 1, {"N", "1"} -> 2, {"1", "N"} -> 
           3, {"N", "2"} -> 4, {"2", "N"} -> 5, {"1", "1"} -> 
           6, {"1", "2"} -> 7, {"2", "1"} -> 8, {"2", "2"} -> 9};
        Table[Replace[Transpose[Partition[founderHaplo[[Flatten[genotypes], i]], 2], {2, 3, 1}], derivedRule, {2}], {i, nchr}]
    ]     

Options[vcftoMagicSNP] = {
    recombinationRate -> 1,
    isFounderInbred ->True,
    outputFileID -> ""
    }
    
vcftoMagicSNP[vcffile_String?FileExistsQ,genoformat:"GT"|"AD"|"GP", founders_,OptionsPattern[]] :=
    Module[ {recomrate,isfounderinbred,outputid,outputfile,res,depth,founders2=founders},
        {recomrate,isfounderinbred,outputid} = OptionValue[{recombinationRate,isFounderInbred,outputFileID}];
        If[ outputid=!="",
            outputfile = outputid<>"_"<>genoformat <> ".csv",
            outputfile = "vcftoMagicSNP_"<>genoformat <> ".csv";
        ];
        depth = genoofVCF[vcffile,genoformat];
        If[Head[founders2]===Integer&&founders2>=1,
        	founders2=depth[[4 ;; 3 + founders2, 1]]
        ];
        res = genotoMagicSNP[depth, recomrate,founders2,isfounderinbred&&genoformat==="GT",outputfile];
        res
    ]
    
vcftoMagicSNP[foundervcf_String?FileExistsQ,offspringvcf_String?FileExistsQ, genoformat : "GT" | "AD"|"GP", OptionsPattern[]] :=
    Module[ {outputfile, isfounderinbred,recomrate, outputid, depths, nfounder, magicsnp},
        {recomrate,isfounderinbred,outputid} = OptionValue[{recombinationRate,isFounderInbred,outputFileID}];
        If[ outputid=!="",
            outputfile = outputid<>"_"<>genoformat <> ".csv",
            outputfile = "vcftoMagicSNP_"<>genoformat <> ".csv";
        ];
        depths = genoofVCF[#, genoformat] & /@ {foundervcf, offspringvcf};
        nfounder = Length[First[depths]] - 3;
        magicsnp = Join[{{"nfounder", nfounder}}, depths[[1]], depths[[2, 4 ;;]]];
        magicsnp[[4, 2 ;;]] *= recomrate 10^(-6.);
        magicsnp[[4, 1]] = "Pos(cM)";
        If[ isfounderinbred && genoformat === "GT",
            magicsnp[[5 ;; nfounder + 4, 2 ;;]] = 
             magicsnp[[5 ;; nfounder + 4, 2 ;;]] /. {11|"11"| "1N" | "N1" -> 1,22|"22"| "2N" | "N2" -> 2, 12 |21|"12"|"21"| "NN" -> "N"}
        ];
        csvExport[outputfile, magicsnp]
    ]    

Options[magicSNPtoVCF] = {
    recombinationRate -> 1,
    tranformTarget ->"All",
    isPhased ->False,
    outputFileID -> ""
    }

magicSNPtoVCF[inputmagicSNP_?(ListQ[#] || StringQ[#] &),OptionsPattern[]] :=
    Module[ {magicsnp = inputmagicSNP, outputfile,target,isphased,recomrate,outputid,nfounder, isfounderGBS, isoffspringGBS,noff,nsnp, data, rule, vcfheader,genoformat},
        {recomrate,target,isphased,outputid} = OptionValue[{recombinationRate,tranformTarget,isPhased,outputFileID}];
        If[ outputid=!="",
            outputfile = outputid <> ".vcf",
            outputfile = "magicSNPtoVCF_Output.vcf";
        ];
        If[ StringQ[magicsnp],
            If[ !FileExistsQ[magicsnp],
                Print["File ", magicsnp," does not exist!"];
                Return[$Failed]
            ];
            magicsnp = Import[magicsnp,"CSV"];
        ];
        nfounder = magicsnp[[1,2]];
        noff = Length[magicsnp[[nfounder + 5 ;;]]];
        nsnp = Length[magicsnp[[2,2;;]]];
        {isfounderGBS, isoffspringGBS} = checkOptionValue[magicsnp, Automatic, Automatic, Automatic, False];
        genoformat = If[ #,
                         "AD",
                         "GT"
                     ] & /@ {isfounderGBS, isoffspringGBS};
        Switch[ToLowerCase[target],
            "all",
            data = Rest[magicsnp];
            If[ Equal @@ genoformat,
                genoformat = First[genoformat],
                Print["magicSNPtoVCF: inconsistent genotypic format for founders and offspring: ",genoformat,"!"];
                Abort[];
            ],
            "founder"|"founders",
            data = magicsnp[[2;;nfounder+4]];
            genoformat = First[genoformat],
            "offspring",
            data = Join[magicsnp[[2;;4]],magicsnp[[nfounder+5;;]]];
            genoformat = Last[genoformat],
            _,
            Print["magicSNPtoVCF: wrong target = ",target,"!"];
            Abort[]
        ];
        If[ genoformat==="GT",
            If[ isphased,
                rule = {"1" | "11" -> "0|0", "2" | "22" -> "1|1", "12" -> "0|1", 
                       "21" -> "1|0", "N" | "NN" -> ".|.", "1N" -> "0|.", "2N" -> "1|.","N1" -> ".|0", "N2" -> ".|1"},
                rule = {"1" | "11" -> "0/0", "2" | "22" -> "1/1", "12" -> "0/1", 
                       "21" -> "1/0", "N" | "NN" -> "./.", "1N" -> "0/.", "2N" -> "1/.","N1" -> "./0", "N2" -> "./1"};
            ];
            (*the first three rows of data are genetic map*)
            data[[4;;, 2 ;;]] = Map[ToString, data[[4;;, 2 ;;]], {2}] /. rule,
            data[[4;;, 2 ;;]] = Map[StringReplace[#, "|" -> ","] &, data[[4;;, 2 ;;]], {2}]
        ];
        data = Transpose[data];
        data = Join[data[[All, {2, 3, 1}]], ConstantArray[".", {nsnp + 1, 6}], data[[All, 4 ;;]], 2];
        data[[1, ;; 9]] = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"};
        data[[2 ;;, 2]] = Round[data[[2 ;;, 2]]/recomrate 10^6];
        data[[2 ;;, 9]] = genoformat;
        data[[2 ;;, 4]] = "A";
        data[[2 ;;, 5]] = "G";
        data[[2 ;;, ;; 3]] = Map[ToString, data[[2 ;;, ;; 3]], {2}];
        vcfheader = {"##fileformat=VCFv4.1", "##source=magicImpute", 
          "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", 
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">"};
        data = Join[vcfheader, data];
        Export[outputfile,data, "TSV"]
    ]
    
getvcfsampleName[vcffile_String] :=
    Module[ {vcfdata, linefield},
        If[ ! FileExistsQ[vcffile],
            Print["getvcfsampleName: vcffile =", vcffile, " does not exist!"];
            Abort[];
        ];
        vcfdata = ReadList[vcffile, String];
        linefield = 
         Position[
           vcfdata, _?(StringTake[#, 2] =!= "##" && 
               StringTake[#, 1] == "#" &), {1}, 1, Heads -> False][[1, 1]];
        Drop[StringSplit[vcfdata[[linefield]], "\t"], 9]
    ]
      
      
genoofVCF[vcffile_String?FileExistsQ,genoformat:"GT"|"AD"|"GP"] :=
    Module[ {vcfdata, linefield, depth, nind, pos, ls,i,bool},
        vcfdata = ReadList[vcffile, String];
        linefield = Position[vcfdata, _?(StringTake[#, 2] =!= "##" && 
               StringTake[#, 1] == "#" &), {1}, 1, Heads -> False][[1, 1]];
        (*columns of depth: {CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,SampleID..}*)
        (*FORMAT, e.g. GT:AD:DP:GQ:PL, where GT=genotype, AD=allelic depth, DP =read depth, GQ=genotype quality, PL=Normalized, Phred-scaled likelihoods*)
        depth = StringSplit[#, "\t"] & /@ vcfdata[[linefield ;;]];
        depth[[2 ;;, ;; 2]] = ToExpression[depth[[2 ;;, ;; 2]]];
        nind = Length[depth[[1, 10 ;;]]];
        depth[[2 ;;, 9 ;;]] = Map[StringSplit[#, ":"] &, depth[[2 ;;, 9 ;;]]];
        depth[[2 ;;, 9]] = Flatten[Position[#, genoformat, {1}, 1, Heads -> False]] & /@ depth[[2 ;;, 9]];
        bool = Min[#] & /@ Map[Length, depth[[2 ;;, 10 ;;]], {2}];
        bool = Thread[bool >= depth[[2 ;;, 9, 1]]];
        If[ (And @@ Union[bool]),
            Do[depth[[i, 10 ;;]] = depth[[i, 10 ;;, depth[[i, 9, 1]]]], {i, 2, Length[depth]}],
            Print["genoofVCF: genotype ",genoformat," data are missing!"];
            Abort[]
        ];
        If[ genoformat==="GT",
            depth[[2 ;;, 10 ;;]] = Map[StringSplit[#, "/" | "|"] &,depth[[2 ;;, 10 ;;]], {2}];
            depth[[2 ;;, 10 ;;]] = Map[StringJoin @@ # &, depth[[2 ;;, 10 ;;]]/. {"."->"N","0" -> "1", "1" | "2" | "3" -> "2"}, {2}],
            (*genoformat==="AD"*)
            depth[[2 ;;, 10 ;;]] = Map[StringSplit[#, ","] &,Replace[depth[[2 ;;, 10 ;;]], {"." -> "0,0"}, {2}],{2}];
            (*depth eleement: {read count for reference allele, and total read count for alternative alleles}*)
            pos = Position[depth[[2 ;;, 10 ;;]], _?(Length[#] >= 3 &), {2}, Heads -> False];
            If[ pos =!= {},
                ls = ToExpression[Extract[depth[[2 ;;, 10 ;;]], pos]];
                ls = Transpose[{ls[[All, 1]], Total[ls[[All, 2 ;;]], {2}]}];
                ls = Map[ToString, ls, {2}];
                depth[[2 ;;, 10 ;;]] = ReplacePart[depth[[2 ;;, 10 ;;]], Thread[pos -> ls]];
            ];
            depth[[2 ;;, 10 ;;]] = Map[StringJoin @@ Riffle[#, "|"] &, depth[[2 ;;, 10 ;;]], {2}];
        ];
        depth[[2 ;;, 9]] = genoformat;
        depth = Join[{First[depth]}, Flatten[SortBy[#, #[[2]] &] & /@ SortBy[SplitBy[depth[[2 ;;]], #[[1]] &], #[[1, 1]] &], 1]];
        (*fill missing marker id*)
        pos = Flatten[Position[depth[[2 ;;, 3]], "."]];
        depth[[pos + 1, 3]] = renameDuplicates[Map[StringJoin["SNP", #[[1]], "_", #[[2]]] &, Map[ToString, depth[[pos + 1, 1 ;; 2]], {2}]]];
        (*keep only some columns, and tranpose them*)
        depth = Transpose[depth[[All, Join[{3, 1, 2}, Range[nind] + 9]]]];
        (*change from ID to MARKER, so that Excel does import it as a SYLK file.*)
        depth[[1,1]] = "MARKER";
        depth
    ]
          
genoofVCF2[vcffile_String?FileExistsQ,genoformat:"GT"|"AD"] :=
    Module[ {vcfdata, linefield, depth, nind, ad, pos, ls,i},
        vcfdata = ReadList[vcffile, String];
        linefield = Position[vcfdata, _?(StringTake[#, 2] =!= "##" && 
               StringTake[#, 1] == "#" &), {1}, 1, Heads -> False][[1, 1]];
        (*columns of depth: {CHROM, POS, ID, REF, ALT, QUAL, FILTER, INFO, FORMAT,SampleID..}*)
        (*FORMAT, e.g. GT:AD:DP:GQ:PL, where GT=genotype, AD=allelic depth, DP =read depth, GQ=genotype quality, PL=Normalized, Phred-scaled likelihoods*)
        depth = StringSplit[#, "\t"] & /@ vcfdata[[linefield ;;]];
        depth[[2 ;;, ;; 2]] = ToExpression[depth[[2 ;;, ;; 2]]];
        nind = Length[depth[[1, 10 ;;]]];
        depth[[2 ;;, 9 ;;]] = Map[StringSplit[#, ":"] &, depth[[2 ;;, 9 ;;]]];
        depth[[2 ;;, 9]] = Flatten[Position[#, genoformat, {1}, 1, Heads -> False]] & /@ depth[[2 ;;, 9]];
        If[ genoformat==="GT",
            Do[ad = Table[{"N", "N"}, {nind}];
               If[ MatchQ[depth[[i, 9]], {_}],
                   pos = Flatten[Position[Length[#] & /@depth[[i, 10 ;;]], _?(# >= depth[[i, 9, 1]] &), {1},Heads -> False]];
                   ad[[pos]] = StringSplit[depth[[i, pos + 9, depth[[i, 9, 1]]]] /. {"./." | ".|." -> "N/N"},"/" | "|"] /. {"0" -> "1", "1"|"2"|"3" -> "2"}
               ];
               depth[[i, 10 ;;]] = ad, {i, 2, Length[depth]}];
            depth[[2 ;;, 10 ;;]] = Map[StringJoin @@ # &, depth[[2 ;;, 10 ;;]], {2}],
            (*genoformat==="AD"*)
            Do[
             ad = Table[{"0", "0"}, {nind}];
             If[ MatchQ[depth[[i, 9]], {_}],
                 pos = Flatten[Position[Length[#] & /@depth[[i, 10 ;;]], _?(# >= depth[[i, 9, 1]] &), {1},Heads -> False]];
                 ad[[pos]] = StringSplit[depth[[i, pos + 9, depth[[i, 9, 1]]]] /. {"." -> "0,0"}, ","]
             ];
             depth[[i, 10 ;;]] = ad, {i, 2, Length[depth]}];
            (*depth eleement: {read count for reference allele, and total read count for alternative alleles}*)
            pos = Position[depth[[2 ;;, 10 ;;]], _?(Length[#] >= 3 &), {2}, Heads -> False];
            If[ pos =!= {},
                ls = ToExpression[Extract[depth[[2 ;;, 10 ;;]], pos]];
                ls = Transpose[{ls[[All, 1]], Total[ls[[All, 2 ;;]], {2}]}];
                ls = Map[ToString, ls, {2}];
                depth[[2 ;;, 10 ;;]] = ReplacePart[depth[[2 ;;, 10 ;;]], Thread[pos -> ls]];
            ];
            depth[[2 ;;, 10 ;;]] = Map[StringJoin @@ Riffle[#, "|"] &, depth[[2 ;;, 10 ;;]], {2}];
        ];
        depth[[2 ;;, 9]] = genoformat;
        depth = Join[{First[depth]}, Flatten[SortBy[#, #[[2]] &] & /@ SortBy[SplitBy[depth[[2 ;;]], #[[1]] &], #[[1, 1]] &], 1]];
        (*fill missing marker id*)
        pos = Flatten[Position[depth[[2 ;;, 3]], "."]];
        depth[[pos + 1, 3]] = renameDuplicates[Map[StringJoin["SNP", #[[1]], "_", #[[2]]] &, Map[ToString, depth[[pos + 1, 1 ;; 2]], {2}]]];
        (*keep only some columns, and tranpose them*)
        depth = Transpose[depth[[All, Join[{3, 1, 2}, Range[nind] + 9]]]];
        (*change from ID to MARKER, so that Excel does import it as a SYLK file.*)
        depth[[1,1]] = "MARKER";
        depth
    ]


(*recombinationRate: centiMorgan per Mb, and it is used to transfer genetic distance (cM) 
of magicad into physical distance (bp) of vcf*)

Options[magicSNPtoLBImpute] = {
    recombinationRate -> 1,
    outputFileID -> ""
    }
magicSNPtoLBImpute[magicadfile_String?FileExistsQ,OptionsPattern[]] :=
    Module[ {outputfile,recomrate,outputid,magicsnp, nind, nsnp, data, callgt, vcfheader},
        {recomrate,outputid} = OptionValue[{recombinationRate,outputFileID}];
        If[ outputid=!="",
            outputfile = outputid <> ".vcf",
            outputfile = "magicSNPtoLBImpute_Output.vcf";
        ];
        magicsnp = Import[magicadfile, "CSV", Path -> Directory[]];
        {nind, nsnp} = Dimensions[Rest[magicsnp]] - {3, 1};
        data = Transpose[Rest[magicsnp]];
        data = Join[data[[All, {2, 3, 1}]], ConstantArray[".", {nsnp + 1, 6}], data[[All, 4 ;;]], 2];
        data[[1, ;; 9]] = {"#CHROM", "POS", "ID", "REF", "ALT", "QUAL","FILTER", "INFO", "FORMAT"};
        data[[2 ;;, 2]] = Round[data[[2 ;;, 2]]/recomrate 10^6];
        data[[2 ;;, 9]] = "GT:AD:GQ";
        data[[2 ;;, 4]] = "A";
        data[[2 ;;, 5]] = "G";
        data[[2 ;;, ;; 3]] = Map[ToString, data[[2 ;;, ;; 3]], {2}];
        callgt = Sign[ToExpression[StringSplit[#, "|"] & /@ data[[2 ;;, 10 ;;]]]];
        callgt = callgt /. {{0, 0} -> "./.", {1, 0} -> "0/0", {0, 1} -> "1/1", {1, 1} -> "0/1"};
        data[[2 ;;, 10 ;;]] = MapThread[StringJoin[#1, ":", #2, ":."] &, {callgt,Map[StringReplace[#, "|" -> ","] &, data[[2 ;;, 10 ;;]], {2}]},2];
        data[[2 ;;, 10 ;;]] = data[[2 ;;, 10 ;;]] /. {"./.:0,0:." -> "./."};
        vcfheader = {"##fileformat=VCFv4.1", "##source=magicImpute", 
          "##FORMAT=<ID=AD,Number=.,Type=Integer,Description=\"Allelic depths for the ref and alt alleles in the order listed\">", 
          "##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">",
          "##FORMAT=<ID=GQ,Number=1,Type=Integer,Description=\"Genotype Quality\">"};
        data = Join[vcfheader, data];
        Export[outputfile, data,"TSV"]
    ]
    
magicSNPtompMap[inputmagicsnp_, outputid_String,chrset_:All,snpset_:All] :=
    Module[ {magicsnp=inputmagicsnp, nfounder, rule, mapdata, founderdata, offdata, outputid2 = outputid},        
        If[ StringQ[magicsnp],
            If[ FileExistsQ[magicsnp],
                magicsnp = Import[magicsnp,"CSV"],
                Print["File ",magicsnp," does not exist!"];
                Abort[];
            ];
        ];
        magicsnp = getsubMagicSNP[magicsnp,chrset,snpset];
        nfounder = magicsnp[[1, 2]];
        rule = {1 | "1" -> 0, 2 | "2" -> 1, 11 | "11" -> 0, 22 | "22" -> 1, 12|21|"12"|"21"->"NA", "N" | "NN" | "1N" | "2N" -> "NA"};
        magicsnp[[5 ;;, 2 ;;]] = magicsnp[[5 ;;, 2 ;;]] /. rule;
        mapdata = Transpose[magicsnp[[2 ;; 4]]];
        founderdata = Join[magicsnp[[{2}]], magicsnp[[5 ;; nfounder + 4]]];
        founderdata[[1,1]] = "";
        offdata = Join[magicsnp[[{2}]], magicsnp[[nfounder + 5 ;;]]];
        offdata[[1,1]] = "";
        If[ outputid2 =!= "",
            outputid2 = outputid2 <> "_"
        ];
        {Export[outputid2 <> "founders.txt", founderdata, "TSV"],
         Export[outputid2 <> "finals.txt", offdata, "TSV"],
         Export[outputid2 <> "map.txt", mapdata, "TSV"]}
    ]
      
mpMaptoMagicSNP[mapMapfilestem_] :=
    Module[ {founderdata, nfounder, offdata, magicsnp},
        founderdata = Import[mapMapfilestem <> ".founder.csv"];
        nfounder = Length[founderdata] - 1;
        offdata = Import[mapMapfilestem <> ".ril.csv"];
        (*remove first column of inserted random phenotypes*)
        offdata = offdata[[All, 2 ;;]];
        magicsnp = Join[offdata[[;; 3]], 
          founderdata[[2 ;;]] /. {"NA" -> "N", 0 -> 1, 1 -> 2}, 
          offdata[[4 ;;]] /. {"NA" -> "NN", 0 -> 11, 1 -> 22}];
        magicsnp[[;; 3, 1]] = {"marker", "chromosome", "cM"};
        magicsnp = Join[{{"nfounder", nfounder}}, magicsnp];
        Export[mapMapfilestem <> "_magicSNP.csv", magicsnp]
    ]
        
genotoMagicSNP[inputdepth_, recomfreq_?NumericQ,founders_List,isfounderinbred_,outputfile_String] :=
    Module[ {magicsnp = inputdepth, pos, ls,nfounder},
        (*recomfreq in unit of cM/Mb*)
        magicsnp[[3, 2 ;;]] *= recomfreq 10^(-6.);
        magicsnp[[3, 1]] = "cM";
        pos = Flatten[Position[magicsnp[[4 ;;, 1]], #, {1}, 1, Heads -> False]] & /@ 
          founders;
        ls = Pick[founders, pos, {}];
        If[ ls =!= {},
            Print["genotoMagicSNP: no founders ", ls];
            Abort[]
        ];
        If[ Length[Flatten[pos]] > Length[pos],
            Print["genotoMagicSNP: founders are not unique!"];
            Abort[];
        ];
        pos = pos[[All, 1]];
        pos = Join[pos, DeleteCases[Range[Length[magicsnp] - 3], Alternatives @@ pos]];
        pos = Join[Range[3], pos + 3];
        magicsnp = Join[{{"nfounder", Length[founders]}}, magicsnp[[pos]]];
        If[ isfounderinbred,
            nfounder = Length[founders];
            magicsnp[[5 ;; nfounder + 4, 2 ;;]] = 
             magicsnp[[5 ;; nfounder + 4, 2 ;;]] /. {11|"11"| "1N" | "N1" -> 1,22|"22"| "2N" | "N2" -> 2, 12 |21|"12"|"21"| "NN" -> "N"}
        ];
        csvExport[outputfile,magicsnp]
    ]
        
caldepthfreq[allelicdepth_] :=
    Module[ {ls},
        ls = Tally[Flatten[allelicdepth]];
        ls[[All, 1]] = Total[ToExpression[StringSplit[ls[[All, 1]], "|"]], {2}];
        ls = {#[[1, 1]], Total[#[[All, 2]]]} & /@ SplitBy[SortBy[ls, First], First];
        ls[[All,2]] = Normalize[ls[[All,2]],Total];
        SortBy[ls, First]
    ]

calMagicReadDepth[inputmagicsnp_,isprint_:True] :=
    Module[ {magicsnp = inputmagicsnp,nfounder, freq,meandepth,lab,indices},
        If[ StringQ[magicsnp],
            If[ FileExistsQ[magicsnp],
                magicsnp = Import[magicsnp,"CSV"],
                Print["File ",magicsnp," does not exist!"];
                Abort[];
            ];
        ];
        nfounder = magicsnp[[1, 2]];
        freq = Table[caldepthfreq[magicsnp[[indices, 2 ;;]]], {indices, {5 ;; nfounder + 4, nfounder + 5 ;; All}}];
        meandepth = N[Dot @@ Transpose[#] & /@ freq];
        If[ isprint,
            lab = {"founder", "offspring"};
            Print[Table[
                ListPlot[freq[[i]], AxesLabel -> {"Depth", "Freq"}, PlotRange -> All, 
                                    ImageSize -> 300, Filling -> Axis, 
                                    PlotLabel -> "Mean " <> lab[[i]] <> " read depth = " <> ToString[meandepth[[i]]]], {i, Length[freq]}]];
        ];
        {meandepth,freq}
    ]

calMagicMissingFraction[inputmagicsnp_, includehalfmissing_: True] :=
    Module[ {magicsnp = inputmagicsnp,nfounder, ls, fracfounder, fracoff},
        If[ StringQ[magicsnp],
            If[ FileExistsQ[magicsnp],
                magicsnp = Import[magicsnp,"CSV"],
                Print["File ",magicsnp," does not exist!"];
                Abort[];
            ];
        ];
        nfounder = magicsnp[[1, 2]];
        ls = Flatten[magicsnp[[5 ;; nfounder + 4, 2 ;;]]];
        fracfounder = N[(Count[ls, "NN" | "N" | "0|0"] + Boole[includehalfmissing] Count[ls, "1N" | "2N" | "N1" | "N2"]/2)/Length[ls]];
        ls = Flatten[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        fracoff = N[(Count[ls, "NN" | "N" | "0|0"] +Boole[includehalfmissing] Count[ls, "1N" | "2N" | "N1" | "N2"]/2)/Length[ls]];
        {fracfounder, fracoff}
    ]
  
randomMagicRead[snpfile_?FileExistsQ, preadlist_List] :=
    Module[ {starttime, magicsnp, nfounder,indices,adsnp, pread, res, countls, rules, 
      outputfile,i},
        starttime = SessionTime[];
        If[ ! (And @@ Thread[Flatten[preadlist] <= 1]),
            Print["Sampling probably of reads must be <=1."];
            Abort[]
        ];
        magicsnp = Import[snpfile, Path -> Directory[]];
        nfounder = magicsnp[[1,2]];
        indices = {5;;nfounder+4,nfounder+5;;All};
        Table[
         PrintTemporary["Time elapsed = " <> ToString[Round[SessionTime[] - starttime, 0.1]] <> " second. Binomial sampling probability of reads: ", pread];
         res = magicsnp;
         If[ !(And@@Thread[pread==1]),
             Do[
                 adsnp = Map[StringSplit[#, "|"] &, magicsnp[[indices[[i]], 2 ;;]]];
                 countls = DeleteCases[Union[Flatten[adsnp]], "0"];
                 rules = MapThread[#1 :> Random[BinomialDistribution[#2, pread[[i]]]] &, {countls,ToExpression[countls]}];
                 res[[indices[[i]], 2 ;;]] = Replace[adsnp, rules, {3}];
                 0,{i,Length[indices]}];
             res[[5 ;;, 2 ;;]] = toDelimitedString[res[[5 ;;, 2 ;;]], "|"];
         ];
         outputfile = StringDrop[snpfile,-4]<>"_pread" <> ToString[N[pread[[1]]]] <> "-" <> ToString[N[pread[[2]]]] <> ".csv";
         csvExport[outputfile, res], {pread, preadlist}]
    ]
    
    
rawGenotypeCall[allelicdepth_?(MatrixQ[#, ArrayQ[#, 2] &] &), qualityscore_, callthreshold_] :=
    Module[ {ad},
        Table[rawGenotypeCall[ad, qualityscore, callthreshold], {ad, allelicdepth}]
    ]

rawGenotypeCall[allelicdepth_?(VectorQ[#, ArrayQ[#, 2] &] &), qualityscore_, callthreshold_] :=
    Module[ {res, span},
        res = rawGenotypeCall[Flatten[allelicdepth, 1], qualityscore, callthreshold];
        span = Accumulate[Prepend[Length[#] & /@ allelicdepth, 0]];
        span = Thread[span[[;; -2]] + 1 ;; span[[2 ;;]]];
        res[[#]] & /@ span
    ]
      
rawGenotypeCall[allelicdepth_?(ArrayQ[#, 2] &), qualityscore_, callthreshold_] :=
    Module[ {baseerrorprob, n1, n2, ntot, p11, p12, p22, sum, pos, p1, p2,res},
        baseerrorprob = 10^(-qualityscore/10.);
        {n1, n2} = Transpose[allelicdepth];
        ntot = n1 + n2;
        p12 = (1/2)^ntot;
        p11 = (1 - baseerrorprob)^n1 baseerrorprob^n2;
        p22 = (baseerrorprob)^n1 (1 - baseerrorprob)^n2;
        sum = p11 + 2 p12 + p22;
        (*assuming that 11,12,21,22 are equal probable*)
        {p11, p12, p22} = #/sum & /@ {p11, p12, p22};
        res = ConstantArray["NN", Length[allelicdepth]];
        res[[Flatten[Position[Thread[p11 > callthreshold], True]]]] = "11";
        res[[Flatten[Position[Thread[p22 > callthreshold], True]]]] = "22";
        res[[Flatten[Position[Thread[2 p12 > callthreshold], True]]]] = "12";
        pos = Flatten[Position[res, "NN"]];
        {p11, p12, p22} = {p11, p12, p22}[[All, pos]];
        {p1, p2} = {p11 + 2 p12, p22 + 2 p12};
        res[[Pick[pos, Thread[And[Thread[p1 > callthreshold], Thread[p1 > p2]]]]]] = "1N";
        res[[Pick[pos, Thread[And[Thread[p2 > callthreshold], Thread[p2 > p1]]]]]] = "2N";
        res
    ]

inbredRawGenotypeCall[allelicdepth_?ArrayQ] :=
    Replace[Sign[allelicdepth], {{0, 0} -> "N", {1, 0} -> "1", {0, 1} ->"2", {1, 1} -> "N"}, {Depth[allelicdepth] - 2}]
    
magicInbredGenotypeCall[magicsnp_] :=
    magicRawGenotypeCall[magicsnp, None,None,True,True]

magicRawGenotypeCall[magicsnp_, qualityscore_, callthreshold_,isfounderinbred_, isoffspringinbred_] :=
    Module[ {thresholds,callsnp = magicsnp, nfounder, ls, i},
        thresholds = Flatten[{callthreshold}];
        nfounder = callsnp[[1, 2]];
        ls = callsnp[[5 ;;, 2 ;;]];
        Monitor[
         Do[ls[[i]] = ToExpression[StringSplit[ls[[i]], "|"]], {i, 
           Length[ls]}], ProgressIndicator[i, {1, Length[ls]}]];
        If[ isfounderinbred,
            ls[[;; nfounder]] = inbredRawGenotypeCall[ls[[;; nfounder]]],
            ls[[;; nfounder]] = rawGenotypeCall[ls[[;; nfounder]], qualityscore, First[thresholds]]
        ];
        If[ isoffspringinbred,
            ls[[nfounder + 1 ;;]] = inbredRawGenotypeCall[ls[[nfounder + 1 ;;]]]/.{"1"->"11","2"->"22","N"->"NN"},
            ls[[nfounder + 1 ;;]] = rawGenotypeCall[ls[[nfounder + 1 ;;]], qualityscore, Last[thresholds]]
        ];
        callsnp[[5 ;;, 2 ;;]] = ls;
        callsnp
    ]
      
(*magicsnp and R/HAPPY inputs*)
(*
PrintTemporary["Importing data chr = " <> ToString[chr]];
        allelefile = "chr" <> ToString[chr] <> ".MAGIC.alleles";
        hapdatafile = "chr" <> ToString[chr] <> ".MAGIC.data";
        *)
gethappyinput[{allelefile_, datafile_}] :=
    Module[ {alleledata, happydata},
        alleledata = Import[allelefile, "Table"];
        alleledata = Split[alleledata, {#1[[1]], #2[[1]]} == {"marker", 
              "allele"} || {#1[[1]], #2[[1]]} == {"allele", "allele"} &];
        happydata = Import[datafile, "Table"];
        happydata = Join[#[[;; 6]], Partition[#[[7 ;;]], 2]] & /@ happydata;
        {alleledata, happydata}
    ]
    
alleletogeno[x_] :=
    {x[[All, 2]],
    Total[Table[x[[i, 3 ;;]] /. {_?Positive -> i, _?PossibleZeroQ -> 0}, {i,Length[x]}]] /. {0 -> "N"}}

happytomagic[{alleledata_, happydata_}] :=
    Module[ {allele, nmark, nstrain, strains, snp, allelecode, fgeno, 
      obsdata, code, magicSNP, ls, genders, pos,mk},
        allele = alleledata;
        {nmark, nstrain} = allele[[1, 1, {2, 4}]];
        strains = allele[[2, 1, 2 ;;]];
        allele = Drop[allele, 2];
        snp = allele[[All, 1]];
        If[ Length[snp] != nmark,
            Print["Inconsistent number of markers!"];
            Beep[]
        ];
        snp = Prepend[snp[[All, {2, 4, 5}]], {"SNP", "Chromosome", "cM"}];
        allele = allele[[All, 2 ;;]];
        allele = DeleteCases[#, _?(#[[2]] == "NA" &)] & /@ allele;
        pos = Position[ArrayQ[#] & /@ allele, False];
        If[ pos =!= {},
            Print["bad SNPs: ", snp[[Flatten[pos] + 1]]];
            Abort[];
        (*allele=Delete[allele,pos];
            snp=Delete[snp,pos+1];
            nmark=Length[allele];*)
        ];
        {allelecode, fgeno} = Transpose[alleletogeno[#] & /@ allele];
        If[ Dimensions[fgeno] != {nmark, nstrain},
            Print["Inconsistent number of strains!"];
            Beep[]
        ];
        fgeno = Prepend[fgeno, strains];
        fgeno = Join[Transpose[snp[[All, ;; 2]]], Transpose[fgeno]];
        (*extract obsdata*)
        obsdata = Transpose[happydata[[All, 7 ;;]]];
        code = Map[Range[Length[#]] &, allelecode];
        ls = Table[obsdata[[mk]] /. Thread[allelecode[[mk]] -> code[[mk]]], {mk,Length[code]}];
        ls = Map[ToString, ls, {3}];
        ls = Map[StringJoin @@ # &, ls, {2}] /. {"NANA" -> "NN"};
        ls = Transpose[Prepend[ls, happydata[[All, 2]]]];
        If[ Union[snp[[2 ;;, 2]]] === {"X"},
            genders = happydata[[All, 5]];
            pos = Flatten[Position[genders, "M"]];
            ls[[pos, 2 ;;]] = Map[StringTake[#, 1] &, ls[[pos, 2 ;;]], {2}];
        ];
        obsdata = Join[fgeno[[;; 2]], ls];
        (*return*)
        allelecode = Join[Transpose[{snp[[All, 1]]}], Prepend[allelecode, Range[2]], 2];
        magicSNP = Prepend[Join[Transpose[snp], fgeno[[3 ;;]], 
           obsdata[[3 ;;]]], {"#founders", Length[fgeno] - 2}];
        {magicSNP, allelecode}
    ]

happytoMagicSNP[happyallelefile_,happydatafile_,outputid_String] :=
    happytoMagicSNP[{{happyallelefile,happydatafile}},outputid]
    
happytoMagicSNP[happyinputfiles_?(ArrayQ[#,2,StringQ]&),outputid_String] :=
    Module[ {res,temp,magic, code},
        If[ ! (And @@ (Length[#] == 2 & /@ ConstantArray[0, {3, 2}])),
            Print["R/HAPPY input files must be a list of {allelefile, datafffile}!"];
            Abort[]
        ];
        temp = Flatten[happyinputfiles];
        temp = Pick[temp,FileExistsQ[#]&/@temp,False];
        If[ temp=!={},
            Print["Files ",temp, " do not exist!"];
            Abort[];
        ];
        res = Map[happytomagic[gethappyinput[#]]&,happyinputfiles];
        magic = res[[All, 1, 2 ;;]];
        magic = Join[magic[[1]], Sequence @@ magic[[2 ;;, All, 2 ;;]], 2];
        magic = Prepend[magic, res[[1, 1, 1]]];
        code = res[[All, 2]];
        code = Join[code[[1]], Sequence @@ code[[2 ;;, 2 ;;]]];
        If[ ! (code[[2 ;;, 1]] == magic[[2, 2 ;;]]),
            Print["Inconsistent marker id between code and magic!"]
        ];
        {Export[outputid<>"_magicSNP.csv", magic],Export[outputid<>"_allelecode.csv", code]}
    ]
    
getJoinMapCPInput[jointmapinputfile_?FileExistsQ] :=
    Module[ {data, indid},
        data = Import[jointmapinputfile, "Table"];
        indid = data[[2, 2 ;;]];
        data = Transpose[data[[3 ;;, 2 ;;]]];
        data[[2 ;; 3]] = Transpose[StringSplit[StringDelete[data[[2]], "<" | ">"], "x"]];
        data[[2 ;;]] = data[[2 ;;]] /. {"--" -> "NN", "hh" -> "11", "hk" -> "12", "kk" -> "22", 
                                        "ll" -> "11", "lm" -> "12", "nn" -> "11","np" -> "12"};
        data = Insert[data, Table["NA", {Length[data[[1]]]}], {{2}, {2}}];
        Join[{{"nfounder", 2}},Join[List /@ Join[{"marker", "chromosome", "cM", "P1", "P2"}, indid], data, 2]]
    ]    

(***************************************************************************************************************)
          
      
End[] (* End Private Context *)

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicDataPreprocess`*"];

EndPackage[]

