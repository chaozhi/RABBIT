(* Mathematica Package *)

(* Created by the Wolfram Workbench Jan 22, 2017 *)

BeginPackage["MagicDataFilter`",{"MagicDefinition`","MagicDataPreprocess`"}]

magicsnpBinning::usage = "magicsnpBinning  "

missingFractionFilter::usage = "missingFractionFilter  "

distortionFilterCP::usage = "distortionFilterCP  "

distortionFilterRIL::usage = "distortionFilterRIL  "

distortionTest::usage = "distortionTest  "

readDuplicatefile::usage = "readDuplicatefile  "

binningDuplicate::usage = "binningDuplicate  "

saveBinnedSNP::usage = "saveBinnedSNP  "

cliquepartition::usage = "cliquepartition  "

magicAllelicCode::usage = "magicAllelicCode  "

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package *)

magicAllelicCode[magicsnp_] :=
    Module[ {nfounder, geno, n22, n12, n11, n2, n1, afreq, bool, pos, res},
        nfounder = magicsnp[[1, 2]];
        geno = Transpose[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        n22 = Count[#, 22] & /@ geno;
        n11 = Count[#, 11] & /@ geno;
        n12 = Count[#, 12 | 21] & /@ geno;
        n1 = Count[#, 1 | "1N"] & /@ geno;
        n2 = Count[#, 2 | "2N"] & /@ geno;
        (*afreq= allelic frequency of allele 2*)
        afreq = Round[N[(2 n22 + n12 + n2)/(2 n22 + 2 n12 + 2 n11 + n2 + n1)],1./nfounder];
        (*set 2 as minor alllele, and set the genotype of founder = 1 or 11 if allelefreq of 2 = 0.5*)
        bool = # == 2 || # == 22 & /@ magicsnp[[5, 2 ;;]];
        pos = Flatten[Position[Boole[bool] Boole[Thread[afreq == 0.5]]+Boole[Thread[afreq > 0.5]], 1]];
        res = magicsnp;
        res[[5 ;;, pos + 1]] = res[[5 ;;, pos + 1]] /. {1 -> 2, 2 -> 1, 11 -> 22, 22 -> 11, 
           12 -> 21, 21 -> 12, "2N" -> "1N", "N2" -> "N1", "1N" -> "2N","N1" -> "N2"};
        res
    ]

Options[magicsnpBinning] = {
  outputFileID -> "",
  isPrintTimeElapsed -> True,
  isRunInParallel -> True
  }

magicsnpBinning[inputmagicSNP_?(ListQ[#] || StringQ[#] &), OptionsPattern[]] :=
    Module[ {magicSNP = inputmagicSNP,magicSNP00,outputid, isprint, isparallel,starttime,
       nsnp, binls,i,outputfiles},
        {outputid, isprint, isparallel} = OptionValue@{outputFileID, isPrintTimeElapsed, isRunInParallel};
        If[ isprint,
            starttime = SessionTime[];
            Print["magicsnpBinning. Start date = ", DateString[]];
        ];
        If[ StringQ[inputmagicSNP],
            If[ !FileExistsQ[inputmagicSNP],
                Print["File ", inputmagicSNP," does not exist!"];
                Return[$Failed]
            ];
            magicSNP = Import[inputmagicSNP,"CSV"];
        ];
        (*magicSNP00 saves original magicSNP, and it will be used for outputing*)
        magicSNP00 = magicSNP;
        magicSNP=magicAllelicCode[magicSNP];
        magicSNP[[3, 2 ;;]] = "NA";
        (*possible genotypes: "N","1","2", "11","12","22","1N","2N","NN"*)
        Do[magicSNP[[i, 2 ;;]] = Map[ToString, magicSNP[[i, 2 ;;]]] /. {"21" -> "12","N1"|"N2"->"NN"},{i,5,Length[magicSNP]}];
        nsnp = Length[magicSNP[[2]]]-1;
        binls = List /@ Range[nsnp];
        (*Put[{magicSNP, binls, isparallel},"D:\\Chaozhi\\MMA Workspace\\MagicDataFilterApp\\temp.txt"];*)
        binls = getSNPbin[magicSNP, binls,isparallel,starttime];
        outputfiles = saveBinnedSNP[magicSNP00, binls, outputid];
        If[ isprint,
            Print["Number of SNPs remain: ",Length[binls]];
            Print["Done. Finished date = ", DateString[], ". Time elapsed = ", Round[SessionTime[] - starttime, 0.01], " seconds."];
        ];
        outputfiles
    ]
    
getrowspan[nsnp_, nseg_] :=
    Module[ {ls, ls2},
        ls = Accumulate[Range[nsnp, 1, -1]];
        ls2 = Round[Range[nseg - 1] Last[ls]/nseg];
        ls = Union[Join[Flatten[Position[ls - #, _?Positive, {1}, 1, Heads -> False] & /@ls2], {1, nsnp+1}]];
        Thread[Most[ls] ;; Rest[ls] - 1]
    ]
          
getSNPbin[magicSNP_,binlist_,isparallel_,starttime_] :=
    Module[ {binlist2 = binlist,bins},
        If[ binlist2 ==Automatic,
            binlist2 = List /@ Range[Length[magicSNP[[2, 2 ;;]]]]
        ];
        bins = If[ Intersection[{"1N", "2N"}, Union[Flatten[magicSNP[[5 ;;, 2 ;;]]]]] ==={},
                   fastgetbin[magicSNP,binlist2,isparallel],
                   0
               ];
        PrintTemporary["{#bin_before,#bin_after} = ", {Length[binlist2],Length[bins]}, ". Time elapsed = ", Round[SessionTime[] - starttime, 0.01]," seconds."];
        bins
    ]    
    
  
findsnpduplicate[snpgeno_, snpidls_, rowspan_, outputfile_] :=
    Module[ {nsnp, outstream, i, ii,jj, ls, res,nmiss},
        nsnp = Dimensions[snpgeno][[2]];
        Quiet[Close[outputfile]];
        Put[outputfile];
        outstream = OpenWrite[outputfile];
        ii = Range[nsnp - 1][[rowspan]];
        If[ First[ii] == 1,
            nmiss = Count[#, 0] & /@ Transpose[snpgeno];
            Write[outstream, Join[{{"Marker-ID","#Missing"}},Transpose[{snpidls,nmiss}]]];
        ];
        Do[
         ls = Sign[snpgeno[[All, i]]*snpgeno[[All, i + 1 ;;]]];
         ls *= Sign[Abs[snpgeno[[All, i]] - snpgeno[[All, i + 1 ;;]]]];
         ls = 1 - Sign[Total[ls]];
         jj = Pick[Range[i + 1, nsnp], ls, 1];
         res = Thread[{i, jj}];
         If[ res =!= {},
             Write[outstream, res]
         ], {i, ii}];
        Close[outstream];
    ]
  
parallelfastgetbin[snpgeno_, snpidls_, outputfile_] :=
    Module[ {nsnp, spanls, filels, countls, count,outstream, i},
        nsnp = Dimensions[snpgeno][[2]];
        spanls = getrowspan[nsnp-1, Max[1, IntegerLength[nsnp-1] - 2] $KernelCount];
        filels = "temporaryfile" <> ToString[#] <> ".txt" & /@ Range[Length[spanls]];
        filels[[1]] = outputfile;
        filels = FileNameJoin[{Directory[], #}] & /@ filels;
        countls = Range[nsnp - 1, 1, -1];
        countls = Total[countls[[#]]] & /@ spanls;
        SetSharedVariable[count];
        count = 0;
        Monitor[ParallelDo[
          findsnpduplicate[snpgeno, snpidls, spanls[[i]], filels[[i]]];
          count += countls[[i]], {i, Length[filels]}, 
          DistributedContexts -> {"MagicDataFilter`", 
            "MagicDataFilter`Private`"}, Method -> "FinestGrained"], 
         ProgressIndicator[count, {0, Total[countls]}]];
        Quiet[Close[filels[[1]]]];
        outstream = OpenAppend[filels[[1]]];
        Do[Write[outstream, #] & /@ ReadList[filels[[i]]], {i, 2,Length[filels]}];
        Close[outstream];
        DeleteFile[#] & /@ Rest[filels];
    ]    

fastgetbin[magicSNP_, binlist_, isparallel_] :=
    Module[ {magicsnp, snpgeno, snpidls,nsnp,outputfile = "temporaryfile_snpduplicate.txt"},
        magicsnp = getsubMagicSNP[magicSNP, All, binlist[[All, 1]]];
        snpgeno = magicsnp[[5 ;;, 2 ;;]] /. {"N" | "NN" | "1N" | "2N" | "N1" | "N2" -> 0, "11" | "1" -> 1, "22" | "2" -> 2, "12" -> 3};
        snpidls = magicsnp[[2, 2;;]];
        nsnp = Length[magicsnp[[2]]]-1;
        If[ isparallel,
            parallelfastgetbin[snpgeno, snpidls, outputfile],
            findsnpduplicate[snpgeno, snpidls, 1;;nsnp - 1, outputfile]
        ];
        binningDuplicate[FileNameJoin[{Directory[], outputfile}]]
    ]   
    
toknngraph[g_, knn_: Automatic] :=
    Module[ {n, k, rules, rules2, ad, ad2},
        If[ knn === Automatic,
            n = VertexCount[g];
            k = Ceiling[Sqrt[n]],
            k = Ceiling[knn];
        ];
        ad = AdjacencyMatrix[g];
        rules = ArrayRules[ad];
        rules2 = SplitBy[rules, #[[1, 1]] &];
        rules2 = Join[Flatten[Take[#, UpTo[k]] & /@ Most[rules2], 1], rules2[[-1]]];
        ad2 = SparseArray[rules2, {n, n}];
        ad2 = Sign[ad2 + Transpose[ad2]];
        AdjacencyGraph[VertexList[g], ad2]
    ]    
    
cliquepartition0[g_] :=
    Module[ {clu,g2,i},
        clu = FindGraphCommunities[g, Method -> "Modularity"];
        clu = Table[
            g2 = Subgraph[g,clu[[i]]];
            If[ VertexCount[g2]>1000,
                g2 = toknngraph[g2]
            ];
            FindGraphCommunities[g2,Method -> "CliquePercolation"],{i,Length[clu]}];
        clu = Flatten[clu,1];
        Sort[#]&/@clu
    ]
    
cliquepartition[g_] :=
    Module[ {ls},
        If[ ConnectedGraphQ[g],
            cliquepartition0[g],
            ls = Flatten[cliquepartition0[#] & /@ ConnectedGraphComponents[g], 1];
            (*each community is represented by the smallest node, 
            and commnunities with shared smallest nodes are merged*)
            ls = SplitBy[SortBy[ls, First], First];
            Union[Flatten[#]] & /@ ls
        ]
    ]

readDuplicatefile[duplicatefile_?FileExistsQ] :=
    Module[ {res,snpidls,nmiss,rule,nsnp},
        res = ReadList[duplicatefile];
        {snpidls, nmiss} = Transpose[Rest[First[res]]];
        nsnp = Length[snpidls];
        res = Flatten[Rest[res], 1];
        (*rule = Join[{{i_, i_} -> 0}, Thread[res -> 1]];*)
        rule = Thread[res -> 1];
        res = SparseArray[rule, {nsnp, nsnp}];
        res += Transpose[res];
        {res,nmiss}
    ]
  
binningDuplicate[duplicatefile_?FileExistsQ] :=
    Module[ {res, nmiss,order, newres, gg, group},
        {res,nmiss} = readDuplicatefile[duplicatefile];
        order = Ordering[nmiss];
        newres = res[[order, order]];
        gg = AdjacencyGraph[newres];
        group = cliquepartition[gg];
        order[[#]] & /@ group
    ]


saveBinnedSNP[magicSNP_,binlist_,outputid_String] :=
    Module[ {bins,cols,magicSNP2 = magicSNP, outputfiles,snpid, binindicator,selectedbin, clusters,j},
        bins = SortBy[binlist, First];
        cols = Join[{1}, 1 + bins[[All, 1]]];
        magicSNP2 = Join[magicSNP[[{1}]], magicSNP[[2 ;;, cols]]];
        snpid = magicSNP[[2, 2 ;;]];
        selectedbin = binindicator = ConstantArray[0, {Length[snpid]}];
        Do[binindicator[[bins[[j]]]] = j;
           selectedbin[[bins[[j, 1]]]] = j, {j, Length[bins]}];
        clusters = Join[{{"Marker-ID", "InputMarkerNo","ReplicateGroup", "Representive"}}, Transpose[{snpid, Range[Length[snpid]],binindicator, selectedbin}]];
        outputfiles = {"magicsnpBinning.csv","binning.csv"};
        If[ outputid=!="",
            outputfiles = outputid<>"_"<>#&/@outputfiles;
        ];
        {csvExport[outputfiles[[1]], magicSNP2],csvExport[outputfiles[[2]], clusters]}
    ]

(***************************************************************************************************************)
     
missingFractionFilter[magicsnp_, maxmissingfraction_, inputtarget_String] :=
    Module[ {target, missingcode, nfounder, noff, nsnp, miss, pos,newmagicsnp},
        missingcode = "1N" | "N1" | "2N" | "N2" | "NN" | "N";
        nfounder = magicsnp[[1, 2]];
        {noff, nsnp} = Dimensions[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        target = ToLowerCase[inputtarget];
        (*"Individual"|"Marker"];*)
        Switch[target,
         "individual",
         miss = Count[#, missingcode] & /@ magicsnp[[nfounder + 5 ;;, 2 ;;]];
         miss = Transpose[{Range[noff], magicsnp[[nfounder + 5 ;;, 1]], miss,nsnp - miss, N[miss/nsnp], Table[0, {noff}]}];
         miss[[All, 6]] = Thread[miss[[All, 5]] > maxmissingfraction];
         miss = SortBy[miss, -#[[5]] &];
         miss = Join[{{"Individual No.", "Individual ID", "#Missing Genotypes", 
             "#Non-missing Genotypes", "Missing fraction", "isDeleted"}}, miss];
         pos = Sort[Pick[miss[[2 ;;]], miss[[2 ;;, 6]], False][[All, 1]]];
         Print[Take[miss, UpTo[Min[20,(noff - Length[pos] + 5)]]] // MatrixForm];
         newmagicsnp = Join[magicsnp[[;; nfounder + 4]],magicsnp[[nfounder + 5 ;;]][[pos]]];
         {newmagicsnp, miss},
         "marker",
         miss = Count[#, missingcode] & /@Transpose[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
         miss = Transpose[{Range[nsnp], magicsnp[[2, 2 ;;]], miss, noff - miss,N[miss/noff], Table[0, {nsnp}]}];
         miss[[All, 6]] = Thread[miss[[All, 5]] > maxmissingfraction];
         miss = SortBy[miss, -#[[5]] &];
         miss = Join[{{"Marker No.", "Marker ID", "#Missing Genotypes", 
             "#Non-missing Genotypes", "Missing fraction", "isDeleted"}},miss];
         pos = Sort[Pick[miss[[2 ;;]], miss[[2 ;;, 6]], False][[All, 1]]];
         Print["Number of deleted individuals: ", nsnp - Length[pos], 
          " with fraction of missing genotype >", maxmissingfraction, "!"];
         Print[Take[miss, UpTo[Min[20,(nsnp - Length[pos] + 5)]]] // MatrixForm];
         newmagicsnp = Join[magicsnp[[{1}]], magicsnp[[2 ;;, Join[{1}, pos + 1]]]];
         {newmagicsnp, miss},
         _,
         Print["missingFractionFilter: wrong target ", target, 
          "! target must be individual or marker!"];
         Abort[]
         ]
    ]
       
(*missingFractionFilter is selection markers that either 1) no missing parental genotypes or 2) missing offspring genotype<maxmissingfraction*)
missingFractionFilter[magicADcsvfile_?FileExistsQ,qualityscore_, genothreshold_, maxmissingfraction_] :=
    Module[ {magicsnp, nfounder, ls0,ls, pos, outfileid,outputfile},
        magicsnp = Import[magicADcsvfile, "CSV"];
        nfounder = magicsnp[[1, 2]];
        Print["Before filtering magicsnp dimensions: ",Dimensions[Rest[magicsnp]]];
        ls0 = magicRawGenotypeCall[magicsnp, qualityscore, genothreshold,False,False];
        ls = Transpose[ls0[[nfounder+5;;, 2 ;;]]];
        ls = (Count[#, "NN" | "1N" | "2N"] & /@ ls)/Last[Dimensions[ls]];
        pos = Flatten[Position[ls - maxmissingfraction, _?Negative]];
        ls = Transpose[ls0[[5 ;; nfounder + 4, 2 ;;]]];
        ls = Count[#, "N"|"NN" | "1N" | "2N"] & /@ ls;
        pos = Union[pos, Flatten[Position[ls, 0]]];
        magicsnp = Join[magicsnp[[{1}]], magicsnp[[2 ;;, Join[{1}, pos + 1]]]];
        Print["After filtering magicsnp dimensions: ", Dimensions[Rest[magicsnp]]];
        outfileid = "_missingFilter" <> ToString[AccountingForm[maxmissingfraction,10]];
        outputfile = StringDrop[magicADcsvfile, -4]<> outfileid<> ".csv";
        csvExport[outputfile, magicsnp]
    ]    
                
distortionTest[genocount_, genoratio_, siglevel_] :=
    Module[ {ratio, df, samplesize, expect, chi2stat, crit},
        ratio = Normalize[genoratio, Total];
        samplesize = Total[genocount, {2}];
        expect = Transpose[samplesize # & /@ ratio];
        chi2stat = Total[(genocount - expect)^2/expect, {2}];
        df = Length[ratio] - 1;
        crit = InverseCDF[ChiSquareDistribution[df], 1 - siglevel];
        Thread[# > crit] & /@ chi2stat
    ]

getParentGenopairs[segretype_String] :=
    Module[ {geno, rule, genopair, genopair2, type, set},
        geno = {"11", "12", "1N", "22", "2N", "NN"};
        rule = Dispatch[{"11" -> {"11"}, "12" -> {"12"}, "22" -> {"22"}, 
          "1N" -> {"11", "12"}, "2N" -> {"12", "22"}, 
          "NN" -> {"11", "12", "22"}}];
        genopair = Flatten[Outer[List, geno, geno], 1];
        genopair2 = Replace[genopair, rule, {2}];
        genopair2 = Flatten[Outer[List, Sequence @@ #], 1] & /@ genopair2;
        genopair2 = Map[Sort, genopair2, {2}];
        genopair2 = Union[#] & /@ genopair2;
        genopair2 = Transpose[{genopair, genopair2}];
        type = Switch[segretype, "S110", {"11", "12"}, "S011", {"12", "22"},
           "S121", {"12", "12"}, _, Abort[]];
        set = Select[genopair2, MemberQ[#[[2]], type] &][[All, 1]];
        {set, Complement[genopair, set]}
    ]
        
distortionFilterCP[inputmagicsnp_?(ListQ[#] || StringQ[#] &), distortionsiglevel_,genoerrorprob_:0.01] :=
    Module[ {magicsnp = inputmagicsnp, nfounder, ls, rule, count, ratiolist,logratiolist, posterprob, segretypepos, newsegretypepos,res, 
       indicator, label, absentgeno, pos, pos2, set,i, offgeno, parentgeno,isimputeparent = True},
        If[ StringQ[magicsnp],
            If[ !FileExistsQ[magicsnp],
                Print["File ", magicsnp," does not exist!"];
                Return[$Failed]
            ];
            magicsnp = Import[magicsnp,"CSV"];
        ];
        magicsnp[[5 ;;, 2 ;;]] = Map[ToString, magicsnp[[5 ;;, 2 ;;]], {2}];
        nfounder = magicsnp[[1, 2]];
        ls = Transpose[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        ls = DeleteCases[#, "NN" | "1N" | "2N" | "N1" | "N2"] & /@ ls;
        ls = ls /. {"21" -> "12"};
        rule = Tally[#] & /@ ls;
        rule = Map[Rule @@ # &, rule, {2}];
        count = ({"11", "12", "22"} /. #) & /@ rule;
        count = count /. {"22" -> 0, "12" -> 0, "11" -> 0};
        (*ordering of segration types in ratiolist is important for the following calculation!*)
        ratiolist = {{1, 1, 0}, {0, 1, 1}, {1, 2, 1}, {1, 0, 0}, {0,1,0}, {0, 0, 1}};
        ratiolist = Normalize[#, Total] & /@ ratiolist;
        ratiolist = # (1 - genoerrorprob) + (1 - #) genoerrorprob/2 & /@ ratiolist;
        logratiolist = Transpose[Log[ratiolist]];
        posterprob = Exp[# - Max[#] & /@ (count.logratiolist)];
        posterprob = Round[Normalize[#, Total] & /@ posterprob, 10^(-5.)];
        segretypepos = Flatten[Ordering[#, -1] & /@ posterprob];
        segretypepos = Table[Flatten[Position[segretypepos, i]], {i, 3}];
        segretypepos =  Table[If[ segretypepos[[i]] === {},
                                  {},
                                  Pick[segretypepos[[i]], distortionTest[count[[segretypepos[[i]]]], ratiolist[[i]], distortionsiglevel], False]
                              ], {i, 3}];
        ls = Length[#] & /@ segretypepos;
        Print[ls, " SNPs out of ", Length[magicsnp[[2]] - 1], 
            " have been identified to be one of segregation types {110, 011, 121} without distortion!"];
        (*Correcting for incompatible offspring allelic depths for segregation types 110 and 011*)
        res = indicator = magicsnp;
        indicator[[5 ;;, 2 ;;]] = 0;
        label = {"S110", "S011", "S121"};
        absentgeno = {"22", "11", "No"};
        offgeno = Transpose[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        pos = Table[
          pos = segretypepos[[i]];
          If[ pos === {},
              {},
              pos2 = Position[offgeno[[pos]], absentgeno[[i]], {2}];
              pos2[[All, 1]] = pos[[pos2[[All, 1]]]];
              Print[Length[pos2], " out of ", Last[Dimensions[offgeno]] Length[pos], 
               " offspring allelic depths that are called to be genotype ", 
               absentgeno[[i]], " are incompatible to segregation type ", 
               label[[i]], ", and they are set to be missing!"];
              pos2
          ], {i, 2}];
        pos = Join @@ pos;
        If[ pos =!= {},
            pos = Reverse[Transpose[pos]];
            pos = Transpose[pos + {nfounder + 4, 1}];
            res = ReplacePart[res, Dispatch[Thread[pos -> "NN"]]];
            indicator = ReplacePart[indicator, Dispatch[Thread[pos -> 1]]];
        ];
        
        (*Removing markers with incompatible parental allelic depths for segregation types 110,011,and 121*)
        parentgeno = Transpose[magicsnp[[5 ;; nfounder + 4, 2 ;;]]];
        newsegretypepos = Table[
           pos = segretypepos[[i]];
           set = Last[getParentGenopairs[label[[i]]]];
           Pick[pos, MemberQ[set, #] & /@ parentgeno[[pos]], False], {i, 3}];
        MapThread[Print[Length[#1] - Length[#2], " out of ", Length[#1]," SNPs at which parental allelic depths are incompatible with segregregation type ", #3, " are removed!"] &, {segretypepos,newsegretypepos, label}];
        segretypepos = newsegretypepos;
        ls = Length[#] & /@ segretypepos;
        Print[ls, " SNPs have been identified to be one of segregation types {110, 011, 121} after removing incompatibile markers!"];
                
        (*Keep only snps with identified segregation types*)
        If[ isimputeparent,
            (*imputing parental genotypes*)
            (*{"S110", "S011", "S121"}*)
            parentgeno = Transpose[magicsnp[[5 ;; nfounder + 4, 2 ;;]]];
            parentgeno[[segretypepos[[1]]]] = Replace[parentgeno[[segretypepos[[1]]]], 
                Dispatch[{{"11", "1N" | "2N" | "NN"} -> {"11", "12"}, {"1N" | "2N" | "NN","11"} -> {"12", "11"},
                 {"12", "1N" | "NN"} -> {"12","11"}, {"1N" | "NN", "12"} -> {"11", "12"}, 
                 {"1N", "2N"} -> {"11", "12"}, {"2N", "1N"} -> {"12", "11"}}], {1}];
            parentgeno[[segretypepos[[2]]]] = Replace[parentgeno[[segretypepos[[2]]]], 
                Dispatch[{{"22","1N" | "2N" | "NN"} -> {"22", "12"}, {"1N" | "2N" | "NN","22"} -> {"12", "22"}, 
                 {"12", "2N" | "NN"} -> {"12","22"}, {"2N" | "NN", "12"} -> {"22", "12"}, 
                 {"1N","2N"} -> {"12", "22"}, {"2N", "1N"} -> {"22", "12"}}], {1}];
            parentgeno[[segretypepos[[3]]]] = Table[{"12", "12"}, {Length[segretypepos[[3]]]}];
            res[[5 ;; nfounder + 4, 2 ;;]] = Transpose[parentgeno];
        ];
        pos = Union[Flatten[segretypepos]];
        res = Join[res[[{1}]], res[[2 ;;, Join[{1}, pos + 1]]]];
        indicator = Join[indicator[[{1}]], indicator[[2 ;;, Join[{1}, pos + 1]]]];
        {res,indicator}
    ]   
    
distortionFilterCP[magicADcsvfile_?FileExistsQ, qualityscore_, callthreshold_, distortionsiglevel_,isoutputcalled_] :=
    Module[ {callthreshold2,magicsnpad,callsnp, nfounder, ls, rule, count, ratiolist,callerrorprob,
       logratiolist, posterprob, segretypepos, newsegretypepos,res, indicator, label, absentgeno, pos, pos2, set,i, offad, parentad,outfileid,outfiles},
        magicsnpad = Import[magicADcsvfile, "CSV"];
        callthreshold2 = Flatten[{callthreshold}];
        callsnp = magicRawGenotypeCall[magicsnpad, qualityscore, callthreshold2,False,False];
        nfounder = callsnp[[1, 2]];
        ls = Transpose[callsnp[[nfounder + 5 ;;, 2 ;;]]];
        ls = DeleteCases[#, "NN" | "1N" | "2N" | "N1" | "N2"] & /@ ls;
        ls = ls /. {"21" -> "12"};
        rule = Tally[#] & /@ ls;
        rule = Map[Rule @@ # &, rule, {2}];
        count = ({"11", "12", "22"} /. #) & /@ rule;
        count = count /. {"22" -> 0, "12" -> 0, "11" -> 0};
        (*ordering of segration types in ratiolist is important for the following calculation!*)
        ratiolist = {{1, 1, 0}, {0, 1, 1}, {1, 2, 1}, {1, 0, 0}, {0,1,0}, {0, 0, 1}};
        ratiolist = Normalize[#, Total] & /@ ratiolist;
        callerrorprob = 1 - Last[callthreshold];
        ratiolist = # (1 - callerrorprob) + (1 - #) callerrorprob/2 & /@ ratiolist;
        logratiolist = Transpose[Log[ratiolist]];
        posterprob = Exp[# - Max[#] & /@ (count.logratiolist)];
        posterprob = Round[Normalize[#, Total] & /@ posterprob, 10^(-5.)];
        segretypepos = Flatten[Ordering[#, -1] & /@ posterprob];
        segretypepos = Table[Flatten[Position[segretypepos, i]], {i, 3}];
        segretypepos = Table[Pick[segretypepos[[i]],distortionTest[count[[segretypepos[[i]]]], ratiolist[[i]], distortionsiglevel], False], {i, 3}];
        ls = Length[#] & /@ segretypepos;
        Print[ls, " SNPs (in total ",Total[ls],") have been identified to be one of segregation types {110, 011, 121} without distortion!"];    
            
        (*Correcting for incompatible offspring allelic depths for segregation types 110 and 011*)
        res = indicator = If[ isoutputcalled,
                              callsnp,
                              magicsnpad
                          ];
        indicator[[5 ;;, 2 ;;]] = 0;
        label = {"S110", "S011", "S121"};
        absentgeno = {"22", "11", "No"};
        offad = Transpose[callsnp[[nfounder + 5 ;;, 2 ;;]]];
        pos = Table[
          pos = segretypepos[[i]];
          pos2 = Position[offad[[pos]], absentgeno[[i]], {2}];
          pos2[[All, 1]] = pos[[pos2[[All, 1]]]];
          Print[Length[pos2], " out of ", Last[Dimensions[offad]] Length[pos], 
           " offspring allelic depths that are called to be genotype ", 
           absentgeno[[i]], " are incompatible to segregation type ", 
           label[[i]], ", and they are set to be missing!"];
          pos2, {i, 2}];
        pos = Reverse[Transpose[Join @@ pos]];
        pos = Transpose[pos + {nfounder + 4, 1}];
        If[ isoutputcalled,
            res = ReplacePart[res, Dispatch[Thread[pos -> "NN"]]],
            res = ReplacePart[res, Dispatch[Thread[pos -> "0|0"]]];
        ];
        indicator = ReplacePart[indicator, Dispatch[Thread[pos -> 1]]];        
        
        (*Removing markers with incompatible parental allelic depths for segregation types 110,011,and 121*)
        parentad = Transpose[callsnp[[5 ;; nfounder + 4, 2 ;;]]];
        newsegretypepos = Table[
           pos = segretypepos[[i]];
           set = Last[getParentGenopairs[label[[i]]]];
           Pick[pos, MemberQ[set, #] & /@ parentad[[pos]], False], {i, 3}];
        MapThread[Print[Length[#1] - Length[#2], " out of ", Length[#1]," SNPs at which parental allelic depths are incompatible with segregregation type ", #3, " are removed!"] &, {segretypepos,newsegretypepos, label}];
        segretypepos = newsegretypepos;
        ls = Length[#] & /@ segretypepos;
        Print[ls, " SNPs (in total ", Total[ls],") have been identified to be one of segregation types {110, 011, 121} after removing incompatibile markers!"];
                
        (*Keep only snps with identified segregation types*)
        pos = Union[Flatten[segretypepos]];
        res = Join[res[[{1}]], res[[2 ;;, Join[{1}, pos + 1]]]];
        indicator = Join[indicator[[{1}]], indicator[[2 ;;, Join[{1}, pos + 1]]]];
        If[ isoutputcalled,
            (*imputing parental genotypes*)
            (*{"S110", "S011", "S121"}*)
            parentad = Transpose[callsnp[[5 ;; nfounder + 4, 2 ;;]]];
            parentad[[segretypepos[[1]]]] = Replace[parentad[[segretypepos[[1]]]], 
                Dispatch[{{"11", "1N" | "2N" | "NN"} -> {"11", "12"}, {"1N" | "2N" | "NN","11"} -> {"12", "11"},
                 {"12", "1N" | "NN"} -> {"12","11"}, {"1N" | "NN", "12"} -> {"11", "12"}, 
                 {"1N", "2N"} -> {"11", "12"}, {"2N", "1N"} -> {"12", "11"}}], {1}];
            parentad[[segretypepos[[2]]]] = Replace[parentad[[segretypepos[[2]]]], 
                Dispatch[{{"22","1N" | "2N" | "NN"} -> {"22", "12"}, {"1N" | "2N" | "NN","22"} -> {"12", "22"}, 
                 {"12", "2N" | "NN"} -> {"12","22"}, {"2N" | "NN", "12"} -> {"22", "12"}, 
                 {"1N","2N"} -> {"12", "22"}, {"2N", "1N"} -> {"22", "12"}}], {1}];
            parentad[[segretypepos[[3]]]] = Table[{"12", "12"}, {Length[segretypepos[[3]]]}];
            res[[5 ;; nfounder + 4, 2 ;;]] = Transpose[parentad];
        ];
        outfileid = "_distortionFilter" <>ToString[AccountingForm[distortionsiglevel, 10]];
        outfileid = {outfileid,outfileid<>"_CorrectingIndicator"};
        outfiles = StringDrop[magicADcsvfile, -4]<> #<> ".csv"&/@outfileid;
        {csvExport[outfiles[[1]], res],csvExport[outfiles[[2]], indicator]}
    ]    
    
(*Assume founders are exchangeable*)
distortionFilterRIL[inputmagicsnp_?(ListQ[#] || StringQ[#] &), sdsiglevel_,genoerrorprob_:0.01] :=
    Module[ {magicsnp = inputmagicsnp, nfounder, founderfreq,i,distortionpos, 
      offspringcount, samplesize, expect, chi2stat, crit,ratiolist,logratiolist,posterprob,maxsegretype,monomorphicpos,keeppos},
        If[ StringQ[magicsnp],
            If[ !FileExistsQ[magicsnp],
                Print["File ", magicsnp," does not exist!"];
                Return[$Failed]
            ];
            magicsnp = Import[magicsnp,"CSV"];
        ];
        magicsnp[[5 ;;, 2 ;;]] = Map[ToString, magicsnp[[5 ;;, 2 ;;]], {2}];
        nfounder = magicsnp[[1, 2]];
        offspringcount = DeleteCases[#, "NN" | "12" | "21"] & /@ 
           Transpose[magicsnp[[nfounder + 5 ;;, 2 ;;]]];
        offspringcount = {Count[#, "11"], Count[#, "22"]} & /@ offspringcount;
        samplesize = Total[offspringcount, {2}];
        If[ sdsiglevel == 0,
            distortionpos = {},
            founderfreq = Transpose[magicsnp[[5 ;; nfounder + 4, 2 ;;]]];
            founderfreq = ToExpression[Distribute[# /. {"N" -> {"1", "2"}}, List] & /@ founderfreq];
            founderfreq -= 1;
            founderfreq = Total[founderfreq, {3}]/nfounder;
            (*to exclude frequency of monomorphic SNPs*)
            founderfreq = Union[DeleteCases[#, 0 | 1]] & /@ founderfreq;
            founderfreq = Map[{1 - #, #} &, founderfreq, {2}];
            Quiet[chi2stat = Table[
                expect = (samplesize[[i]] founderfreq[[i]]);
                Total[(offspringcount[[i]] - #)^2/#] & /@ expect, {i,Length[offspringcount]}];];
              (*Indeterminate\[Rule]monomorphic markers*)
            chi2stat = chi2stat /. {Indeterminate -> Infinity, ComplexInfinity -> Infinity};
            chi2stat = Min[#] & /@ chi2stat;
            crit = InverseSurvivalFunction[ChiSquareDistribution[1], sdsiglevel];
            distortionpos = Flatten[Position[Thread[# >= crit] & /@ chi2stat, True]];
        ];
        If[ TrueQ[genoerrorprob==0],
            monomorphicpos = {},
            (*drop monomorphic SNPs*)
            ratiolist = Table[{i/nfounder, 1 - i/nfounder}, {i, 0, nfounder}];
            ratiolist = # (1 - genoerrorprob) + (1 - #) genoerrorprob/2 & /@ ratiolist;
            logratiolist = Transpose[Log[ratiolist]];
            posterprob = Exp[# - Max[#] & /@ (offspringcount.logratiolist)];
            posterprob = Round[Normalize[#, Total] & /@ posterprob, 10^(-5.)];
            maxsegretype = Flatten[Ordering[#, -1] & /@ posterprob];
            monomorphicpos = Flatten[Position[maxsegretype, 1 | (nfounder + 1)]];
        ];
        keeppos = Complement[Range[Length[offspringcount]], distortionpos,monomorphicpos];
        If[ TrueQ[genoerrorprob==0],
            Print[Length[distortionpos], " out of ", Length[offspringcount], " SNPs are detected to be under segregation distortion at significance level = ", sdsiglevel],
            Print[Length[monomorphicpos], " out of ", Length[offspringcount], " SNPs are detected to be monomorphic assuming genotype error probability = ", genoerrorprob, "."];
            Print[Length[distortionpos], " out of ", Length[offspringcount], " SNPs are detected to be under segregation distortion at significance level = ", sdsiglevel, 
                ". Among them ",Length[Intersection[distortionpos, monomorphicpos]], " SNPs are also monomporhic."];
        ];
        (*Print["In total, ", Length[offspringcount] - Length[keeppos], " out of ", Length[offspringcount], " SNPs are dropped!"];*)
        {Join[magicsnp[[{1}]], magicsnp[[2 ;;, Join[{1}, 1 + keeppos]]]], keeppos, chi2stat, crit}
    ]
    
(***************************************************************************************************************)    



End[]

EndPackage[]

