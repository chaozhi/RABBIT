
(* Created by the Wolfram Workbench Nov 30, 2015 *)

(* :Author: Chaozhi Zheng <chaozhi@gmail.com>*)
(* :Description: A package for simulating genotypic data (SNPArry or sequcing-based read data) in experimental populations*)
(* :KeyFunctions: simPedigreeInfor, magicGeneDropping*)

BeginPackage["MagicGeneDropping`",{"MagicDefinition`","MagicSimulate`"}]
(* Exported symbols added here with SymbolName::usage *) 

magicGeneDropping::usage = "magicGeneDropping [founderhaplo, popdesign] simulates genomic data in the population with the breeding design specified by popdesign. The founderhaplo specifies the founder haplotypes and genetic map over a set of markers or the input filename. The popdesign speficies the breeding design information in several possible ways: a list of mating schemes from founder population to the last generation, a list of values denoting the junction distribution, or a filename for population pedigree information. "

founderReadDepth::usage = "founderReadDepth is an option of magicGeneDropping to specify the type of founder genotypes. A postive value denotes the coverage depth of sequencing-based-genotyping and genotypes are given in terms of allelic depths; The value 0 denotes the SNP array genotypes."

offspringReadDepth::usage = "offspringReadDepth is an option of magicGeneDropping to specify the type of offspring genotypes. A postive value denotes the coverage depth of sequencing-based-genotyping and genotypes are given in terms of allelic depths; The value 0 denotes the SNP array genotypes."

founderGenotypeMissing::usage = "founderGenotypeMissing is an option of magicGeneDropping to specify the fraction of missing genotypes among founders. "

offspringGenotypeMissing::usage = "offspringGenotypeMissing is an option of magicGeneDropping to specify the fraction of missing genotypes among offspring. "

phredQualScore::usage = "phredQualScore is an option of magicGeneDropping to specify the Phred quality score"

interferenceStrength::usage = "interferenceStrength is an option to specifiy the strength of chiasma interference. "

isObligate::usage = "isObligate is an option to specify whether there must be at least one crossover event. "

sampleGender::usage = "sampleGender is an option to specify the gender of samples. The option value must be \"Hermaphrodite\",\"Female\", \"Male\", or \"Random\". "

funnelFunction::usage = "funnelFunction is an option to specify the fuction of defining funnel code for each sample. "

depthDistribution::usage = "depthDistribution is an option to specify the probability of read depth (total read count at a site), and it must be \"Poisson\" or \"NegativeBinomial-n\", where n is a positive shape parameter!"

dropMonomorphic::usage = "dropMonomorphic is an option to specify whether to drop monomorphic markers. "

splitPedigreeInfor::usage = "splitPedigreeInfor[pedigreeInfor] splits the pedigreeInform into description, pedigree, and sampleInfor. Here sampleInfor specifies for each sample the memeberID in pedigree and funnel code. "

mergePedigreeInfor::usage = "mergePedigreeInfor[pedigree, sampleInfor] merges pedigree and sampleInfor  into pedigreeInfor. "

simPedigreeInfor::usage = "simPedigreeInfor[pedigree] returns pedigreeInfor with sampleInfor being simulated from pedigree. "

expandPedigreeInfor::usage = "expandPedigreeInfor returns a single large pedigree by combining pedgree and sampleInfor. "

Begin["`Private`"]
(* Implementation of the package *)

sameSetQ[list1_, list2_] :=
    Complement[list1, list2] === Complement[list2, list1] === {}   
    
(*apply only to funneled based pedigreeInfor with inbred founders and discrete generations; unordered pair of {MotherID, FatherID} in 4th column.*)
(*pedigree members are labeled by natural number starting from 1...*)
expandPedigreeInfor[pedigreeinfor_] :=
    Module[ {pedigree, sampleinfor, funnelcode, noffspring, nfounder,pedsize, res, ls, rule,i},
        {pedigree, sampleinfor} = Rest[splitPedigreeInfor[pedigreeinfor]];
        funnelcode = sampleinfor[[2 ;;, 3]];
        {noffspring, nfounder} = Dimensions[funnelcode];
        pedsize = Length[pedigree] - 1;
        res = Table[pedigree, {noffspring}];
        Do[
         res[[i, 2 ;;, {2, 4}]] =res[[i, 2 ;;, {2, 4}]] /.Thread[Range[nfounder] -> funnelcode[[i]]];
         res[[i, 2 ;;, {2, 4}]] = res[[i, 2 ;;, {2, 4}]] /.Thread[Range[nfounder + 1, pedsize] -> 
             Range[pedsize (i - 1) + nfounder + 1, pedsize i]], {i, noffspring}];
        ls = res;
        Do[ls[[i]] = SplitBy[Rest[ls[[i]]], First], {i, Length[res]}];
        ls = Flatten[ls[[All, 2, All, {2, 4}]], 1];
        ls[[All, 2]] = Sort[#] & /@ ls[[All, 2]];
        ls = SortBy[#, First] & /@ GatherBy[ls, Last];
        ls = ls[[All, All, 1]];
        rule = Flatten[Thread[Rest[#] -> First[#]] & /@ ls];
        res[[All, 2 ;;, {2, 4}]] = res[[All, 2 ;;, {2, 4}]] /. rule;
        res = Join[res[[1]], SortBy[Flatten[res[[2 ;;, (nfounder + 2) ;;]], 1], First]];
        res[[2 ;;, 4]] = Sort[#] & /@ res[[2 ;;, 4]];
        res = Join[res[[{1}]], SortBy[DeleteDuplicates[res[[2 ;;]]],First]];
        res[[2 ;; nfounder + 1]] = SortBy[res[[2 ;; nfounder + 1]], #[[2]] &];
        res[[2 ;;, {2, 4}]] = res[[2 ;;, {2, 4}]] /.Thread[res[[2 ;;, 2]] -> Range[Length[res] - (1)]];
        res
    ]    
    
splitPedigreeInfor[pedigreeInfor_] :=
    Module[ {key,infor,description,pedigree, sampleInfor,nFounder,founders},
        key = pedigreeInfor[[1, 1]];
        infor = Partition[Split[pedigreeInfor, #1[[1]] != key && #2[[1]] != key &], 2];
        description = infor[[All,1]];
        {pedigree, sampleInfor} = infor[[All, 2]];
        pedigree = Join[pedigree[[All, ;; -3]], Transpose[{pedigree[[All, -2 ;;]]}], 2];
        sampleInfor[[2 ;;, -1]] = ToString[#] & /@ sampleInfor[[2 ;;, -1]];
        sampleInfor[[2 ;;, -1]] = StringSplit[#, "-"] & /@ sampleInfor[[2 ;;, -1]];
        nFounder = Count[pedigree[[2 ;;, -1]], {0, 0}];
        founders = Map[ToString,pedigree[[2 ;; nFounder + 1, 2]]];
        If[ ! (And @@ Union[sameSetQ[founders, #] & /@ sampleInfor[[2 ;;, -1]]]),
            Print["splitPedigreeInfor: inconsistent founder IDs between pedigree and funnelcode of pedigreeInfor!"];
            Abort[]
        ];
        sampleInfor[[2 ;;, -1]] = sampleInfor[[2 ;;, -1]] /. Thread[founders -> Range[nFounder]];
        {description, pedigree, sampleInfor}
    ]

mergePedigreeInfor[pedigree_,sampleInfor_,key_:"Pedigree-Information" ] :=
    Module[ {infor = sampleInfor,nfounder},
        nfounder = Max[infor[[2, -1]]];
        infor[[2 ;;, -1]] = StringJoin[Riffle[#, "-"]] & /@ Map[ToString,infor[[2 ;;, -1]],{2}];
        Join[{{key, "The in-depth pedigree with funnelcode 1-2..."<>ToString[nfounder]}}, Flatten[#] & /@ pedigree, 
         {{key, "The memberID and funnelcode for each sampled individual"}},infor]
    ] 
    
Options[simPedigreeInfor] = {
    sampleSize -> 100,
    sampleGender->"Hermaphrodite",
    funnelFunction ->(RandomSample[#]&),
    outputFileID -> ""
}          
       
simPedigreeInfor[pedigree_List, OptionsPattern[]] :=
    Module[ {samplesize,samplegender,funnelfun,outputid,gender,idlist, pos, memberlist, nFounder,funnelcode, sampleinfor,pedinfor},
        {samplesize,samplegender,funnelfun,outputid} = OptionValue@{sampleSize,sampleGender,funnelFunction,outputFileID};
        If[ outputid=!="",
            outputid = outputid<>"_"
        ];
        gender = Switch[ToLowerCase[samplegender],
                    "hermaphrodite"|"notavailable", 0,
                    "female", 1,
                    "male", 2,
                    "random", _,
                    _, Print["simPedigreeInfor: sampleGender must be \"Hermaphrodite\",\"Female\", \"Male\", or \"Random\"!"];
                       Abort[]];
        nFounder = Count[pedigree[[2 ;;, -1]], {0, 0}];
        idlist = "ProgenyLine" <> ToString[#] & /@ Range[samplesize];
        pos = Flatten[Position[pedigree[[2 ;;, 1]], pedigree[[-1, 1]]]]+1;
        pos = Pick[pos, pedigree[[pos, 3]], gender];
        If[ pos==={},
            Print["simPedigreeInfor: there are no individuals having the given gender in the last generation!"];
            Abort[];
        ];
        If[Length[pos]>=samplesize,
        	pos = Sort[RandomSample[pos, samplesize]],
        	pos = RandomChoice[pos, samplesize];
        ];
        memberlist = pedigree[[pos, 2]];
        (*funnelcode = Table[RandomSample[Range[nFounder]], {samplesize}];*)
        (*funnelcode = (SplitBy[SortBy[Transpose[{Range[nFounder], pedigree[[2 ;; nFounder + 1, 3]]}], Last], Last])[[All, All, 1]];
        funnelcode = Table[Flatten[RandomSample[#] & /@ funnelcode], {samplesize}][[All, Flatten[funnelcode]]],*)
        funnelcode = Range[nFounder];
        funnelcode = Table[funnelfun[funnelcode], {samplesize}];
        sampleinfor = Join[{{"ProgenyLine", "MemberID", "Funnelcode"}}, Transpose[{idlist, memberlist, funnelcode}]];
        pedinfor = mergePedigreeInfor[pedigree, sampleinfor];
        Export[outputid<> "PedigreeInfor.csv",pedinfor,"CSV"]
    ]

simPedigreeInfor[popcode_String, opts : OptionsPattern[]] :=
 Module[{samplegender, isoogamy, ped},
  samplegender = OptionValue[sampleGender];
  isoogamy = !MatchQ[ToLowerCase[samplegender],"hermaphrodite"|"notavailable"];
  ped = compilePopPedigree[popcode, isoogamy];
  simPedigreeInfor[ped, opts]
  ]
  
simPedigreeInfor[inipop_List, mateScheme_, opts:OptionsPattern[]] :=
    simPedigreeInfor[simPedigree[inipop, mateScheme], opts]   
    
simPedigreeInfor[nFounder_Integer,mateScheme_, opts:OptionsPattern[]] :=
    Module[ {samplegender,isoogamy},
        samplegender = OptionValue[sampleGender];
        isoogamy = !MatchQ[ToLowerCase[samplegender],"hermaphrodite"|"notavailable"];
        simPedigreeInfor[setIniPop[nFounder,isoogamy],mateScheme, opts]
    ]
    
checkfounderHaplo[simSetup_, founderHaplo_] :=
    Module[ {popped,iniPopFGL,nFounder,chrLength,temp},
        {popped, iniPopFGL} = simSetup[[{1,2}]];
        chrLength = iniPopFGL[[-1, -1, All, 1, -1, 1]];
        nFounder = Count[popped[[2 ;;, 4]], {0, 0}];
        (*nFgl = Length[Union[Flatten[iniPopFGL[[2 ;;, -1, 1, All, 1, 2]]]]];*)
        If[ TrueQ[nFounder != Length[founderHaplo] - 3],
            Print["checkfounderHaplo: The values for the number of Founders between popFgl and founderHaplo are not consistent!"];
            Abort[];
        ];
        temp = #[[{1, -1}, 2]] & /@ SplitBy[Transpose[founderHaplo[[2 ;; 3, 2 ;;]]], First];
        If[ Length[temp] != Length[chrLength],
            Print["checkfounderHaplo: The values for the number of chromosomes between popFgl and founderHaplo are not consistent!"];
            Abort[];
        ];
        If[ ! (And @@ Thread[temp[[All, 2]] <= chrLength]) && 
          And @@ Thread[temp[[All, 1]] > 0],
            Print["checkfounderHaplo: SNP marker locations are out of the range!"];
            Abort[];
        ];
        True
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
    
geneDropping[popfgl_, founderHaplo_] :=
    Module[ {popfgldiplo, grid, popfgldiplogrid, nFounder, posmalefounder,fhaplo, truefgldiplo, isdiploid,dropdiplo,haplo,posmale,colX,i,j},
        checkfounderHaplo[popfgl[[1]], founderHaplo];
        isdiploid = !(And @@ (Equal @@ # & /@ popfgl[[1, 2]][[2 ;;, -1, 1, All, 1, 2]]));
        popfgldiplo = calFglDiplo[popfgl];
        grid = Transpose[founderHaplo[[{2, 3}, 2 ;;]]];
        grid = #[[All, 2]] & /@ SplitBy[grid, First];
        popfgldiplogrid = calFglDiploGrid[popfgldiplo, grid];
        posmale = Flatten[Position[popfgldiplogrid[[2 ;;,3]],2]];
        colX = 1 + Flatten[Position[founderHaplo[[2, 2 ;;]], "X" | "x"]];
        (*to check consistency of founder genders with pedigree*)
        posmalefounder = Flatten[Position[getgenders[founderHaplo[[4 ;;, colX]],isdiploid],"Male"]];
        truefgldiplo = Flatten[#, 1] & /@ popfgldiplogrid[[2 ;;, -1, All, All, 2]];
        truefgldiplo = Join[Transpose[{"Sample" <> ToString[#] & /@ Range[Length[truefgldiplo]]}], truefgldiplo, 2];
        fhaplo = founderHaplo;
        (*assuming that haplotype alleles 1 and 2 appear only in male founders*)
        fhaplo[[posmalefounder+3,colX]] = fhaplo[[posmalefounder+3,colX]]/.{"1"->"1Y","2"->"2Y","N"->"NY"};
        truefgldiplo = Join[fhaplo, truefgldiplo];
        dropdiplo = truefgldiplo;
        nFounder = Length[fhaplo] - 3;
        Do[haplo = Flatten[Characters[#] & /@ dropdiplo[[4 ;; nFounder + 3, i]]];
           haplo = Table[Extract[haplo, dropdiplo[[nFounder + 4 ;;, i, {j}]]], {j, 2}];
           dropdiplo[[nFounder + 4 ;;, i]] = StringJoin @@ # & /@ Transpose[haplo], {i, 2, Dimensions[dropdiplo][[2]]}];
        dropdiplo[[nFounder + 4 ;;,2;;]] = dropdiplo[[nFounder + 4 ;;,2;;]] /. Thread[{"1N", "2N", "N1", "N2"} -> "NN"];
        (*gender\[Equal]2===male,only X haplotype*)
        (*Take only X chromsome for males*)
        truefgldiplo[[posmalefounder+3,colX]] = Map[StringTake[#, 1] &, truefgldiplo[[posmalefounder+3,colX]], {2}];
        dropdiplo[[posmalefounder+3,colX]] = Map[StringTake[#, 1] &, dropdiplo[[posmalefounder+3,colX]], {2}];
        dropdiplo[[posmale+nFounder+3,colX]] = Map[StringTake[#, 1] &, dropdiplo[[posmale+nFounder+3,colX]], {2}];
        truefgldiplo[[posmale+nFounder+3,colX]] = truefgldiplo[[posmale+nFounder+3,colX,1]];
        {Join[{{"nFounder", nFounder}}, dropdiplo], Join[{{"nFounder", nFounder}}, truefgldiplo]}
    ]
    
applygenoerr[geno_, eps_] :=
    Module[ {rules, err, pos, obs, res = geno,ii},
          (*haplotypes "1","2","N" reference to X of male XY sex chromosomes*)
        rules = {{"1", "2", "N", "11", "12", "1N", "21", "22", "2N", "N1", 
           "N2", "NN"}, {"2", "1", "N", "21", "22", "2N", "11", "12", "1N", 
           "N1", "N2", "NN"}, {"1", "2", "N", "12", "11", "1N", "22", "21", 
           "2N", "N2", "N1", "NN"}, {"2", "1", "N", "22", "21", "2N", "12", 
           "11", "1N", "N2", "N1", "NN"}};
        (*The rules[[ii]] for ii=1,2,3,4,refer to 0,1,1,2 allelic errors respectively.*)
        rules = Thread[rules[[1]] -> #] & /@ rules;
        err = RandomChoice[{(1 - eps)^2, eps (1 - eps), (1 - eps) eps, 
            eps^2} -> Range[4], Dimensions[res]];
        Do[
         pos = Position[err, ii];
         obs = Extract[res, pos] /. rules[[ii]];
         res = ReplacePart[res, Thread[pos -> obs]], {ii, 2, 4}];
        res
    ]

applyrandommissing[geno_, missingfraction_,isGBS_] :=
    Module[ {pos, missinggeno},
        pos = Position[RandomChoice[{1 - missingfraction, missingfraction} -> {0, 1}, Take[Dimensions[geno],2]], 1];
        If[ isGBS,
            missinggeno = Table["0|0", {Length[pos]}],
            missinggeno = StringLength[Extract[geno, pos]] /. {1 -> "N", 2 -> "NN"}
        ];
        ReplacePart[geno, Thread[pos -> missinggeno]]
    ]
  
 
checkfounderSNP[founderHaplo_,isfounderinbred_,ismissing_] :=
    Module[ {ls,set},
        ls = Union[Flatten[founderHaplo[[4 ;;, 2 ;;]]]];
        set = If[ ismissing,
                  If[ isfounderinbred,
                      {"1", "2", "N"},
                      {"1", "2", "N", "11", "12", "1N", "21", "22", "2N", "N1", "N2", "NN"}
                  ],
                  If[ isfounderinbred,
                      {"1", "2"},
                      {"1", "2", "11", "12",  "21", "22"}
                  ]
              ];
        If[ !SubsetQ[set,ls],
            Print["The possible alleles of founder genotypes must be ",set,"; the alleles ", Complement[ls,set], " are not allowed!"];
            Abort[];
        ];
        True
    ]


simulatepopfgl[pedigree_, sampleInfor_,isfounderinbred_,isOogamy_,chrLength_,interferStrength_, isObligate_] :=
    Module[ {isFunnelScheme,nFounder,samplesize,popfgl,
        ind,lineid,founderFGL,inipop,iniPopFGL,pedfgl,memberlist},
        (*founders are coded as natural numbers starting from 1*)
        nFounder = Count[pedigree[[2 ;;, -1]], {0, 0}];
        memberlist = sampleInfor[[2 ;;, 2]];
        isFunnelScheme = Length[Union[memberlist]]<Length[memberlist];
        samplesize = Length[sampleInfor] - 1;
        If[ isFunnelScheme,
            Monitor[  
            popfgl = Table[
              lineid = sampleInfor[[ind + 1, 2]];
              If[ isfounderinbred,
                  founderFGL = Transpose[{sampleInfor[[ind + 1, -1]], sampleInfor[[ind + 1, -1]]}],
                  founderFGL = Transpose[{2 sampleInfor[[ind + 1, -1]] - 1, 2 sampleInfor[[ind + 1, -1]]}]
              ];
              inipop = setIniPop[nFounder, isOogamy];
              iniPopFGL = setIniPopFGL[inipop, founderFGL, chrLength];
              pedfgl = simMagicFgl[pedigree, iniPopFGL, interferStrength, isObligate, isOogamy, pedigreeMemberList -> All];
              First[Select[Rest[pedfgl], #[[2]] == lineid &, 1]], {ind, samplesize}];
              , ProgressIndicator[ind, {1, samplesize}]];
            popfgl = Join[pedfgl[[{1}]], popfgl],
            If[ isfounderinbred,
                founderFGL = Transpose[{Range[nFounder], Range[nFounder]}],
                founderFGL = Partition[Range[1, 2 nFounder], 2]
            ];
            inipop = setIniPop[nFounder, isOogamy];
            iniPopFGL = setIniPopFGL[inipop, founderFGL, chrLength];
            popfgl = simMagicFgl[pedigree, iniPopFGL, interferStrength, isObligate, isOogamy, pedigreeMemberList -> All];            
            If[ !(Complement[memberlist, popfgl[[2 ;;, 2]]] === {}),
                Print["magicGeneDropping: population size in the last generation is smaller than samplesize!"];
                Abort[]
            ];
            popfgl = popfgl[[Join[{1}, 1 + Flatten[Position[popfgl[[2 ;;, 2]], #, {1}, 1, Heads -> False] & /@ memberlist]]]];
        ];
        popfgl
    ]

          

          
checkreaddepth[founderdepth_,offspringdepth_,nFounder_,nOffspring_,nSNP_] :=
    Module[ {message,isfounderGBS,isoffspringGBS},
        message = "checkreaddepth: the option founderReadDepth must be either a positive real value or a matrix of non-negative intergers with dimensions  {nFounder,nSNP} = ";
        Switch[Dimensions[founderdepth],
                {},
                If[ ! (founderdepth >= 0),
                    Print[message, {nFounder, nSNP}, "."];
                    Abort[];
                ];
                isfounderGBS = Positive[founderdepth],
                {nFounder, nSNP},
                If[ ! (And @@ Map[IntegerQ[#] && # >= 0 &, Flatten[founderdepth]]),
                    Print[message, {nFounder, nSNP}, "."];
                    Abort[];
                ];
                isfounderGBS = True,
                _,
                Print[message,  {nFounder, nSNP}, "."];
                Abort[];
            ];
        message = "checkreaddepth: the option offspringdepth must be either a positive real value or a matrix of non-negative intergers with dimensions  {nOffspring,nSNP} = ";
        Switch[Dimensions[offspringdepth],
            {},
            If[ ! (offspringdepth >= 0),
                Print[message,  {nOffspring, nSNP}, "."];
                Abort[];
            ];
            isoffspringGBS = Positive[offspringdepth],
            {nOffspring, nSNP},
            If[ ! (And @@ Map[IntegerQ[#] && # >= 0 &, Flatten[offspringdepth]]),
                Print[message, {nOffspring, nSNP}, "."];
                Abort[];
            ];
            isoffspringGBS = True,
            _,
            Print[message,  {nOffspring, nSNP}, "."];
            Abort[];
        ];
        {isfounderGBS,isoffspringGBS}
    ]          

Options[magicGeneDropping] = {
    founderAllelicError -> 0.005,
    offspringAllelicError -> 0.005, 
    isFounderInbred -> True,
    dropMonomorphic->False,
    interferenceStrength ->0,
    isObligate ->False,
    phredQualScore ->30,
    founderGenotypeMissing -> 0,
    founderReadDepth -> 0,
    offspringGenotypeMissing ->0,
    offspringReadDepth ->0,
    depthDistribution ->"NegativeBinomial-1",
    outputFileID -> "",
    isPrintTimeElapsed ->True
}
                    
magicGeneDropping[inputfounderHaplo_?(ListQ[#]||StringQ[#]&), inputpedigreeInfor_?(VectorQ[#, ListQ]||StringQ[#]&), opts : OptionsPattern[]] :=
    Module[ {founderHaplo = inputfounderHaplo,pedigreeInfor = inputpedigreeInfor,pos,outputid,isprint,isfounderinbred,phredscore,isOogamy,chrLength,nFgl,
        founderdepth,offspringdepth, founderfraction,offspringfraction,isdropmono,starttime,pedigree, sampleInfor, nFounder, nOffspring,nSNP,founderid,popfgl,
        truediplo, truefgldiplo, obsdiplo, genderlist,truefgl,filelist,rule,sampleid,temp,isfounderGBS,isoffspringGBS,
        depthdist,epsF,eps,interferencestrength,isobligate},
        {phredscore,founderdepth,offspringdepth,founderfraction,offspringfraction,isdropmono} = OptionValue@{
         phredQualScore,founderReadDepth,offspringReadDepth,founderGenotypeMissing,offspringGenotypeMissing,dropMonomorphic};
        {depthdist,epsF,eps,isfounderinbred,interferencestrength,isobligate,outputid,isprint} = OptionValue@{
            depthDistribution,founderAllelicError,offspringAllelicError,isFounderInbred,interferenceStrength,isObligate,outputFileID,isPrintTimeElapsed};
        If[ !StringMatchQ[depthdist, "Poisson" | ("NegativeBinomial-" ~~ DigitCharacter ..)],
            Print["magicGeneDropping: unknow depthdist = ",depthdist, ". It must be \"Poisson\" or \"NegativeBinomial-n\"!, where n is a positive shape parameter!"];
            Abort[];
        ];
        If[ outputid=!="",
            outputid = outputid<>"_"
        ];
        If[ !(1>=#>=0),
            Print["magicGeneDropping: missing fraction of genotypes must be in the range between 0 and 1."];
            Abort[]
        ]&/@{founderfraction,offspringfraction};
        If[ StringQ[pedigreeInfor],
            If[ !FileExistsQ[pedigreeInfor],
                Print["magicGeneDropping: File ", pedigreeInfor," does not exist!"];
                Abort[];
            ];
            pedigreeInfor = Import[pedigreeInfor, "CSV",Path->Directory[]]
        ];
        {pedigree, sampleInfor} = Rest[splitPedigreeInfor[pedigreeInfor]];       
        nFounder = Count[pedigree[[2 ;;, 4]], {0, 0}]; 
        nOffspring = Length[sampleInfor] - 1;
        If[ StringQ[founderHaplo],
            If[ !FileExistsQ[founderHaplo],
                Print["magicGeneDropping: File ", founderHaplo," does not exist!"];
                Abort[];
            ];
            founderHaplo = Import[founderHaplo, "CSV",Path->Directory[]];
        ];
        founderHaplo[[4 ;;, 2 ;;]] = Map[ToString, founderHaplo[[4 ;;, 2 ;;]], {2}];     
        nSNP = Length[First[founderHaplo]] - 1;
        (*If there are too many founder haplotypes, take only the beginning*) 
        If[ MatchQ[Union[Flatten[founderHaplo[[4 ;;, 2 ;;]]]], {"1", "2", "N"} | {"1", "2"}],
        	(*If founderhaplo contains inbred founders or outbred founders?*)
            If[ Length[founderHaplo]-3<nFounder (1 + Boole[! isfounderinbred]),
                Print["magicpopGeneDropping: no enough input founders!"];
                Abort[];
            ];
            founderHaplo = Take[founderHaplo, 3 + nFounder (1 + Boole[! isfounderinbred])];
            If[ !isfounderinbred,
                founderHaplo = toOutbredFounder[founderHaplo];
            ],
            If[ Length[founderHaplo]-3<nFounder,
                Print["magicpopGeneDropping: no enough input founders!"];
                Abort[];
            ];
            founderHaplo = Take[founderHaplo, 3 + nFounder];
        ];   
        If[isdropmono,
		  (*drop markers that are observed to be definitely monomorphic in founders!*)
		  pos = MemberQ[{{"1"}, {"2"}, {"11"}, {"22"}}, Union[#]] & /@ Transpose[founderHaplo[[4 ;;, 2 ;;]]];
		  pos = Flatten[Position[pos, False]];
		  founderHaplo = founderHaplo[[All, Join[{1}, 1 + pos]]];
		 ];
        {isfounderGBS,isoffspringGBS} = checkreaddepth[founderdepth,offspringdepth,nFounder,nOffspring,nSNP];
        checkfounderSNP[founderHaplo,isfounderinbred,isfounderGBS&&isoffspringGBS];
        starttime = SessionTime[];
        If[ isprint,
            Print["magicGeneDropping. Start date = ",DateString[]];
        ];
        isOogamy = StringMatchQ[ToString[founderHaplo[[2, -1]]], "X" | "x"];
        (*chrLength = Ceiling[SplitBy[Transpose[founderHaplo[[2 ;; 3, 2 ;;]]], First][[All, -1, -1]]/10] 10;*)
        chrLength = Ceiling[SplitBy[Transpose[founderHaplo[[2 ;; 3, 2 ;;]]], First][[All, -1, -1]]];
        If[ isprint,
            Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \tStart simulating offspring genomes in terms of FGL lists."];
        ];
        popfgl = simulatepopfgl[pedigree, sampleInfor,isfounderinbred,isOogamy,chrLength,interferencestrength,isobligate];
        If[ isprint,
            Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \tStart dropping founder alleles on the offspring FGL lists."];
        ];
        {truediplo, truefgldiplo} = geneDropping[popfgl, founderHaplo];
        If[ isprint,
            Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \tStart simulating allelic errors in founders and offspring."];
        ];
        obsdiplo = simObservedGenotype[truediplo,phredscore,isfounderinbred,{epsF,eps},{isfounderGBS,isoffspringGBS},{founderfraction,offspringfraction},{founderdepth,offspringdepth},depthdist];
        If[ isprint,
            temp = Flatten[obsdiplo[[5 ;; nFounder + 4, 2 ;;]]];
            Print["The realized fraction of missing genotypes in founders: ", Round[Count[temp, "N" | "NN" | "0|0"]/Length[temp], 10^(-3.)]];
            temp = Flatten[obsdiplo[[nFounder + 5 ;;, 2 ;;]]];
            Print["The realized fraction of missing genotypes in offspring: ",Round[Count[temp, "N" | "NN" | "0|0"]/Length[temp], 10^(-3.)]];
        ];
        genderlist = sampleInfor[[2 ;;, 2]] /. (Rule @@ # & /@ pedigree[[2 ;;, {2, 3}]]);
        (*Put[obsdiplo,genderlist,"temp.txt"];*)
        truefgl = Transpose[{sampleInfor[[2 ;;, 1]], popfgl[[2 ;;, 1]],popfgl[[2 ;;, 2]], genderlist, sampleInfor[[2 ;;, 3]], popfgl[[2 ;;, -1]]}];
        truefgl = Join[{Join[{"ProgenyLineID"}, pedigree[[1, ;; 3]], sampleInfor[[1, {3}]], {"genome-FGL"}]}, truefgl];
        (*gender of individuals is encoded into theirs ids.*)
        If[ MemberQ[Union[founderHaplo[[2, 2 ;;]]], "X" | "x"]&&(!isfounderinbred)&&isfounderGBS,
            founderid = MapThread[StringJoin[#1, "_", #2] &, {truediplo[[5 ;; nFounder + 4, 1]],
                        Pick[pedigree[[2 ;;, 3]], pedigree[[2 ;;, 4]], {0, 0}] /. {1 -> "Female", 2 -> "Male",0 -> "hermaphrodite"}}];
            Print[founderid,"l;; here"];
            truediplo[[5 ;; nFounder + 4, 1]] = obsdiplo[[5 ;; nFounder + 4, 1]] = truefgldiplo[[5 ;; nFounder + 4, 1]] =  founderid;
        ];
        If[ MemberQ[Union[founderHaplo[[2, 2 ;;]]], "X" | "x"]&&isoffspringGBS,
            sampleid = MapThread[StringJoin[#1, "_", #2] &, {sampleInfor[[2 ;;, 1]], genderlist /. {1 -> "Female", 2 -> "Male", 0 -> "hermaphrodite"}}],
            sampleid = sampleInfor[[2 ;;, 1]]
        ];
        truediplo[[nFounder + 5 ;;, 1]] = obsdiplo[[nFounder + 5 ;;, 1]] = truefgldiplo[[nFounder + 5 ;;, 1]] =  sampleid;
        (*fgl are coded as natural numbers starting from 1*)
        nFgl = nFounder (1+Boole[!isfounderinbred]);
        rule = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        rule = MapIndexed[#1 -> First[#2] &, rule];
        truefgldiplo[[nFounder + 5 ;;, 2 ;;]] = truefgldiplo[[nFounder + 5 ;;, 2 ;;]] /. rule;
        filelist = {"TrueValues_Fgl.txt","TrueValues_FglDiplotype.csv",
                    "TrueValues_Diplotype.csv","ObservedGenotype.csv"};
        filelist = outputid<> #&/@filelist;
        If[ isprint,
            Print["Saving true values and simulated data. ",DateString[],"\n\t", Transpose[{{"True FGL", "True FGL diplotype", 
            "True diplotype", "Observed genotype"},filelist}]//TableForm];
        ];
        Put[Sequence @@ truefgl, First[filelist]];
        MapThread[csvExport[#1,#2]&,{Rest[filelist],{truefgldiplo,truediplo,obsdiplo}}];
        If[ isprint,
            Print["Done! Finished date =",DateString[], ". \tTime elapsed = ", Round[SessionTime[] - starttime,0.1], " Seconds."]
        ];
        filelist
    ]

simAlleleRead[alleles_, baseerrorprob_, readdepth_,depthdist_] :=
    Module[ {reads, pos, counts, temp,dim,n,mean},
        dim = Take[Dimensions[alleles], 2];
        If[ MatrixQ[readdepth],
            If[ Dimensions[readdepth]!=dim,
                Print["simAlleleRead: Inconsistent dimensions between alleles and readdepth!"];
                Abort[]
            ];
            reads = Map[{#,RandomInteger[BinomialDistribution[#, 1/2]]} &, readdepth, {2}];
            reads[[All,All,1]]-= reads[[All,All,2]],
            (*readdeth/2, where denotes ploidy level n=2*)
            mean = readdepth/2;
            Which[
                StringMatchQ[depthdist,"Poisson"],
                reads = RandomInteger[PoissonDistribution[mean], Append[dim, 2]],
                StringMatchQ[depthdist, ("NegativeBinomial-" ~~ DigitCharacter ..)],
                n = ToExpression[Last[StringSplit[depthdist, "-"]]];
                If[ TrueQ[Positive[n]],
                    reads = RandomInteger[NegativeBinomialDistribution[n,n/(n+mean)], Append[dim, 2]],
                    Print["simAlleleRead: shape parameter n for NegativeBinomial distribution must be positive. n= ",n];
                    Abort[];
                ];
            ];
        ];
        pos = Position[Map[Length, alleles, {2}], 1];
        reads = ReplacePart[reads, Thread[pos -> Extract[reads, pos][[All, {1}]]]];
        pos = Position[reads, _?Positive, {3}, Heads -> False];
        counts = Extract[reads, pos];
        temp = RandomInteger[BinomialDistribution[#, 1 - baseerrorprob]] & /@ counts;
        counts = Transpose[{temp, counts - temp}];
        temp = Flatten[Position[Extract[alleles, pos], "2"]];
        counts[[temp]] = counts[[temp, {2, 1}]];
        reads = ReplacePart[reads, Thread[pos -> counts]];
        reads = Replace[reads, {0 -> {0, 0}}, {3}];
        reads = Map[Total, reads, {2}];
        toDelimitedString[reads,"|"]
    ]

applyallelicerror[alleles_, errorprob_] :=
    Module[ {pos},
        pos = Replace[alleles, _ :>RandomChoice[{errorprob, 1 - errorprob} -> {1, 0}], {3}];
        pos = Position[pos, 1, {3}, Heads -> False];
        ReplacePart[alleles, Thread[pos -> (Extract[alleles, pos] /. {"1" -> "2", "2" -> "1"})]]
    ]
  
simObservedGenotype[popdiplo_,phredscore_,isfounderinbred_,errorprobs_, isGBS_,missingfractions_,depths_,depthdist_] :=
    Module[ {nFounder, baseerrorprob,alleles, res = popdiplo},
        nFounder = popdiplo[[1, 2]];
        baseerrorprob = 10^(-phredscore/10.);
        If[ isGBS[[1]],
            alleles = Map[Characters, popdiplo[[5 ;; nFounder + 4, 2 ;;]], {2}];
            alleles = applyallelicerror[alleles, errorprobs[[1]]];
            If[ isfounderinbred,
                alleles = Replace[alleles, {x_} :> {x, x}, {2}]
            ];
            res[[5 ;; nFounder + 4, 2 ;;]] = simAlleleRead[alleles, baseerrorprob, depths[[1]],depthdist],
            res[[5 ;; nFounder + 4, 2 ;;]] = applygenoerr[res[[5 ;; nFounder + 4, 2 ;;]], errorprobs[[1]]]/. Thread[{"1N", "2N", "N1", "N2"} -> "NN"]
        ];
        res[[5 ;; nFounder + 4, 2 ;;]] = applyrandommissing[res[[5 ;; nFounder + 4, 2 ;;]], missingfractions[[1]],isGBS[[1]]];
        If[ isGBS[[2]],
            alleles = Map[Characters, popdiplo[[nFounder + 5 ;;, 2 ;;]], {2}];
            alleles = applyallelicerror[alleles, errorprobs[[2]]];
            res[[nFounder + 5 ;;, 2 ;;]] = simAlleleRead[alleles, baseerrorprob, depths[[2]],depthdist],
            res[[nFounder + 5 ;;, 2 ;;]] = applygenoerr[res[[nFounder + 5 ;;, 2 ;;]], errorprobs[[2]]] /. Thread[{"1N", "2N", "N1", "N2"} -> "NN"]
        ];
        res[[nFounder + 5 ;;, 2 ;;]] = applyrandommissing[res[[nFounder + 5 ;;, 2 ;;]], missingfractions[[2]],isGBS[[2]]];
        res
    ]


toOutbredFounder[founderHaplo_] :=
    Module[ {fhaplo, genderlist, posMale, posX},
        fhaplo = Transpose[#] & /@Partition[Map[ToString, founderHaplo[[4 ;;, 2 ;;]], {2}], 2];
        fhaplo = Join[List /@ founderHaplo[[4 ;; ;; 2, 1]],Map[StringJoin @@ # &, fhaplo, {2}], 2];
        fhaplo = Join[founderHaplo[[;; 3]], fhaplo];
        genderlist = 
         If[ OddQ[#],
             "Female",
             "Male"
         ] & /@ Range[Length[fhaplo] - 3];
        posMale = Flatten[Position[genderlist, "Male"]] + 3;
        posX = Flatten[Position[founderHaplo[[2, 2 ;;]], "X" | "x"]] + 1;
        fhaplo[[posMale, posX]] = Map[StringTake[#, 1] &, fhaplo[[posMale, posX]], {2}];
        fhaplo
    ]
      
End[]

EndPackage[]


