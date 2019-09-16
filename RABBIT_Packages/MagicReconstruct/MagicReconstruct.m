(* Mathematica Package *)

(* Created by the Wolfram Workbench Mar 18, 2014 *)

BeginPackage["MagicReconstruct`",{"MagicDefinition`","ContinuousTimeHmm`","MagicDataPreprocess`",
	"MagicDataFilter`",
    "MagicReconstruct`MagicModel`","MagicReconstruct`MagicLikelihood`",
    "MagicReconstruct`MagicAlgorithm`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicReconstruct`*"];
ClearAll @@ Names["MagicReconstruct`*"];

magicReconstruct::usage = "magicReconstruct[magicsnp, model, popdesign] performs haplotype reconstruction in multi-parental populations. The magicsnp specifies the input genotypic data matrix or filename. The model specifies whether the maternally and patermally derived chromosomes are indepdent (\"indepModel\"), completely dependent (\"depModel\"), or  modeled jointly (\"jointModel\"). The popdesign speficies the breeding design information in three possible ways: a list of mating schemes from founder population to the last generation, a list of values denoting the junction distribution, or a filename for population pedigree information. "

reconstructAlgorithm::usage = "reconstructAlgorithm is an option to specify the alogrithm for haplotype reconstruction, it must be be one of \"origPathSampling\", \"origPosteriorDecoding\", or \"origViterbiDecoding\". "

outputSmooth::usage = "outputSmooth is an option to smooth the posterior probabilities when Method->origPosteriorDecoding by grouping nearby snps so that the length of snp clusters < outputSmooth (in cM) except for singleton clusters. By default, outputSmooth ->0. "

saveAsSummaryMR::usage = "saveAsSummaryMR[resultFile, summaryFile] produces a summary of the resultFile produced by magicReconstruct, and saves in summaryFile."

saveCondProb::usage = "saveCondProb[resultFile, summaryFile,probtype] extracts conditional probabilites of type probtype from the resultFile produced by magicReconstruct, and saves in summaryFile."

getSummaryMR::usage = "getSummaryMR[summaryFile] imports the summaryFile produced by saveAsSummaryMR."

getCondProb::usage = "getCondProb[summaryFile] gives conditional probabilities of L^2 diplotypes for each offspring, where L is the number of FGLs, and summaryFile is the summaryFile produced by saveAsSummaryMR."

toGenoProb::usage = "toGenoProb[diploprob] transfers conditional probabilities from L^2 diplotypes into L(L+1)/2 genotypes for each offspring, where L is the number of FGLs."

toHaploProb::usage = "toHaploProb[genoprob] transfers conditional probabilities from L(L+1)/2 genotypes into L haplotypes for each offspring, where L is the number of FGLs."

toIBDProb::usage = "toIBDProb[genoprob] transfers conditional probabilities from L(L+1)/2 genotypes into one IBD probability for each offspring, where L is the number of FGLs."

calOrigLogl::usage = "calOrigLogl[magicSNP, model, epsF, eps, popDesign] calculates the marginal likelihood. "

calOrigGeneration::usage = "calOrigGeneration[magicSNP, model, epsF, eps, Fmin, popDesignMax] calculates the posterior probabilities of the sample generation. "

calAddKinship::usage = "calAddKinship[haploprobfile, ibdprobfile, outputkinshipfile] calculates kinship matrix from the inputfile of haplotype probability and the inputfile of ibd probability. "

plotAncestryProbGUI::usage = "plotAncestryProbGUI[summaryfile, trueFGLdiplofile] or plotAncestryProbGUI[summaryfile] returns an animate for visualizing posterior probability. The summaryfile is the outputfile returned by saveAsSummaryMR.  "

isPlotGenoProb::usage = "isPlotGenoProb is an option to specify whether to plot genotype probability (True) or haplotype probability (False). "

plotCondProb::usage = "plotCondProb  "

plotAncestryProb::usage = "plotAncestryProb  "


Begin["`Private`"]
(* Implementation of the package *)
     
Options[magicReconstruct] = {
    founderAllelicError -> 0.005,
    offspringAllelicError -> 0.005, 
    isFounderInbred -> True,
    reconstructAlgorithm -> "origPathSampling",
    sampleSize -> 1000,
    sequenceDataOption ->{
        isOffspringAllelicDepth -> Automatic,
        minPhredQualScore -> 30
        },
    outputSmooth ->0,
    outputFileID -> "", 
    isPrintTimeElapsed ->True    
}
            
magicReconstruct::wrongMethod :=
    "The option reconstructAlgorithm has to take \"origPathSampling\", \"origPosteriorDecoding\", or \"origViterbiDecoding\"."

magicReconstruct::wrongModel :=
    "The 2nd model input parameter has to take \"jointModel\", \"indepModel\", or \"depModel\"."

magicReconstruct[inputmagicSNP_?(ListQ[#] ||StringQ[#]&),model_String,inputpopDesign_?(ListQ[#] ||StringQ[#]&), opts : OptionsPattern[]] :=
    Module[ {magicSNP = inputmagicSNP, popDesign = inputpopDesign,isfounderdepth = False,epsF,eps,isfounderinbred,sampleLabel,sampleMarkovProcess, 
        outputfile,outputid,algorithmname,samplesize, outputsmooth,isprint,isoffspringdepth,minphredscore,starttime,deltd, 
        nFounder,nFgl,posA,posX,foundergender,offspringgender,founderid,sampleid, snpMap,haploMap,outstream, founderHaplo, 
        derivedGeno,obsGeno, res,ind, diplotypes,haplotodiplo,printstep,weights,clusters,clusterhaploMap,delay},
        {epsF,eps} = OptionValue@{founderAllelicError,offspringAllelicError};
        {isoffspringdepth,minphredscore} = OptionValue[Thread[sequenceDataOption -> {isOffspringAllelicDepth,minPhredQualScore}]];
        {isfounderinbred,algorithmname, samplesize,outputsmooth,outputid,isprint} = OptionValue@{
            isFounderInbred,reconstructAlgorithm, sampleSize,outputSmooth,outputFileID,isPrintTimeElapsed};
        If[ outputid=!="",
            outputfile = outputid<>"_magicReconstruct.txt",
            outputfile = "magicReconstruct.txt"
        ];
        If[ isprint,
            starttime = SessionTime[];
            Print["magicReconstruct. Start date = ",DateString[], "\toutputfile = ",outputfile];
        ];
        If[ isprint&&(!isfounderinbred),
            PrintTemporary["Important: the genotypes of outbred founders are assumed to be phased! If unphased, perform parental linkage phasing using magicImpute!"];
        ];
        If[ !MemberQ[{"jointModel","indepModel","depModel"},model],
            Print["magicImputeFounder: model has to take \"jointModel\", \"indepModel\", or \"depModel\"."];
            Return[$Failed]
        ];
        If[ MemberQ[{"origPathSampling", "origPosteriorDecoding", "origViterbiDecoding"},algorithmname],
            origalgorithm = Symbol[algorithmname],
            Message[magicReconstruct::wrongMethod];
            Return[$Failed]
        ];
        If[ StringQ[magicSNP],
            If[ !FileExistsQ[magicSNP],
                Print["File ", magicSNP," does not exist!"];
                Return[$Failed]
            ];
            magicSNP = Import[magicSNP,"CSV"];
            magicSNP = DeleteCases[magicSNP, {}];
        ];
        If[ StringQ[popDesign],
            If[ !FileExistsQ[popDesign],
                Print["File ", popDesign," does not exist!"];
                Return[$Failed]
            ];
            popDesign = Import[inputpopDesign,"CSV"];
            popDesign = DeleteCases[popDesign, {}];
        ];
        {isfounderdepth,isoffspringdepth} = checkOptionValue[magicSNP, isfounderdepth, isoffspringdepth,isfounderinbred,isprint];
        SNPValidation[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        (*founderHapo dimensions {nfounder,nchr,nsnp}*)
        (*obsGeno dimensions     {noffspring,nchr,nsnp}*)
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,
            founderid,sampleid} = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        (*Put[nFounder, popDesign, isfounderinbred,model, posA, posX, offspringgender, sampleid,deltd,"tempwheat.txt"];
        Abort[];*)
        {sampleLabel, sampleMarkovProcess} = sampleDiscretePriorProcess[nFounder, popDesign, isfounderinbred,model, posA, posX, offspringgender, sampleid,deltd];
        {clusters,weights,clusterhaploMap} = getClusterhaploMap[haploMap,outputsmooth];
        Quiet[Close[outputfile]];
        Put[outputfile];
        outstream = OpenAppend[outputfile];
        Write[outstream,{model,epsF,eps,popDesign,outputfile,algorithmname, samplesize,outputsmooth,isprint,isfounderinbred,isoffspringdepth,
            foundergender,offspringgender,posA,posX,founderid,sampleid,snpMap,founderHaplo, obsGeno,clusters,weights,clusterhaploMap}];
        nFgl = nFounder (1 + Boole[! isfounderinbred]);
        diplotypes = origGenotype[nFgl][[1, 2]];
        haplotodiplo = Flatten[Position[Equal @@ # & /@ diplotypes, True]];
        printstep = 10^Max[0, IntegerLength[Length[sampleid]] - 2];
        (*derivedGeno dimensions {nchr,nsnp,ngenotype,nfounder}*)
        derivedGeno = getDerivedGenotype[founderHaplo, ToLowerCase[model]==="depmodel",True];
        delay = If[isprint,0,10^6.];
        Monitor[Do[
             If[ isprint&&Mod[ind,printstep]==0,
                 PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \tStart analyzing "<>sampleid[[ind]]<>"--"<>ToString[ind]<>"th of "<>ToString[Length[obsGeno]]];
             ];
             res = individualReconstruct[ind,model,founderHaplo, derivedGeno,obsGeno,epsF, eps,minphredscore, posA, posX, offspringgender,
                                        sampleLabel,sampleMarkovProcess,clusters,weights,haplotodiplo,diplotypes,
                                        isoffspringdepth,algorithmname,outputsmooth,samplesize];
             Write[outstream,res], {ind, Length[obsGeno]}];
             0, ProgressIndicator[ind, {1, Length[obsGeno]}],delay];
        Close[outstream];
        If[ isprint,
            Print["Done! Finished date =",DateString[], ". \tTime elapsed = ", Round[SessionTime[] - starttime,0.1], " Seconds."];
        ];
        outputfile
    ]
         
individualReconstruct[ind_,model_,founderHaplo_, derivedGeno_,obsGeno_,epsF_, eps_, minphredscore_,posA_, posX_, offspringgender_,
    sampleLabel_,sampleMarkovProcess_,clusters_,weights_,haplotodiplo_,diplotypes_,
    isoffspringdepth_,algorithmname_,outputsmooth_,samplesize_] :=
    Module[ {haplocode,diplocode,startProb,tranProb,dataProb,res,temp,prob,i,j},
    (*haplocode is a permutation of founder genome labels, natural integers starting from 1.*)
             (*haplocode refers to reordering of haploid states; diplocode refers to reordering of diploid states*)
             (*ind+1, +1 refers to the column titles ={sampleid, priorlabel(gender/memberID), funnelcode,haplocode,diplocode}*)
        {haplocode,diplocode} = sampleLabel[[ind+1,4;;5]];
        {startProb,tranProb} = sampleLabel[[ind+1,2]]/.sampleMarkovProcess;
        If[ !OrderedQ[haplocode],
            {startProb,tranProb} = relabelMarkovProcess[{startProb,tranProb}, haplocode,diplocode];
        ];
        If[ isoffspringdepth,
            dataProb = lineMagicLikelihoodGBS[model,founderHaplo, derivedGeno,obsGeno[[ind]], epsF, eps,minphredscore, posA, posX, offspringgender[[ind]]],
            dataProb = lineMagicLikelihood[model,founderHaplo, obsGeno[[ind]], epsF, eps, posA, posX, offspringgender[[ind]]];
        ];
        res = If[ algorithmname === "origPathSampling",
                  origalgorithm[startProb, tranProb, dataProb, samplesize],
                  origalgorithm[startProb, tranProb, dataProb]
              ];
        res = Normal[res];
        If[ algorithmname === "origPosteriorDecoding"&&outputsmooth>0,
            temp = res[[All, 2]];
            res[[All, 2]] = Table[Total[temp[[i, clusters[[i, j]]]] weights[[i, j]]], {i, Length[clusters]}, {j,Length[clusters[[i]]]}];
            res[[All, 2]] = N[Round[res[[All, 2]], 10^(-5)]];
        ];
        If[ offspringgender[[ind]] == "Male"&&(ToLowerCase[model]=!="depmodel"),
            Switch[
             algorithmname,
             "origPosteriorDecoding",
             prob = ConstantArray[0, ReplacePart[Dimensions[res[[posX, 2]]], -1 -> Length[diplotypes]]];
             prob[[All, All, haplotodiplo]] = Normal[res[[posX, 2]]];
             res[[posX, 2]] = prob,
             "origViterbiDecoding",
             res[[posX, 2, ;; -2, 2]] = haplotodiplo[[#]] & /@ res[[posX, 2, ;; -2, 2]],
             "origPathSampling",
             res[[posX, 2, All, ;; -2, 2]] = Map[haplotodiplo[[#]] &, res[[posX, 2, All, ;; -2, 2]], {2}],
             _,
             Print["Wrong reconstructAlgorithm " <> algorithmname <> "!"];
            ]
        ];
        res
    ]    
      
calOrigLogl2[inputmagicSNP_List, model_String, epsF_?NonNegative, eps_?NonNegative, inputpopDesign_List] :=
    Module[ {magicSNP = inputmagicSNP, popDesign = inputpopDesign,isfounderinbred = True,isfounderdepth = False,isoffspringdepth = False,
        deltd, founderHaplo, obsGeno,foundergender,offspringgender,posA,posX,
        nFounder,founderid,sampleid,snpMap,haploMap,logl,sampleLabel, sampleMarkovProcess,startProb,tranProb,haplocode,diplocode,dataProb,ind},
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,
            founderid,sampleid} = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        {sampleLabel, sampleMarkovProcess} = sampleDiscretePriorProcess[nFounder,popDesign, isfounderinbred, model, posA, posX, offspringgender, sampleid, deltd];
        logl = Table[
          (*haplocode is a permutation of founder genome labels, natural integers starting from 1.*)
          (*haplocode refers to reordering of haploid states; diplocode refers to reordering of diploid states*)
          (*ind+1, +1 refers to the column titles ={sampleid, priorlabel(gender/memberID), funnelcode,haplocode,diplocode}*)          
          {haplocode,diplocode} = sampleLabel[[ind+1,4;;5]];
          {startProb,tranProb} = sampleLabel[[ind+1,2]]/.sampleMarkovProcess;
          If[ !OrderedQ[haplocode],
              {startProb,tranProb} = relabelMarkovProcess[{startProb,tranProb}, haplocode,diplocode];
          ];
          dataProb = lineMagicLikelihood[model,founderHaplo, obsGeno[[ind]], epsF, eps, posA, posX, offspringgender[[ind]]];
          origLogLiklihood[startProb, tranProb, dataProb], {ind,Length[obsGeno]}];
        logl
    ]

calOrigLogl[inputmagicSNP_List, model_String, epsF_?NonNegative, eps_?NonNegative, inputpopDesign_List] :=
    Module[ {magicSNP = inputmagicSNP, popDesign = inputpopDesign,isfounderinbred = True,isfounderdepth = False,isoffspringdepth = False,
        deltd, founderHaplo, obsGeno,foundergender,offspringgender,posA,posX,continuedMarkovProcess,markov,chrposA,chrposX,chr,
        nFounder,founderid,sampleid,snpMap,haploMap,sampleLabel,startProb,tranProb,haplocode,diplocode,dataProb,ind},
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,
            founderid,sampleid} = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        {sampleLabel, continuedMarkovProcess} = sampleContinuedPriorProcess[nFounder, popDesign, isfounderinbred, 
                model, posA, posX, offspringgender, sampleid];
        Sum[
            markov = continuedMarkovProcess;
            markov[[All, 2]] = markov[[All, 2, All, {chr}]];
            markov = toDiscreteMarkovProcess[markov, deltd[[{chr}]]];
            chrposA = If[ MemberQ[posA, chr],
                          {1},
                          {}
                      ];
            chrposX = If[ MemberQ[posX, chr],
                          {1},
                          {}
                      ];
            Table[
             {haplocode, diplocode} = sampleLabel[[ind + 1, 4 ;; 5]];
             {startProb, tranProb} = sampleLabel[[ind + 1, 2]] /. markov;
             If[ !OrderedQ[haplocode],
                 {startProb,tranProb} = relabelMarkovProcess[{startProb,tranProb}, haplocode,diplocode];
             ];
             dataProb = lineMagicLikelihood[model, founderHaplo[[All, {chr}]],obsGeno[[ind, {chr}]], epsF, eps, chrposA, chrposX, offspringgender[[ind]]];
             origLogLiklihood[startProb, tranProb, dataProb], {ind, Length[obsGeno]}], {chr, Length[deltd]}]
    ]
    
calOrigGeneration[magicSNP_List, model_String, epsF_?NonNegative, 
  eps_?NonNegative, Fmin_Integer, popDesignMax_?(VectorQ[#, StringQ] &)] :=
    Module[ {Fmax, schemelist,i,lls},
        Fmax = Length[popDesignMax];
        schemelist = Table[Take[popDesignMax, i], {i, Fmin, Fmax}];
        lls = calOrigGeneration[magicSNP, model, epsF, eps, schemelist];
        lls[[All, 2, 1]] = "Generation";
        lls[[All, 3 ;;, 1]] += Fmin - 1;
        lls
    ]

calOrigGeneration[magicSNP_List, model_String, epsF_?NonNegative, 
  eps_?NonNegative, schemelist_?(VectorQ[#, ListQ] &),isSameGeneration_:False] :=
    Module[ {loglls, lls, nFounder, id, i},
        Monitor[
              loglls = Table[calOrigLogl[magicSNP, model, epsF, eps, schemelist[[i]]], {i, Length[schemelist]}], 
                 ProgressIndicator[i, {1, Length[schemelist]}]];
        nFounder = magicSNP[[1, 2]];
        id = magicSNP[[nFounder + 5 ;;, 1]];
        If[ isSameGeneration,
            lls = Total[loglls, {2}];
            lls = Normalize[Exp[lls - Max[lls]], Total];
            lls = Sort[
              Transpose[{Range[Length[schemelist]], lls}], #1[[2]] > #2[[2]] &];
            lls = Join[{{"SchemeIndex", "PosteriorProb"}}, lls],
            lls = Normalize[Exp[# - Max[#]], Total] & /@ Transpose[loglls];
            lls = Sort[Transpose[{Range[Length[schemelist]], #}], #1[[2]] > #2[[2]] &] & /@ lls;
            lls = Transpose[Join[{id, Table[{"SchemeIndex", "PosteriorProb"}, {Length[id]}]}, Transpose[lls]]]
        ]
    ]

toIBDProb[genoprob_] :=
    Module[ {depth, len, nfgl, genotypes, pos,inbred,x,vec},
        depth = Depth[genoprob];
        len = Union[Flatten[Map[Last[Dimensions[#]] &, genoprob, {depth - 3}]]];
        nfgl = Solve[x (x + 1)/2 == First[len] && x > 0, x][[1, 1, 2]];
        If[ ! (Length[len] == 1 && IntegerQ[nfgl]),
            Print["toIBDProb: the argument is not an array of unphased genotype probabilities!"];
            Abort[]
        ];
        genotypes = origGenotype[nfgl][[1, 1]];
        pos = Flatten[Position[genotypes, {x_, x_}]];
        vec = SparseArray[{}, {Length[genotypes]}];
        vec[[pos]] = 1;
        inbred = Map[vec.Transpose[#] &, genoprob, {depth - 3}];
        N[Round[inbred, 10^(-5)]]
    ]  

toHaploProb[genoprob_] :=
    Module[ {depth, len, nfgl, genotypes, mtx,haploprob,x,i,j},
        depth = Depth[genoprob];
        len = Union[Flatten[Map[Last[Dimensions[#]] &, genoprob, {depth - 3}]]];
        nfgl = Solve[x (x + 1)/2 == First[len] && x > 0, x][[1, 1, 2]];
        If[ ! (Length[len] == 1 && IntegerQ[nfgl]),
            Print["toGenoProb: the argument is not an array of unphased genotype probabilities!"];
            Abort[]
        ];
        genotypes = origGenotype[nfgl][[1, 1]];
        mtx = ConstantArray[0, {Length[genotypes], nfgl}];
        Do[mtx[[i, genotypes[[i, j]]]] += 1, {i, Length[genotypes]}, {j,Dimensions[genotypes][[2]]}];
        mtx = Transpose[mtx]/Dimensions[genotypes][[2]];
        haploprob = Map[Transpose[mtx.Transpose[#]] &, genoprob, {depth - 3}];
        N[Round[haploprob,10^(-5)]]
    ]
  
(*diploprob[[i,j]]: a list of diplotype probabilities at marker j of linkage group i, for a given indvidual*)
(*also works if there are no level of linkage groups, also works if there are additonal level of individuals*)    
toGenoProb[diploprob_] :=
    Module[ {depth,len,nfgl, geno2diplo, i, mtx},
        depth = Depth[diploprob];
        len = Union[Flatten[Map[Last[Dimensions[#]] &, diploprob, {depth - 3}]]];
        nfgl = Sqrt[First[len]];
        If[ ! (Length[len] == 1 && IntegerQ[nfgl]),
            Print["toGenoProb: the argument is not an array of diplotype probabilities!"];
            Abort[]
        ];
        geno2diplo = origGenotype[nfgl][[2, 1]];
        mtx = SparseArray[{}, {nfgl (nfgl + 1)/2, nfgl^2}];
        Do[mtx[[i, geno2diplo[[i]]]] = 1, {i, Length[mtx]}];
        Map[Transpose[mtx.Transpose[#]] &, diploprob, {depth - 3}]
    ]

(*return diploprob with dimensions: {noffspring, nlinkagegroup,nsnp,nstate=nfgl^2, nfgl (depmodel/maleX)}*)
getCondProb[resultFile_String?FileExistsQ] :=
    Module[ {res},
        res = ReadList[resultFile];
        res[[2 ;;, All, 2]]
    ]

saveAsSummaryMR[resultFile_String?FileExistsQ, summaryFile_String] :=
    Module[ {res, model,epsF, eps,fglid, popDesign,outputid,algorithmname, samplesize,outputsmooth,isprint,isfounderinbred,isoffspringdepth,foundergender,offspringgender,posA,posX,
        founderid,sampleid,snpMap,founderHaplo, obsGeno,clusters,weights,clusterhaploMap,clustersnpMap, nfgl,chrid,logl,logl2,genotypes, diplotypes, genoID, genotypes2, haploID, haplotypes2, 
      rowID, genoprob2, haploprob2,ibdprob2, diploID, diplotypes2, path, summary, diplotohaplo,key = "magicReconstruct-Summary",outfile},
        PrintTemporary["Time used in reading "<>resultFile<>": "
            <>ToString[Round[AbsoluteTiming[
         res = ReadList[resultFile];][[1]],0.01]]<>" seconds. "<>DateString[]];
        {model,epsF,eps,popDesign,outputid,algorithmname, samplesize,outputsmooth,isprint,isfounderinbred,isoffspringdepth,
            foundergender,offspringgender,posA,posX,founderid,sampleid,snpMap,founderHaplo, obsGeno,clusters,weights,clusterhaploMap} = res[[1]];
        clustersnpMap = clusterhaploMap[[;;3]];        
        fglid = If[isfounderinbred, founderid, Flatten[{# <> "_m", # <> "_p"} & /@ founderid]];
        nfgl = Length[fglid];
        {genotypes, diplotypes} = origGenotype[nfgl][[1]];
        diplotohaplo = Replace[diplotypes, {{x_, x_} :> x, {_, _} -> 0}, {1}];
        haploID = "haplotype" <> ToString[#] & /@ Range[nfgl];
        haplotypes2 = Transpose[{haploID, Range[nfgl], fglid}];
        haplotypes2 = Join[{{"Haplotype", "Code", "founder"}}, haplotypes2];
        genoID = "genotype" <> ToString[#] & /@ Range[Length[genotypes]];
        genotypes2 = Transpose[{genoID,
           toDelimitedString[genotypes,"|"],
           toDelimitedString[Map[fglid[[#]] &, genotypes, {2}],"|"]}];
        genotypes2 = Join[{{"Genotype", "Code", "founder"}}, genotypes2];
        diploID = "diplotype" <> ToString[#] & /@ Range[Length[diplotypes]];
        diplotypes2 = Transpose[{diploID,
           toDelimitedString[diplotypes,"|"],
           toDelimitedString[Map[fglid[[#]] &, diplotypes, {2}],"|"]}];
        diplotypes2 = Join[{{"Diplotype", "Code", "founder"}}, diplotypes2];
        chrid = Split[clusterhaploMap[[2, 2 ;;]]][[All, 1]];
        logl = res[[2 ;;, All, 1]];
        logl2 = Join[Transpose[{Prepend[sampleid, "Lines/Chr"]}],Prepend[logl, chrid], 2];
        Switch[
             algorithmname, 
             "origPosteriorDecoding",
             If[ MatchQ[ToLowerCase[model], "depmodel"],
                 haploprob2 = res[[2 ;;, All, 2]];
                 rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, haploID]];
                 haploprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ haploprob2, 1];
                 haploprob2 = Join[clustersnpMap,Join[Transpose[{rowID}], haploprob2, 2]];
                 ibdprob2 = genoprob2 = {None};,
                 PrintTemporary["Time used in transforming from diplotprob to genoprob: "
                     <>ToString[Round[AbsoluteTiming[genoprob2 = toGenoProb[res[[2 ;;, All, 2]]];][[1]],0.01]]<>" seconds. "<>DateString[]];
                 ClearAll[res];
                 PrintTemporary["Time used in transforming from genoprob to haploprob: "
                     <>ToString[Round[AbsoluteTiming[haploprob2 = toHaploProb[genoprob2];][[1]],0.01]]<>" seconds. "<>DateString[]];
                 PrintTemporary["Time used in transforming from haploprob to ibdprob: "
                     <>ToString[Round[AbsoluteTiming[ibdprob2 = toIBDProb[genoprob2];][[1]],0.01]]<>" seconds. "<>DateString[]];
                 (*to add row/col labels for genoprob2,haploprob2,and ibdprob2*)
                 rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, genoID]];
                 genoprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ genoprob2, 1];
                 genoprob2 = Join[clustersnpMap,Join[Transpose[{rowID}], genoprob2, 2]];
                 rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, haploID]];
                 haploprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ haploprob2, 1];
                 haploprob2 = Join[clustersnpMap,Join[Transpose[{rowID}], haploprob2, 2]];
                 rowID = #1 <> "_IBDProbability" & /@ sampleid;
                 ibdprob2[[Flatten[Position[offspringgender, "Male"]], posX]] *= 0;
                 ibdprob2 = Flatten[#] & /@ ibdprob2;
                 ibdprob2 = Join[clustersnpMap, Join[Transpose[{rowID}], ibdprob2, 2]];
             ];
             (*export*)
             summary = Join[{{key, "founderhaplo","Genetic map and founder hapolotypes"}}, clusterhaploMap,
               {{key, "logl","Ln marginal likelihood"}}, logl2,
               {{key, "genotype","Genotypes in order"}}, genotypes2,
               {{key, "genoprob","Conditonal genotype probability"}}, genoprob2,
               {{key, "haplotype","haplotypes in order"}}, haplotypes2,
               {{key, "haploprob","Conditonal haplotype probability"}}, haploprob2,
               {{key, "ibdprob","Conditonal IBD probability"}}, ibdprob2];
             (* Export[summaryFile, ExportString[summary, "CSV"], "Table"],*)
             (*MapThread[PrintTemporary["Time used in exporting "<>#1<>":"
                 <>ToString[Round[AbsoluteTiming[csvExport[#1,#2];][[1]],0.01]]
                 <>" seconds. "<>DateString[]]&, 
                     {StringDrop[summaryFile, -4] <> # <> ".csv" & /@ {"_founderhaplo","_genoprob","_haploprob", "_ibdprob"},
                     {clusterhaploMap,genoprob2,haploprob2,ibdprob2}}];*)
             PrintTemporary["Time used in exporting "<>summaryFile<>":"
                 <>ToString[Round[AbsoluteTiming[outfile = csvExport[summaryFile,summary];][[1]],0.01]]
                 <>" seconds. "<>DateString[]],
             "origViterbiDecoding",
             path = res[[2 ;;, All, 2]];
             ClearAll[res];
             path = Map[CtStringFormat[#, "-"] &, path, {2}];
             path = StringJoin[Riffle[#, "||"]] & /@ path;
             path = Join[{{"Lines", "ViterbiPath"}}, Transpose[{sampleid, path}]];
             (*export*)
             summary = Join[{{key,"founderhaplo","Genetic map and founder hapolotypes"}}, clusterhaploMap, 
                 {{key,"logl", "Ln marginal likelihood"}}, logl2,
                 {{key, "haplotype","haplotypes in order"}}, haplotypes2,
                 {{key, "diplotype","diplotypes in order"}}, If[ model === "depModel",
                                                                 {None},
                                                                 diplotypes2
                                                             ],
               {{key, "viterbipath","Viterbi path of "<>If[ model === "depModel",
                                                            "haplotypes",
                                                            "diplotypes"
                                                        ]}}, path];
             (* Export[summaryFile, ExportString[summary, "CSV"], "Table"],*)
             outfile = csvExport[summaryFile,summary],
             "origPathSampling",
             path = res[[2 ;;, All, 2]];
             ClearAll[res];
             path = Map[CtStringFormat[#, "-"] &, path, {3}];
             path = Transpose[#] & /@ path;
             path = Map[StringJoin[Riffle[#, "|"]] &, path, {2}];
             path = Map[StringJoin[Riffle[#, "||"]] &, path, {1}];
             path = Join[{{"Lines", "SampledPaths"}}, Transpose[{sampleid, path}]];
             (*export*)
             summary = Join[{{key,"founderhaplo","Genetic map and founder hapolotypes"}}, clusterhaploMap, 
               {{key,"logl", "Ln marginal likelihood"}}, logl2,
               {{key, "haplotype","haplotypes in order"}}, haplotypes2,
               {{key, "diplotype","diplotypes in order"}}, If[ model === "depModel",
                                                               {None},
                                                               diplotypes2
                                                           ],
               {{key, "viterbipath","Sampled paths of "<>If[ model === "depModel",
                                                             "haplotypes",
                                                             "diplotypes"
                                                         ]}}, path];
             (* Export[summaryFile, ExportString[summary, "CSV"], "Table"],*)
             outfile = csvExport[summaryFile,summary],
             _,
             Print["Wrong " <> resultFile <> "!"];
             outfile = Abort[]
         ];
        outfile
    ]    
    
saveCondProb[resultFile_String?FileExistsQ, summaryFile_String, probtype_String:"genotype"] :=
    Module[ {res,model,epsF,eps,fglid,popDesign,outputid,algorithmname, samplesize,outputsmooth,isprint,isfounderinbred,isoffspringdepth,
            foundergender,offspringgender,posA,posX,founderid,sampleid,snpMap,founderHaplo, obsGeno,clusters,weights,clusterhaploMap,
            clustersnpMap,nfgl,genotypes, diplotypes,diplotohaplo,haploID,haplotypes2,genoID,genotypes2,diploID,diplotypes2,chrid,
            logl,logl2,condprob,diploprob2,genoprob2,haploprob2,rowID,summary,key = "magicReconstruct-Summary",outfile},
        If[ ! MatchQ[probtype, "diplotype" | "genotype" | "haplotype"],
            Print["The probtype must be diplotype, genotype, or haplotype! wrong probtype of ", probtype, "!"];
        ];
        PrintTemporary["Time used in reading "<>resultFile<>": "
            <>ToString[Round[AbsoluteTiming[
         res = ReadList[resultFile];][[1]],0.01]]<>" seconds. "<>DateString[]];
        {model,epsF,eps,popDesign,outputid,algorithmname, samplesize,outputsmooth,isprint,isfounderinbred,isoffspringdepth,
            foundergender,offspringgender,posA,posX,founderid,sampleid,snpMap,founderHaplo, obsGeno,clusters,weights,clusterhaploMap} = res[[1]];
        If[ ! MatchQ[algorithmname, "origPosteriorDecoding"],
            Print["Conditional probabilities can be extracted only for the case of origPosteriorDecoding, but not for ", algorithmname, "!"];
            Abort[]
        ];
        If[ MatchQ[ToLowerCase[model], 
            "depmodel"] && (! MatchQ[probtype, "haplotype"]),
            Print["Conditional ", probtype, " probabilities are not available for ", model];
            Abort[];
        ];
        clustersnpMap = clusterhaploMap[[;;3]];        
        fglid = If[isfounderinbred, founderid, Flatten[{# <> "_m", # <> "_p"} & /@ founderid]];
        nfgl = Length[fglid];
        {genotypes, diplotypes} = origGenotype[nfgl][[1]];
        diplotohaplo = Replace[diplotypes, {{x_, x_} :> x, {_, _} -> 0}, {1}];
        haploID = "haplotype" <> ToString[#] & /@ Range[nfgl];
        haplotypes2 = Transpose[{haploID, Range[nfgl], fglid}];
        haplotypes2 = Join[{{"Haplotype", "Code", "founder"}}, haplotypes2];
        genoID = "genotype" <> ToString[#] & /@ Range[Length[genotypes]];
        genotypes2 = Transpose[{genoID,
           toDelimitedString[genotypes,"|"],
           toDelimitedString[Map[fglid[[#]] &, genotypes, {2}],"|"]}];
        genotypes2 = Join[{{"Genotype", "Code", "founder"}}, genotypes2];
        diploID = "diplotype" <> ToString[#] & /@ Range[Length[diplotypes]];
        diplotypes2 = Transpose[{diploID,
           toDelimitedString[diplotypes,"|"],
           toDelimitedString[Map[fglid[[#]] &, diplotypes, {2}],"|"]}];
        diplotypes2 = Join[{{"Diplotype", "Code", "founder"}}, diplotypes2];
        chrid = Split[clusterhaploMap[[2, 2 ;;]]][[All, 1]];
        logl = res[[2 ;;, All, 1]];
        logl2 = Join[Transpose[{Prepend[sampleid, "Lines/Chr"]}],Prepend[logl, chrid], 2];
        Switch[probtype,
         "diplotype",
         diploprob2 = res[[2 ;;, All, 2]];
         ClearAll[res];
         rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, diploID]];
         diploprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ diploprob2, 1];
         diploprob2 = Join[clustersnpMap, Join[Transpose[{rowID}], diploprob2, 2]];
         summary = Join[{{key,"Genetic map and founder hapolotypes"}}, clusterhaploMap, 
             {{key, "Ln marginal likelihood"}}, logl2,
             {{key, "Diplotypes in order"}}, diplotypes2, 
             {{key, "Conditonal diplotype probability"}},diploprob2],
         "genotype" | "haplotype",
         condprob = res[[2 ;;, All, 2]];
         ClearAll[res];
         If[ MatchQ[probtype, "haplotype"],
             rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, haploID]];
             haploprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ condprob, 1];
             ClearAll[condprob];
             haploprob2 = Join[clustersnpMap, Join[Transpose[{rowID}], haploprob2, 2]];
             summary = Join[{{key,"founderhaplo","Genetic map and founder hapolotypes"}}, clusterhaploMap, 
                 {{key, "logl","Ln marginal likelihood"}}, logl2,
                 {{key, "halotype", "haplotypes in order"}},haplotypes2, 
                 {{key, "haploprob","Conditonal haplotype probability"}},haploprob2],
             (*genotypes*)
             PrintTemporary["Time used in transforming from diplotprob to genoprob: "
                     <>ToString[Round[AbsoluteTiming[genoprob2 = toGenoProb[condprob];][[1]],0.01]]<>" seconds. "<>DateString[]];
             ClearAll[condprob];
             rowID = Flatten[Outer[StringJoin[#1, "_", #2] &, sampleid, genoID]];
             genoprob2 = Flatten[Transpose[Flatten[#, 1]] & /@ genoprob2, 1];
             genoprob2 = Join[clustersnpMap, Join[Transpose[{rowID}], genoprob2, 2]];
             summary = Join[{{key,"founderhaplo","Genetic map and founder hapolotypes"}}, clusterhaploMap, 
                 {{key, "logl","Ln marginal likelihood"}}, logl2,
                 {{key, "genotype","Genotypes in order"}}, genotypes2, 
                 {{key, "genoprob", "Conditonal genotype probability"}},genoprob2]
         ];
         ];
        PrintTemporary[
         "Time used in exporting " <> summaryFile <> ":" <> 
          ToString[Round[AbsoluteTiming[
                If[ Length[First[clusterhaploMap]]<5000,
                    outfile = csvExport[summaryFile, summary],
                    outfile = csvExport2[summaryFile, summary]
                ];
              ][[1]], 0.01]] <>" seconds. " <> DateString[]];
        outfile
    ]    

getSummaryMR[summaryFile_String?FileExistsQ] :=
    Module[ {res, description,key = "magicReconstruct-Summary"},
        res = readTable[summaryFile, "CSV"];
        res = Partition[Split[res, #1[[1]] != key && #2[[1]] != key &], 2];
        description = Flatten[#] & /@ res[[All, 1]];
        res = Join[{description}, res[[All, 2]]];
        res
    ]
    
calAddKinship[haploprobfile_String?FileExistsQ, ibdprobfile_String?FileExistsQ, outputkinshipfile_String] :=
    Module[ {starttime,haploprob, ibdprob, noffspring, nfgl, kinship, offspringid, i, j},
        starttime = SessionTime[];
        PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. Reading " <> haploprobfile];
        haploprob = readTable[haploprobfile,"CSV"];
        PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. Reading " <> ibdprobfile];
        ibdprob = readTable[ibdprobfile,"CSV"];
        offspringid = StringJoin @@ # & /@ (StringSplit[ibdprob[[4 ;;, 1]],"_"][[All, ;; -2]]);
        noffspring = Length[ibdprob] - 3;
        nfgl = (Length[haploprob] - 3)/noffspring;
        (*weights = Normalize[Flatten[Differences[#] & /@ SplitBy[Transpose[haploprob[[2 ;; 3, 2 ;;]]], First][[All, All,2]]], Total];*)
        (*weights = Last[Dimensions[haploprob]] - 1 - Length[Split[haploprob[[2, 2 ;;]]]];*)
        (*weights = Table[1/weights, {weights}];*)
        PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. Transforming ibdprob..."];
        (*The resulting ibdprob dimensions: {nchr,noffspring,nsnp}*)
        ibdprob = Transpose[#] & /@ (SplitBy[Transpose[ibdprob[[All, 2 ;;]]], #[[2]] &][[All, All, 4 ;;]]);
        PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. Transforming haploprob..."];
        (*The resulting haploprob dimensions: {nchr,noffspring,nsnp,nfgl}*)
        haploprob = Map[Transpose, Partition[Transpose[#], nfgl] & /@ (SplitBy[Transpose[haploprob[[All, 2 ;;]]], #[[2]] &][[All, All,4 ;;]]), {2}];
        PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. Calculating kinship..."];
        kinship = ConstantArray[0, {noffspring, noffspring}];
        Monitor[
         Do[kinship[[i, j]] = kinship[[j, i]] = 
            Mean[Flatten[Total[haploprob[[All, i]] haploprob[[All, j]], {3}][[All, ;; -2]]]], {i, noffspring}, {j,i + 1, noffspring}];
         Do[kinship[[i, i]] = (Mean[Flatten[ibdprob[[All, i, ;; -2]]]] + 1)/2, {i,noffspring}];
         0, ProgressIndicator[i, {1, noffspring}]];
        kinship *= 2;
        kinship = Join[Transpose[{Prepend[offspringid, "IndividualID"]}],Join[{offspringid}, kinship], 2];
        csvExport[outputkinshipfile,kinship];
    ]  

Options[plotAncestryProb] = Join[{isPlotGenoProb -> Automatic,linkageGroupSet->All},Options[ListPlot]]
Options[plotAncestryProbGUI] = Join[Options[plotAncestryProb],Options[ListAnimate]];

plotAncestryProbGUI[summaryfile_String?FileExistsQ, opts : OptionsPattern[]] :=
    plotAncestryProbGUI[summaryfile,None,opts]

getchrsetdata[data_, chrset_] :=
    If[ data === None,
        None,
        Rest[getsubMagicSNP[Join[{{None, None}}, data], chrset, All]]
    ]

plothaploprob[condprob_, truefgldiplo_,offspringls_,opts_] :=
    Module[ {diplorule, nfounder, offprob, truegeno, gg, label,line,ls,g1, g2,g3},        
        If[ truefgldiplo =!= None,
        	nfounder = truefgldiplo[[1, 2]];
            diplorule = Flatten[Outer[List, Range[nfounder], Range[nfounder]], 1];
            diplorule = Thread[Range[Length[diplorule]] -> diplorule];
            truegeno = truefgldiplo[[nfounder + 5 ;;, 2 ;;]] /. diplorule;
        ];
        gg = Table[
          offprob = condprob[[line]];
          label = "Offspring " <> ToString[line]<>": "<>ToString[offspringls[[line]]] <>". GrayLevel=Probability, ";
              If[ truefgldiplo =!= None,
                  label = label <>"\!\(\*
					StyleBox[\"red\",\nFontColor->RGBColor[1, 0, 0]]\)\!\(\*
					StyleBox[\" \",\nFontColor->RGBColor[1, 0, 0]]\)\!\(\*
					StyleBox[\"X\",\nFontColor->RGBColor[1, 0, 0]]\) =true IBD genotype, "<>"\!\(\*
					StyleBox[\"blue\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
					StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
					StyleBox[\"O\",\nFontColor->RGBColor[0, 0, 1]]\)\!\(\*
					StyleBox[\" \",\nFontColor->RGBColor[0, 0, 1]]\)=true non-IBD genotype.",
                  label = label <> "vertical line=Chr boundary.";
              ];
          g1 = MatrixPlot[Reverse[offprob], FrameLabel -> {"Haplotype", "SNP index"}, 
            ColorFunctionScaling -> False, MaxPlotPoints -> Infinity, 
            AspectRatio -> 1/2.5, DataRange -> {{1, Dimensions[offprob][[2]]}, {1, Length[offprob]}}];
          g2 = If[ truefgldiplo =!= None,
                   {ls = truegeno[[line]] /. {i_, i_} :> i;
                    ls = Transpose[{Range[Length[ls]], ls}];
                    ListPlot[Select[ls, MatchQ[#, {_, _Integer}] &],PlotRange -> {{1, Dimensions[offprob][[2]]}, {1,Length[offprob]}}, 
                        PlotRangePadding -> 0,PlotMarkers -> {"\[Times]", 10}, PlotStyle -> Red],
                    ls = Select[ls, MatchQ[#, {_, _List}] &];
                    ls = Flatten[Thread[#] & /@ ls, 1];
                    ListPlot[ls, PlotRange -> {{1, Dimensions[offprob][[2]]}, {1,Length[offprob]}}, PlotRangePadding -> 0, 
                     PlotMarkers -> {"\[EmptyCircle]", 14}, PlotStyle -> Blue]},
                   ListPlot[{1}, PlotMarkers -> {Automatic, 0}]
               ];
          If[ truefgldiplo =!= None,
          	g3=ListLinePlot[Thread[{#, {0.3, Length[offprob] + 0.7}}] & /@Accumulate[Tally[truefgldiplo[[3, 2 ;;]]][[All, 2]]],
                PlotStyle -> Directive[Dashed, Thick, Brown]];
          	gg={g1,g2,g3},
          	gg={g1,g2}
          ];
          Show[Sequence@@gg,opts,ImageSize -> 1000, Frame -> True, PlotLabel -> label, 
           LabelStyle -> Directive[FontSize -> 14, Black, FontFamily -> "Helvetica"]], {line, Length[condprob]}];
        gg
    ]
    
plotgenoprob[condprob_, truefgldiplo_,offspringls_,opts_] :=
    Module[ {nfounder, x,types, truegeno, gg, offprob,label, line,g1, g2, g3},
    	nfounder = x /. Flatten[Solve[x (x + 1)/2 == Dimensions[condprob][[2]] && x > 0, x]];
    	types = origGenotype[nfounder];
        If[ truefgldiplo =!= None,
        	If[nfounder == truefgldiplo[[1, 2]],Print["plotgenoprob: inconsistent nfounder"]];        	
            truegeno = truefgldiplo[[nfounder + 5 ;;, 2 ;;]] /.Thread[Range[Length[types[[2, 2]]]] -> types[[2, 2]]];
        ];
        gg = Table[
          offprob = condprob[[line]];
          label = "Offspring " <> ToString[line]<>": "<>ToString[offspringls[[line]]] <>". GrayLevel=Probability, ";
          If[ truefgldiplo =!= None,
              label = label <> "red X =true genotype.",
              label = label <> "vertical line=Chr boundary.";
          ];
          g1 = MatrixPlot[Reverse[offprob], FrameLabel -> {"Genotype", "SNP index"}, 
            ColorFunctionScaling -> False, MaxPlotPoints -> Infinity, 
            AspectRatio -> 1/2.5, DataRange -> {{1, Dimensions[offprob][[2]]}, {1, 
               Length[offprob]}}];
          g2 = If[ truefgldiplo =!= None,
                   ListPlot[Transpose[{Range[Length[truegeno[[line]]]], truegeno[[line]]}],
                     PlotRange -> {{1, Dimensions[offprob][[2]]}, {1, Length[offprob]}}, PlotRangePadding -> 0, 
                     PlotMarkers -> {"\[Times]", 10}, PlotStyle -> Red],
                   ListPlot[{1}, PlotMarkers -> {Automatic, 0}]
               ];          
          If[ truefgldiplo =!= None,
          	g3 = ListLinePlot[Thread[{#, {0.3, Length[offprob] + 0.7}}] & /@Accumulate[Tally[truefgldiplo[[3, 2 ;;]]][[All, 2]]], 
            PlotStyle -> Directive[Dashed, Thick, Brown]];
          	gg={g1,g2,g3},
          	gg={g1,g2}
          ];          
          Show[Sequence@@gg, opts,FrameTicks -> {{Transpose[{Range[Length[types[[1, 1]]]], types[[1, 1]]}], None}, {Automatic, None}},
           ImageSize -> 1000, Frame -> True, PlotLabel -> label, 
           LabelStyle -> Directive[FontSize -> 14, Black, FontFamily -> "Helvetica"]], {line, Length[condprob]}];
        gg
    ]
    
plotAncestryProb[summaryfile_String?FileExistsQ, inputtruefgldiplo_?(ListQ[#] ||StringQ[#]||#===None&),
	opts : OptionsPattern[]] :=
    Module[ {summary, nfounder, nfgl,noffspring, truefgldiplo = inputtruefgldiplo, 
      isgeno,chrsubset, gg, pos, condprob, types, line,isfounderinbred,offspringls},
        {isgeno,chrsubset} = OptionValue@{isPlotGenoProb,linkageGroupSet};        
        If[truefgldiplo=!=None,
	        If[ StringQ[truefgldiplo],
	            If[ ! FileExistsQ[truefgldiplo],
	                Print["File ", truefgldiplo, " does not exist!"];
	                Return[$Failed]
	            ];
	            truefgldiplo = Import[truefgldiplo, "CSV", Path -> Directory[]];
	        ];
	        truefgldiplo = getsubMagicSNP[truefgldiplo,chrsubset,All];
        ];        
        (*extract summary file*)
        summary = getSummaryMR[summaryfile];
        (*{"Genetic map and founder hapolotypes", "Conditonal genotype probability", 
        "Conditonal haplotype probability", "Conditonal IBD probability"}*)
        If[isgeno === Automatic, isgeno = ! MatchQ[summary[[5, 1, 1]], "None" | None]];
        If[ MatchQ[summary[[5, 1, 1]], "None" | None] && TrueQ[isgeno],
            Print["Genotype probabilities do not exist!"];
            Return[$Failed]
        ];
        If[isgeno,
    		summary[[{2, 5, 7, 8}]] =getchrsetdata[#, chrsubset] & /@ summary[[{2, 5, 7, 8}]],
    		summary[[{2, 7}]] =getchrsetdata[#, chrsubset] & /@ summary[[{2, 7}]]
    	];
        nfounder = Length[summary[[6, 2 ;;, 2]]];
        (*noffspring = (Length[summary[[7]]] - 3)/nfounder;*)
        offspringls = StringDelete[summary[[7]][[4 ;; ;; nfounder, 1]], "_haplotype1"];
        noffspring = Length[offspringls];
        (*to infer isfounderinbred from truefgl*)    
        isfounderinbred = True;        
        nfgl = nfounder(1 + Boole[! isfounderinbred]); 
        If[ isgeno,
        	(*genoprob*)        	
        	types = origGenotype[nfounder];        	
        	condprob = Table[
        		pos = Span @@ {4 + Length[types[[1, 1]]] (line - 1), 3 + Length[types[[1, 1]]] line};
        		summary[[5, pos, 2 ;;]], {line, noffspring}],        	
        	(*haploprob*)
        	condprob = Table[
        		pos = Span @@ {4 + nfgl (line - 1), 3 + nfgl line};
        		summary[[7, pos, 2 ;;]], {line, noffspring}];        	
        ];        
        If[ isgeno,
        	(*genoprob*)        	
        	gg = plotgenoprob[condprob, truefgldiplo,offspringls,FilterRules[{opts}, Options[ListPlot]]],
        	(*haploprob*)
        	gg = plothaploprob[condprob, truefgldiplo,offspringls,FilterRules[{opts}, Options[ListPlot]]];
        ];         
        gg
    ]

plotAncestryProbGUI[summaryfile_String?FileExistsQ, truefgldiplo_?(ListQ[#] ||StringQ[#]||#===None&),
	opts : OptionsPattern[]] :=
    Module[ {gg},
        gg = plotAncestryProb[summaryfile,truefgldiplo,FilterRules[{opts}, Options[plotAncestryProb]]];      
        ListAnimate[gg, Sequence @@ FilterRules[{opts}, Options[ListAnimate]], 
         AnimationRunning -> True, AnimationDirection -> ForwardBackward, DefaultDuration -> 2 Length[gg]]
    ]
    
Options[plotCondProb] = Join[{linkageGroupSet->All},Options[ListPlot],Options[ListAnimate]];

plotCondProb[probfile_String?FileExistsQ, opts : OptionsPattern[]] :=
    plotCondProb[probfile,None,opts]
    
plotCondProb[probfile_String?FileExistsQ, inputtruefgldiplo_?(ListQ[#] ||StringQ[#]||#===None&),
	opts : OptionsPattern[]] :=
    Module[ {summary, nfounder, nfgl,noffspring, truefgldiplo = inputtruefgldiplo, 
      chrsubset, gg, condprob,nstate,type,prob,isfounderinbred},
        {chrsubset} = OptionValue@{linkageGroupSet};
        (*true*)
        If[ StringQ[truefgldiplo],
            If[ ! FileExistsQ[truefgldiplo],
                Print["File ", truefgldiplo, " does not exist!"];
                Return[$Failed]
            ];
            truefgldiplo = Import[truefgldiplo, "CSV", Path -> Directory[]];
        ];
        truefgldiplo = getsubMagicSNP[truefgldiplo,chrsubset,All];
        nfounder = truefgldiplo[[1, 2]];
        (*to infer isfounderinbred from truefgl*)    
        isfounderinbred = True;        
        nfgl = nfounder(1 + Boole[! isfounderinbred]);
        noffspring = Length[truefgldiplo] - nfounder - 4;
        (*condprob*)
        summary = Import[probfile, "CSV", Path -> Directory[]];
		nstate = (Length[summary] - 3)/noffspring;
		If[Head[nstate] =!= Integer,
		  Print["Inconsistent dimensions. nrow=", Length[summary]];
		  Abort[]
		  ];
		(*Print["{nfgl,noffspring,nstate} = ", {nfgl, noffspring, nstate}];*)
		type = Which[nstate == nfgl*nfgl, "diplotype", 
		   nstate == nfgl*(nfgl + 1)/2, "genotype", nstate == nfgl, 
		   "haplotype",
		   _,
		   Print["inconsistent dimenions. {nfgl,noffspring,nstate} = ", {nfgl,
		      noffspring, nstate}]; Abort[]
		   ];
		prob = SplitBy[Transpose[summary[[All, 2 ;;]]], #[[2]] &];
		prob = Flatten[prob[[chrsubset]], 1][[All, 4 ;;]];
		condprob = Partition[Transpose[prob], nstate];
		Which[ 
			type=="genotype",
        	(*genoprob*)        	
        	gg = plotgenoprob[condprob, truefgldiplo],
        	type=="haplotype",
        	(*haploprob*)
        	gg = plothaploprob[condprob, truefgldiplo],
        	type=="haplotype",
        	Abort[];
        ]; 
		(**)
        ListAnimate[gg, Sequence @@ FilterRules[{opts}, Options[ListAnimate]], 
         AnimationRunning -> True, AnimationDirection -> ForwardBackward, DefaultDuration -> 2 noffspring]
    ]
          
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicReconstruct`*"];

EndPackage[]

