(* Mathematica Package *)

(* Created by the Wolfram Workbench Feb 4, 2016 *)

(* :Author: Chaozhi Zheng <chaozhi@gmail.com>*)
(* :Description: A package for genotype imputation in experimental populations*)
(* :KeyFunctions: magicImpute*)

BeginPackage["MagicImpute`",{"ContinuousTimeHmm`","MagicDefinition`","MagicGeneDropping`",
    "MagicReconstruct`","MagicDataPreprocess`", "MagicReconstruct`MagicModel`","MagicReconstruct`MagicLikelihood`"}]

magicImpute::usage = "magicImpute[magicsnp, model, popdesign] imputes missing genotypes and corrects erroneous genotypes in multi-parental populations. The magicsnp specifies the input genotypic data matrix or filename. The model specifies whether the maternally and patermally derived chromosomes are indepdent (\"indepModel\"), completely dependent (\"depModel\"), or  modeled jointly (\"jointModel\"). The popdesign speficies the breeding design information in several possible ways: a list of mating schemes from founder population to the last generation, a list of values denoting the junction distribution, or a filename for population pedigree information. "

magicImputeFounder::usage = "magicImputeFounder is used only internally."

magicImputeOffspring::usage = "magicImputeFounder is used only internally."

imputingTarget::usage = "imputingTarget is an option to specify the imputing target: \"Founders\", \"Offspring\", or \"All\". "

calImputeAccuracy::usage = "calImputeAccuracy[obsmagicsnp, truemagicsnp, estmagicsnp,excludehalfimpute] returns {{nImputedGeno_Founder, nMissingGeno_Founder}, {nCorrectlyImputedGeno_Founder, nImputedGeno_Founder}, {nImputedGeno_Offspring, nMissingGeno_Offspring}, {nCorrectlyImputedGeno_Offspring, nImputedGeno_Offspring}}, corresponding to the calculations of imputation fraction of founder genoytpes, imputation accuracy of founder genotypes, imputation fraction of offspring genotypes, and imputation accracy of offspring genotypes, respectively. If there are multiple linage groups, each pair such as {nImputedGeno_Founder, nMissingGeno_Founder} is replaced by a list of pairs{{nImputedGeno_Founder_LG1, nMissingGeno_Founder_LG1}, {nImputedGeno_Founder_LG2, nMissingGeno_Founder_LG2}, ...} for all linkage groups. By default, excludehalfimpute = True so that half imputed genotypes such as 1N and N2 are regarded as non-imputed."

calErrorDetection::usage = "calErrorDetection[obsmagicsnp, truemagicsnp, estmagicsnp] or calErrorDetection[obsmagicsnp, truemagicsnp, estmagicprob,detectingthreshold] returns {ntruedetecting, ntruecorrecting, ndetecting, nobserror, nobsgeno}, denoting number of true detected errors, number of true corrected errors, number of observed genotypes that are different from estimated genotypes, number of observed genotypes that are different from true genotypes, and number of observed genotypes. "

calPhasingAccuracy::usage = "calPhasingAccuracy[truemagicsnp, estmagicsnp] returns {{{nNonSwitchError_Founder, nChance_Founder}, {nCorrectingPhasing_Founder, nChance_Founder}}, {{nNonSwitchError_Offspring, nChance_Offspring}, {nCorrectingPhasing_Offspring, nChance_Offspring}}}, corresponding to the calculations of switch phasing accuracy in founders, accuracy of phasing heterozygotes in founders, switch phasing accuracy in offspring, and accuracy of phasing heterozygotes in offspring. For a given individual in a given linkage group, the number of chances is given by the number of heterozygotes minus 1 if there are at least one heterozygoes, and 0 otherwise. "

calHeteroAccuracy::usage = "calHeteroAccuracy[truemagicsnp, estmagicsnp] returns {{nCorrectlyEstimatedHetero, nEstimatedHetero}, {nWronglyEstimatedHetero, nTrueHetero}}, corresponding to the calculations of the accuracy of estimated heterogenotypes (unphased) and the fraction of true heterogenotypes that are undetected. "

genoErrorPattern::usage = "genoErrorPattern[obsmagicsnp, estmagicsnp,truemagicsnp] returns the status of each genotype estimation. The three arguments correspond to observed data, estimated genotypes, and true genotypes, respectively. There are 6 labeled statuses: \"TrueCorrect\" (genotype errors that are changed correctly), \"TrueDetect\"(genotype errors that are changed wrongly), \"FalseNegative\" (genotype errors that are not detected), \"FalsePositive\" (Non genotype errors that are wrongly corrected), and \"FalseImpute\" (wrongly imputed genotypes). \n genoErrorPattern[obsmagicsnp, estmagicsnp] returns results with three possible statuses: \"NonImputation\" (missing genoytypes are not imputed), \"Imputation\" (missing genoytypes are imputed), and \"Correction\" (observed genoytypes are changed). "

plotErrorPatternGUI::usage = "plotErrorPatternGUI[obsmagicsnp, estmagicsnp,truemagicsnp] plots and summarizes the status of each genoetype estimation. The three arguments correspond to observed data, estimated genotypes, and true genotypes, respectively. There are 6 labeled statuses: \"TrueCorrect\" (genotype errors that are changed correctly), \"TrueDetect\"(genotype errors that are changed wrongly), \"FalseNegative\" (genotype errors that are not detected), \"FalsePositive\" (Non genotype errors that are wrongly corrected), and \"FalseImpute\" (wrongly imputed genotypes). \n plotErrorPatternGUI[obsmagicsnp, estmagicsnp] plots and summarizes the results with three possible statuses: \"NonImputation\" (missing genoytypes are not imputed), \"Imputation\" (missing genoytypes are imputed), and \"Correction\" (observed genoytypes are changed). "

parentForwardCalculation::usage = "parentForwardCalculation  "

parentBackwardMaximize::usage = "parentBackwardMaximize  "



Begin["`Private`"]
(* Implementation of the package *)
  
(*Total[Total[origprob priorhaploweight], {2}] == Total[fworigscale priorhaploweight]*)
calPosteriorProb[prob_, priorhaploweight_] :=
    Module[ {pmax,origscale, phaseprob, origprob},
    	pmax = Map[Max, prob, {2}];
    	origscale = Total[prob/pmax, {3}] pmax;
        (*origscale = Total[prob, {3}];*)
        origprob = prob/origscale;
        phaseprob = Total[Log[origscale], {2}] + Log[priorhaploweight];
        phaseprob -= Max[phaseprob];
        (*phaseprob = Replace[phaseprob-Max[phaseprob], x_ /; x < -50 -> -Infinity, {1}];*)
        phaseprob = Normalize[Exp[phaseprob], Total];
        {phaseprob, origprob}
    ]      
        
parentForwardCalculation[model_,samplelabel_, markovprocess_, inputfhaploset_, obsgeno_, epsF_, eps_, minphredscore_,isoffspringdepth_, ismaleX_,isprint_] :=
    Module[ {isdepModel,fhaploset = inputfhaploset,nseq, startprob,tinyprob,fwphaseindex,fwphaseprob, fworigprob, sitedataprob, tranprob, samplecode,t,delay},
        nseq = Length[fhaploset];
        startprob = sampleStartprob[samplelabel, markovprocess];
        fwphaseprob = fworigprob = fwphaseindex = Table[0, {nseq}];
        (*sitedataprob dimensions:{nphase[[1]],noffspring,nstate}*)
        If[ isoffspringdepth,
            sitedataprob = siteMagicLikelihoodGBS[model,fhaploset[[1, All, 1]], obsgeno[[1]], epsF,eps, minphredscore, ismaleX],
            sitedataprob = siteMagicLikelihood[model,fhaploset[[1, All, 1]], obsgeno[[1]], epsF,eps, ismaleX];
        ];
        (*Put[isoffspringdepth,model,fhaploset, obsgeno, epsF,eps, minphredscore, ismaleX,sitedataprob,"tempimpute.txt"];
        Abort[];*)
        isdepModel = ToLowerCase[model]==="depmodel";
        samplecode = sampleCode[samplelabel, isdepModel, ismaleX];              
        fworigprob[[1]] = Map[# startprob &, sitedataprob, {1}];
        {fwphaseprob[[1]], fworigprob[[1]]} = calPosteriorProb[fworigprob[[1]], fhaploset[[1, All, 2]]];
        fwphaseindex[[1]] = Range[Length[fhaploset[[1]]]];        
        (*Put[model,fhaploset, obsgeno, epsF,eps, ismaleX,sitedataprob,startprob,samplelabel, markovprocess,samplecode,fwphaseprob, fworigprob, fwphaseindex,"imputetemp.txt"];
        Abort[];*)     
        delay = If[isprint, 0, 10^6.];     
        Monitor[Do[                   
          tranprob = sampleTransition[samplelabel, markovprocess,samplecode, t - 1];
          If[ isoffspringdepth,
              sitedataprob = siteMagicLikelihoodGBS[model,fhaploset[[t, All, 1]], obsgeno[[t]],epsF, eps, minphredscore, ismaleX],
              sitedataprob = siteMagicLikelihood[model,fhaploset[[t, All, 1]], obsgeno[[t]], epsF,eps, ismaleX];
          ];          
          fworigprob[[t]] = MapThread[#1.#2 &, {fwphaseprob[[t-1]].fworigprob[[t-1]], tranprob}];
          fworigprob[[t]] = Table[fworigprob[[t]],{Length[sitedataprob]}] sitedataprob;  
                  
          {fwphaseprob[[t]], fworigprob[[t]]} = calPosteriorProb[fworigprob[[t]], fhaploset[[t, All, 2]]];          
		  tinyprob = Max[fwphaseprob[[t]]] 10^(-10.);
		  (*tinyprob = Max[fwphaseprob[[t]]] 10^(-20.);*)
          fwphaseindex[[t]] = Flatten[Position[Sign[Round[fwphaseprob[[t]], tinyprob]], 1]];
          {fwphaseprob[[t]], fworigprob[[t]], fhaploset[[t]]} = {fwphaseprob[[t]], fworigprob[[t]], fhaploset[[t]]}[[All, fwphaseindex[[t]]]];
          fwphaseprob[[t]] = Normalize[fwphaseprob[[t]],Total];
          0, {t, 2, nseq}], ProgressIndicator[t, {2, nseq}],delay];
        {fwphaseprob, fworigprob,fwphaseindex}
    ] 
            
(*posMax[list_] :=
    First[FirstPosition[list, Max[list]]]*)
    
(*posMax[list_] :=
    First[Ordering[list, -1]];*)
    
posMax[list_] :=
    RandomChoice[Flatten[Position[list, Max[list]]]]  
        
parentBackwardMaximize[fwphaseprob_, fworigprob_,isdepModel_,ismaleX_,samplelabel_, samplemarkovprocess_,backtmin_,isprint_] :=
    Module[ {nseq, noffspring, fphase, orig, phaseprob,tranprob,ls, weight, ind, t,samplecode,tmin,delay},
        nseq = Length[fwphaseprob];
        noffspring = Length[samplelabel]-1;
        fphase = orig = phaseprob = Table[0,{nseq}];
        phaseprob[[-1]] = Normalize[fwphaseprob[[-1]],Total];
        fphase[[-1]] = posMax[fwphaseprob[[-1]]];
        ls = Normal[fworigprob[[-1]]];
        orig[[-1]] = Table[posMax[ls[[fphase[[-1]], ind]]], {ind, noffspring}];
        samplecode = sampleCode[samplelabel, isdepModel, ismaleX];
        tmin = Min[nseq-1,Max[backtmin,1]];
        delay = If[isprint, 0, 10^6.];
        Monitor[Do[
              tranprob = sampleTransition[samplelabel, samplemarkovprocess,samplecode, t];
              ls = Table[tranprob[[ind, All, orig[[t + 1, ind]]]], {ind, noffspring}];
              ls = Normal[ls # & /@ fworigprob[[t]]];
              weight = Log[Total[ls, {3}]];
              weight = Log[fwphaseprob[[t]]] + Total[weight, {2}];
              weight -= Max[weight];
              (*weight = Replace[weight-Max[weight], x_ /; x < -50 -> -Infinity, {1}];*)
              phaseprob[[t]] =  Normalize[Exp[weight],Total];
              fphase[[t]] = posMax[phaseprob[[t]]];
              orig[[t]] = Table[posMax[ls[[fphase[[t]], ind]]], {ind, noffspring}];
              0, {t, nseq - 1, tmin, -1}], ProgressIndicator[t, {1, nseq - 1}],delay];
        {fphase, phaseprob,orig}
    ]  

(*set epsF=0 if all parental genotypes are missing.*)
calparentepsF[founderHaplo_, isfounderdepth_, epsF_] :=
    Module[ {missingcode, ls},
        missingcode = If[ isfounderdepth,
                          {0, 0},
                          "N" | "NN"
                      ];
        ls = Map[MatchQ[#, missingcode] &, founderHaplo, {3}];
        ls = Boole[And @@ Union[Flatten[#]] & /@ ls];
        epsF (1 - ls)
    ]

magicImputeFounder[magicSNP_, model_, epsF_, eps_, popDesign_,minphredscore_,maxfoundererror_,genothreshold_,
    isfounderinbred_,isfounderdepth_,isoffspringdepth_,starttime_,isprint_,isextrareverse_] :=
    Module[ {isdepModel,deltd, founderHaplo, epsFls,epsF2,obsGeno, foundergender,offspringgender, posA, posX, nFounder,
          founderid, sampleid, snpMap,haploMap,samplelabel, sampleMarkovProcess,markovprocess,chrobsgeno,fhaploset,posfoundermale,
          ismaleX,isfoundermaleX,phase, ch,imputedmagicsnp,revmarkovprocess,fwphaseprob, fworigprob,revphaseprob, revorigprob,
          fwphaseindex,fwphase,revphaseindex,nn,fhaploset2,backtmin},
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,
            foundergender,offspringgender,founderid,sampleid}  = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        (*founderHaplo dimensions change from {nfounder,nchr,nsnp} into {nchr, nsnp, nfounder}*)
        founderHaplo = Transpose[#] & /@ Transpose[founderHaplo];
        epsFls = calparentepsF[founderHaplo, isfounderdepth, epsF];
        (*sampleMarkovProcess is a list of rules for each group (according  memeberID or gender):
          MemberID/gender -> {startProb0,tranProb0}, not accounting for the ordering of funnelcode*)
        (*startProb0=a list of initial distribution for each linkage group;
        tranProb0= a list of transition prob matrices between consecutive markers for each linkage group*)        
        (*Put[nFounder,popDesign, isfounderinbred,model, posA, posX, offspringgender, sampleid,deltd,"impute0.txt"];*)
        {samplelabel, sampleMarkovProcess} = sampleDiscretePriorProcess[nFounder,popDesign, isfounderinbred,model, posA, posX, offspringgender, sampleid,deltd];
        If[Length[offspringgender] =!= Length[samplelabel] - 1,
  			Print["Inconsistent number of offspring between genodata and popdesign!"]; Abort[]
  		];
        isdepModel = ToLowerCase[model]==="depmodel";
        phase = Table[
            epsF2 = epsFls[[ch]];
            markovprocess = sampleMarkovProcess;
            markovprocess[[All, 2]] = sampleMarkovProcess[[All, 2, All, ch]];
            chrobsgeno = Transpose[obsGeno[[All, ch]]];
            If[ ! MemberQ[posX, ch],
                ismaleX = Table[False, {Length[sampleid]}];
                isfoundermaleX = Table[False, {Length[foundergender]}],
                ismaleX = MatchQ[#, "Male" | "male"] & /@ offspringgender;
                isfoundermaleX = MatchQ[#, "Male" | "male"] & /@ foundergender
            ];
            If[ isfounderdepth,
                fhaploset = Map[siteParentHaploPriorGBS[#, epsF2,minphredscore,maxfoundererror,genothreshold,isfoundermaleX,isfounderinbred] &, founderHaplo[[ch]]],
                fhaploset = Map[siteParentHaploPrior[#, epsF2,isfoundermaleX,isfounderinbred,maxfoundererror] &, founderHaplo[[ch]]];
            ];               
            If[Union[Length[#] & /@ fhaploset] === {1},
            	phase = fhaploset[[All, 1, 1]],
	            If[ isprint,
	                Print["Time elapsed = " <>ToString[Round[SessionTime[] - starttime, 0.1]] 
	                <> " Seconds. \t Start imputing founder linkagegroup " <> ToString[ch]<>" out of "<>ToString[Length[founderHaplo]]];	                
	                Print["#markers: ", Length[fhaploset],". #phases per locus: " <> ToString[Round[Mean[Length[#] & /@ fhaploset],0.01]]," out of "<>ToString[2^(nFounder (1+Boole[!isfounderinbred]))]];
	            ];	            
	            nn = Round[Length[fhaploset]/2];	            	            
	            (*Put[model,samplelabel, markovprocess, fhaploset, chrobsgeno, 
	               epsF2, eps, minphredscore, isoffspringdepth, ismaleX,isprint,"imputetemp.txt"];
	            Abort[];*)
	            {fwphaseprob, fworigprob, fwphaseindex} = parentForwardCalculation[model,samplelabel, markovprocess, fhaploset, chrobsgeno, 
	               epsF2, eps, minphredscore, isoffspringdepth, ismaleX,isprint];
	            (*Put[fwphaseprob, fworigprob, fwphaseindex,model,samplelabel, markovprocess, fhaploset, chrobsgeno, epsF2, eps, minphredscore, isoffspringdepth, ismaleX,isprint,"imputetemp2.txt"];
	            Abort[];*)            
	            backtmin = If[ isextrareverse,
	                           nn+1,
	                           1
	                       ];	                     
	            {fwphase,fwphaseprob} = Most[parentBackwardMaximize[fwphaseprob, fworigprob, isdepModel,ismaleX,samplelabel, markovprocess,backtmin,isprint]];	            	            
	            fhaploset = MapThread[#1[[#2]] &, {fhaploset,fwphaseindex}];
	            If[ isprint,
	                PrintTemporary["Mean number of phases per locus after forward calculation: " <> ToString[Round[Mean[Length[#] & /@ fwphaseindex],0.01]]," out of "<>ToString[2^(nFounder (1+Boole[!isfounderinbred]))]];
	            ];
	            If[ !isextrareverse,
	                phase = fwphase,
	                (*fixing phasing of the second half of the chromosome*)
	                fhaploset2 = fhaploset;
	                fhaploset2[[nn + 1 ;;]] = List /@ MapThread[#1[[#2]] &, {fhaploset[[nn + 1 ;;]],fwphase[[nn + 1 ;;]]}];
	                fhaploset2[[nn + 1 ;;, All, 2]] = 1;
	                (*reverse chromosome direction*)
	                revmarkovprocess = markovprocess;	                	
	                revmarkovprocess[[All, 2, 2]] = Reverse[#] & /@ revmarkovprocess[[All, 2, 2]];
	                revmarkovprocess[[All, 2, 2]] = Map[Transpose, revmarkovprocess[[All, 2, 2]], {2}];
	                {revphaseprob, revorigprob, revphaseindex} = parentForwardCalculation[model,samplelabel, revmarkovprocess, 
	                    Reverse[fhaploset2], Reverse[chrobsgeno], epsF2, eps, minphredscore, isoffspringdepth, ismaleX,isprint];
	                (*Put[fhaploset,revphaseindex,revphaseprob, revorigprob, isdepModel,ismaleX, samplelabel, revmarkovprocess,"tempphase0_"<>ToString[$KernelID]<>".txt"];*)	                	
	                {phase,revphaseprob} = Most[parentBackwardMaximize[revphaseprob, revorigprob, isdepModel,ismaleX, samplelabel, revmarkovprocess,Length[fhaploset]+1-nn,isprint]];
	                {revphaseindex,phase,revphaseprob} = Reverse[#]&/@{revphaseindex,phase,revphaseprob};
	                ClearAll[revorigprob];
	                phase = MapThread[#1[[#2]] &, {revphaseindex, phase}];
	                phase[[nn+1;;]] = fwphase[[nn+1;;]];
	            ];
	            phase = MapThread[#1[[#2]] &, {fhaploset[[All, All, 1]], phase}];
	            phase
            ], {ch, Length[founderHaplo]}];
        If[ ! isfounderinbred,        	        	
            posfoundermale = Flatten[Position[foundergender, "Male"]];
            phase[[posX, All, 2 posfoundermale]] = "";
            phase = Map[StringJoin, Map[Partition[#, 2] &, phase, {2}], {3}];
        ];        	
        imputedmagicsnp = magicSNP;
        imputedmagicsnp[[5 ;; nFounder + 4, 2 ;;]] = Transpose[Flatten[phase,1]];
        (*DeleteFile["temporary_phase_chr" <> ToString[#] <> ".txt"] & /@ Range[Length[founderHaplo]];*)         	
        imputedmagicsnp
    ]
  

(*genoprob with dimensions:{noffspring,nlinkagegroup,nsnp,nstate=36},irregular because nsnp depdens on linkagegroup*)
(*return genoprob with dimensions {noffspring,nlinkagegroup,nsnp,nstate=4, because of biallelic snp at diploid loci}*)
calSNPGenoProb[reconstructoutputfile_,model_,founderHaplo_, epsF_, eps_, offspringgender_, posX_] :=
    Module[ {condprob,nFounder, nchr, posA, genotypes, derGeno, malepos, nonmalepos, haplopos, 
        derivedhaplo, i,isdepmodel,offspring},
        {nFounder, nchr} = Take[Dimensions[founderHaplo], 2];
        posA = Complement[Range[nchr], posX];
        genotypes = origGenotype[nFounder][[1, 1]];
        isdepmodel = ToLowerCase[model]==="depmodel";
        derGeno = getDerivedGenotype[founderHaplo, isdepmodel,False];
        If[ isdepmodel,
            (*depModel*)
            condprob = getCondProb[reconstructoutputfile];
            derivedhaplo = derGeno /. {6 -> 2, 9 -> 3};
            condprob = posteriorTrueHaplo[condprob, derivedhaplo,epsF, eps],
            (*indepModel or jointModel*)
            condprob = toGenoProb[getCondProb[reconstructoutputfile]];
            malepos = Flatten[Position[offspringgender, "Male"]];
            nonmalepos = Complement[Range[Length[offspringgender]], malepos];
            haplopos = Flatten[Position[genotypes, #, {1}, Heads -> False] & /@Table[{i, i}, {i, nFounder}]];
            If[ posA=!={},
                offspring = Partition[Range[Length[condprob]],UpTo[50]];
                Do[condprob[[i, posA]] = posteriorTrueDiplo[condprob[[i, posA]], derGeno[[posA]], genotypes, epsF, eps],{i,offspring}]
            ];
            If[ (nonmalepos =!= {}) && (posX =!= {}),
                condprob[[nonmalepos, posX]] = posteriorTrueDiplo[condprob[[nonmalepos, posX]],derGeno[[posX]], genotypes, epsF,eps];
            ];
            If[ (malepos =!= {}) && (posX =!= {}),
                derivedhaplo = derGeno[[posX, All, haplopos]] /. {6 -> 2, 9 -> 3};
                condprob[[malepos, posX]] = posteriorTrueHaplo[condprob[[malepos, posX, All, haplopos]],derivedhaplo, epsF, eps];
            ];
        ];
        condprob
    ]
    
imputeBestSNPGeno[isdepmodel_,obsGeno_, snpgenoprob_, offspringgender_, posX_, errorgenobound_, imputingbound_,isprint_] :=
    Module[ {posNN, posN, posnonmiss, best, malepos, bestprob, prob,rule,geno,detectedpos,
      nonimputeposNN, nonimputeposN, errorpos, bestgeno,errorcount,nonimputecount,avgimputeprob,temp,pos,missingcode},
        best = maxIndexPair[snpgenoprob];
        (*Print[Dimensions[best],"; ", Dimensions[best[[All,1]]]];*)
        bestprob = best[[All, All, All, 1]];
        bestgeno = best[[All, All, All, 2,1]];
        ClearAll[best];
        (*set imputation*)
        posNN = Position[obsGeno, "NN"|"1N"|"2N"|"N1"|"N2"];
        posN = Position[obsGeno, "N"];
        posnonmiss = Complement[Position[obsGeno, _, {3}, Heads -> False], posNN, posN];
        If[ isdepmodel,
            bestgeno = bestgeno /. Thread[Range[2] -> {"11", "22"}],
            bestgeno = bestgeno /. Thread[Range[4] -> {"11", "12", "21", "22"}];
        ];
        If[ posX =!= {},
            malepos = Flatten[Position[offspringgender, "Male"]];
            bestgeno[[malepos, posX]] = bestgeno[[malepos, posX]] /. {"11" -> "1", "12"|"22" -> "2"};
        ];
        nonimputeposNN = Pick[posNN, Thread[Extract[bestprob, posNN] <= imputingbound]];
        (*bestgeno = ReplacePart[bestgeno, Thread[nonimputeposNN -> "NN"]];*)
        If[nonimputeposNN =!= {},
	        If[ isdepmodel,
	            bestgeno = ReplacePart[bestgeno, nonimputeposNN -> "NN"];
	            bestprob = ReplacePart[bestprob, nonimputeposNN -> "NA"],            
	            temp = Transpose[Extract[snpgenoprob, nonimputeposNN]];
	            (*{11,12,21,22}  <->  {{1,2},{3,4},{1,3},{2, 4}} <-> {"1N","2N","N1","N2"}*)
	            missingcode = {"1N", "2N", "N1", "N2"};
	            temp = Transpose[Total[temp[[#]]] & /@ {{1, 2}, {3, 4}, {1, 3}, {2, 4}}];
	            temp = maxIndexPair[temp];
	            temp[[All, 2]] = missingcode[[temp[[All, 2, 1]]]];
	            pos = Flatten[Position[Boole[Thread[temp[[All, 1]] > imputingbound]] Boole[Thread[(Length[#] & /@ temp[[All, 2]]) <= 1]], 0]];
	            temp[[pos, 1]] = "NA";
	            temp[[pos, 2]] = "NN";
	            bestgeno = ReplacePart[bestgeno, Thread[nonimputeposNN -> temp[[All, 2]]]];
	            bestprob = ReplacePart[bestprob, Thread[nonimputeposNN -> temp[[All, 1]]]];
	        ];
        ];
        If[ posN =!= {},
            nonimputeposN = Pick[posN, Thread[Extract[bestprob, posN] <= imputingbound]];
            bestgeno = ReplacePart[bestgeno, Thread[nonimputeposN -> "N"]];
        ];
        (*error detection*)
        errorpos = Pick[posnonmiss,MapThread[!MemberQ[#1, #2] &, ({Replace[Extract[obsGeno, posnonmiss], {"12" | "21" -> {"12", "21"},x_ :> {x}}, {1}], 
            Extract[bestgeno, posnonmiss]})]];
        detectedpos = Flatten[Position[Thread[Extract[bestprob, errorpos] > errorgenobound], True]];
        If[ isdepmodel,
            rule = Join[Thread[{"11", "22"} -> Range[2]],Thread[{"1", "2"} -> Range[2]], {"12" | "21" -> 0}],
            rule = Join[Thread[{"11", "12", "21", "22"} -> Range[4]],Thread[{"1", "2"} -> Range[2]]];
        ];
        (*set detected errorenous genotype as the best genotype, and the rest non-missing genotype as the observed genotypes*)
        pos = Complement[errorpos, errorpos[[detectedpos]]];
        If[ pos=!={},
            geno = Replace[Extract[obsGeno, pos], {"12" | "21" -> {"12", "21"},x_ :> {x}}, {1}];
            prob = Replace[MapThread[#1[[#2]] &, {Extract[snpgenoprob, pos], geno /. rule}],List -> 0, {2}];
            temp = Flatten[Ordering[#, -1] & /@ prob];
            geno = MapThread[#1[[#2]] &, {geno, temp}];
            prob = MapThread[#1[[#2]] &, {prob, temp}];
            bestprob = ReplacePart[bestprob, Thread[pos -> prob]];
            bestgeno = ReplacePart[bestgeno, Thread[pos -> geno]];
        ];
        (*If[ posX =!= {},
            bestgeno[[malepos, posX]] = bestgeno[[malepos, posX]] /. {"NN" -> "N"};
        ];*)
        errorpos = errorpos[[detectedpos]];
        nonimputecount = {N[Count[Flatten[Characters[Extract[bestgeno, nonimputeposNN]]], "N"]/2] + Length[nonimputeposN],
            N[Count[Flatten[Characters[Extract[obsGeno, posNN]]], "N"]/2] + Length[posN]};
        errorcount = {Length[errorpos],Length[posnonmiss]};
        avgimputeprob = Mean[Extract[bestprob, Complement[posNN,nonimputeposNN]]];
        {bestgeno, bestprob,errorpos,errorcount,nonimputecount,avgimputeprob}
    ]

imputeBestSNPGenoGBS[isdepmodel_,obsGeno_,minphredscore_,snpgenoprob_, offspringgender_, posX_, errorgenobound_,imputingbound_,isprint_] :=
    Module[ {rule,best, bestprob, bestgeno, malepos, posnonimpute, posimpute,errorpos,calledgeno,errorcount,nonimputecount,avgimputeprob,
        missingcode,temp,pos,rawwcallthreshold = 0.95},
        If[ isdepmodel,
            rule = Thread[Range[2] -> {"11", "22"}],
            rule = Thread[Range[4] -> {"11", "12", "21", "22"}];
        ];
        best = maxIndexPair[snpgenoprob];
        bestprob = best[[All, All, All, 1]];
        bestgeno = best[[All, All, All, 2, 1]] /. rule;
        ClearAll[best];
        If[ posX =!= {},
            malepos = Flatten[Position[offspringgender, "Male"]];
            bestgeno[[malepos, posX]] = bestgeno[[malepos, posX]] /. {"11" -> "1", "12"|"22" -> "2"};
        ];
        posnonimpute = Position[bestprob - imputingbound, _?NonPositive, {3}, Heads -> False];
        posimpute = Position[bestprob - imputingbound, _?Positive, {3}, Heads -> False];
        If[ isdepmodel,
            bestgeno = ReplacePart[bestgeno, posnonimpute ->"NN"];
            bestprob = ReplacePart[bestprob, posnonimpute -> "NA"],
            temp = Transpose[Extract[snpgenoprob, posnonimpute]];
            (*{11,12,21,22}<->{{1,2},{3,4},{1,3},{2,4}}<->{"1N","2N","N1","N2"}*)
            missingcode = {"1N", "2N", "N1", "N2"};
            temp = Transpose[Total[temp[[#]]] & /@ {{1, 2}, {3, 4}, {1, 3}, {2, 4}}];
            temp = maxIndexPair[temp];
            temp[[All, 2]] = missingcode[[temp[[All, 2, 1]]]];
            pos = Flatten[Position[Boole[Thread[temp[[All, 1]] > imputingbound]] Boole[Thread[(Length[#] & /@ temp[[All, 2]]) <= 1]], 0]];
            temp[[pos, 1]] = "NA";
            temp[[pos, 2]] = "NN";
            bestgeno = ReplacePart[bestgeno, Thread[posnonimpute ->temp[[All, 2]]]];
            bestprob = ReplacePart[bestprob, Thread[posnonimpute ->temp[[All, 1]]]];
        ];
        If[ posX =!= {},
            bestgeno[[malepos, posX]] = bestgeno[[malepos, posX]] /. {"NN" -> "N", "1N" -> "1", "2N" -> "2", "N1" -> "N", "N2" -> "N"};
        ];
        (*error detection*)
        {errorpos,calledgeno} = geterrorposGBS[bestgeno,obsGeno,minphredscore,malepos, posX,rawwcallthreshold];
        errorpos = errorpos[[Flatten[Position[Thread[Extract[bestprob, errorpos] > errorgenobound], True]]]];
        (*The detected wrong genotypes have been corrected in bestgeno*)
        (*bestprob = ReplacePart[bestprob, Thread[errorpos -> "NA"]];
        bestgeno = ReplacePart[bestgeno, Thread[errorpos -> "NN"]];
        If[ posX =!= {},
            bestgeno[[malepos, posX]] = bestgeno[[malepos, posX]] /. {"NN" -> "N"};
        ];*)
        nonimputecount = {N[Count[Flatten[Characters[Extract[bestgeno, posnonimpute]]], "N"]/2],Length[posimpute] + Length[posnonimpute]};
        errorcount = {Length[errorpos],Length[posimpute] + Length[posnonimpute]};
        avgimputeprob = Mean[Extract[bestprob, posimpute]];
        {bestgeno, bestprob,errorpos,calledgeno,errorcount,nonimputecount,avgimputeprob}
    ]   
            
geterrorposGBS[bestgeno_,observedAD_,minphredscore_,haploidIndpos_,haploidchrpos_,rawcallthreshold_] :=
    Module[ {calledgeno,posnonmiss,estnonmissgeno,callednonmissgeno,errorpos,rule},
        calledgeno = rawGenotypeCall[observedAD, minphredscore,rawcallthreshold];
        If[ (haploidIndpos=!={})&&(haploidchrpos=!={}),
            calledgeno[[haploidIndpos, haploidchrpos]] = calledgeno[[haploidIndpos, haploidchrpos]] /. {"11"|"1N" | "N1" -> "1","22"|"2N" | "N2" -> "2", "NN" -> "N"};
        ];
        posnonmiss = Position[observedAD, {0, 0}, {3}, Heads -> False];
        posnonmiss = Complement[Position[observedAD, _, {3}, Heads -> False], posnonmiss];
        rule = {"1N" -> {"11", "12"}, "2N" -> {"21", "22"}, 
        "N1" -> {"11", "21"}, "N2" -> {"12", "22"}, 
        "NN" -> {"11", "21", "12", "22"}, "N" -> {"1", "2"}, x_ :> {x}};
        estnonmissgeno = Replace[Extract[bestgeno, posnonmiss], rule, {1}] /. {"21" -> "12"};
        callednonmissgeno = Replace[Extract[calledgeno, posnonmiss], rule, {1}] /. {"21" -> "12"};
        errorpos = Pick[posnonmiss, MapThread[Intersection[#1, #2] === {} &, {callednonmissgeno,estnonmissgeno}]];
        {errorpos,calledgeno}
    ]    
               
magicImputeOffspringSplit[magicSNP_, model_, epsF_, eps_, popDesign_, minphredscore_, errorgenobound_,imputingbound_,
  isfounderinbred_, isfounderdepth_, isoffspringdepth_, starttime_, isprint_] :=
    Module[ {subset, subsnp, imputedmagicsnp, imputedsnpgenoprob, 
      truegenoposterior, errorgeno, errorcount, nonimputecount,avgimputeprob, pos, ch},
     (*subset[[ch]]: {chrname, column indices including column 1*)
        subset = magicSNP[[3, 2 ;;]];
        subset = {#[[1, 2]], Join[{1}, #[[All, 1]] + 1]} & /@ SplitBy[Transpose[{Range[Length[subset]], subset}], Last];
        {imputedmagicsnp, imputedsnpgenoprob, truegenoposterior, errorgeno,errorcount, nonimputecount, avgimputeprob} = Transpose[Table[
           If[ isprint,
               PrintTemporary["Time elapsed = " <> ToString[Round[SessionTime[] - starttime, 0.1]] <> 
                  " Seconds. \t Start imputing offspring linkagegroup " <> ToString[ch]<>" out of "<>ToString[Length[subset]]];
           ];
           subsnp = Join[magicSNP[[{1}]], magicSNP[[2 ;;, subset[[ch, 2]]]]];
           magicImputeOffspring[subsnp, model, epsF, eps, popDesign, minphredscore, errorgenobound, imputingbound, 
            isfounderinbred, isfounderdepth, isoffspringdepth, starttime, False], {ch, Length[subset]}]];
        {imputedmagicsnp, imputedsnpgenoprob, truegenoposterior} = Join[First[#], Sequence @@ #[[2 ;;, All, 2 ;;]], 2] & /@ {imputedmagicsnp, imputedsnpgenoprob, truegenoposterior};
        imputedmagicsnp[[1]] = imputedmagicsnp[[1, ;; 2]];
        imputedsnpgenoprob[[1]] = imputedsnpgenoprob[[1, ;; 2]];
        truegenoposterior[[1]] = truegenoposterior[[1, ;; 2]];
        pos = ToExpression[Map[StringSplit[#, "|"] &, errorgeno[[All, 2 ;;, 1]]]];
        pos[[All, All, 1]] *= Range[Length[pos]];
        errorgeno[[All, 2 ;;, 1]] = toDelimitedString[pos, "|"];
        errorgeno = Join[First[errorgeno],Flatten[errorgeno[[2 ;;, 2 ;;]], 1]];
        {errorcount, nonimputecount} = Total[#] & /@ {errorcount, nonimputecount};
        avgimputeprob = Mean[avgimputeprob];
        {imputedmagicsnp, imputedsnpgenoprob, truegenoposterior, errorgeno,errorcount, nonimputecount, avgimputeprob}
    ]
  
               
magicImputeOffspring[magicSNP_, model_, epsF_, eps_, popDesign_,
    minphredscore_,errorgenobound_,imputingbound_, isfounderinbred_,isfounderdepth_,isoffspringdepth_,starttime_,isprint_] :=
    Module[ {outputid,temporary,reconstructoutputfile,deltd, founderHaplo, obsGeno, foundergender,offspringgender, posA, posX, nFounder, 
          founderid, sampleid, snpMap,haploMap,snpgenoprob,bestgeno, bestprob,calledgeno, errorpos,imputedmagicsnp, imputedsnpgenoprob,
          errorgeno,truegenoposterior,i,errorcount,nonimputecount,avgimputeprob,isdepmodel},
        outputid = ToString[Unique[temporary]];
        reconstructoutputfile = magicReconstruct[magicSNP, model,popDesign,
            founderAllelicError -> epsF,offspringAllelicError -> eps,   
            reconstructAlgorithm -> "origPosteriorDecoding",isFounderInbred -> isfounderinbred,
            sequenceDataOption ->{isOffspringAllelicDepth -> isoffspringdepth, minPhredQualScore -> minphredscore},
            outputFileID->outputid,isPrintTimeElapsed -> False];
        If[ isprint,
            PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \t Start calculating posterior probability of true genotypes! "];
        ];
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,
            founderid,sampleid}  = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];        
        (*Put[reconstructoutputfile,founderHaplo, obsGeno, epsF, eps,minphredscore, offspringgender, posX,isoffspringdepth,errorgenobound,imputingbound,isprint,"temp.txt"];
        Print["calSNPGenoProb"];*)
        (*Abort[];*)        
        snpgenoprob = calSNPGenoProb[reconstructoutputfile,model,founderHaplo, epsF, eps,offspringgender, posX];        
        DeleteFile[reconstructoutputfile];
        snpgenoprob = Map[Normalize[#, Total] &, snpgenoprob, {3}];
        snpgenoprob = Round[snpgenoprob, 10^(-5.)];
        If[ isprint,
            PrintTemporary["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \t Start genotype error detecting and imputing! "];
        ];
        isdepmodel = ToLowerCase[model] == "depmodel";
        (*Put[isdepmodel,obsGeno,minphredscore,snpgenoprob, offspringgender, posX, errorgenobound,imputingbound,isprint,"temp.txt"];*)
        (*Put[isdepmodel,obsGeno, snpgenoprob, offspringgender, posX, errorgenobound,imputingbound,isprint,"temp.txt"];*)
        If[ isoffspringdepth,
            {bestgeno, bestprob,errorpos,calledgeno,errorcount,nonimputecount,avgimputeprob} = imputeBestSNPGenoGBS[isdepmodel,obsGeno,minphredscore,snpgenoprob, 
                offspringgender, posX, errorgenobound,imputingbound,isprint],
            {bestgeno, bestprob,errorpos,errorcount,nonimputecount,avgimputeprob} = imputeBestSNPGeno[isdepmodel,obsGeno, snpgenoprob, 
                offspringgender, posX, errorgenobound,imputingbound,isprint];
            calledgeno = obsGeno;
        ];
        (*Put[bestgeno, bestprob,errorpos,errorcount,nonimputecount,avgimputeprob,isdepmodel,obsGeno, snpgenoprob, 
                offspringgender, posX, errorgenobound,imputingbound,isprint,"imputetemp.txt"];
        Abort[];*)
        imputedmagicsnp = imputedsnpgenoprob = truegenoposterior = magicSNP;
        Do[
            imputedmagicsnp[[nFounder + 4+i, 2 ;;]] = Flatten[bestgeno[[i]]];
            imputedsnpgenoprob[[nFounder + 4+i, 2 ;;]] = Flatten[bestprob[[i]]];
            truegenoposterior[[nFounder + 4+i, 2 ;;]] = toDelimitedString[Flatten[snpgenoprob[[i]],1],"|"];
            0,{i,Length[bestgeno]}];
        (*set errorgeno*)
        errorgeno = geterrorpos[errorpos,obsGeno,calledgeno, bestgeno,snpMap,sampleid,offspringgender,isoffspringdepth,snpgenoprob];
        {imputedmagicsnp, imputedsnpgenoprob,truegenoposterior,errorgeno,errorcount,nonimputecount,avgimputeprob}
    ]

geterrorpos[errorpos_,obsGeno_,calledgeno_, bestgeno_,snpMap_,individualid_,individualgender_,isoffspringdepth_,snpgenoprob_:None] :=
    Module[ {errorgeno,map,chrs,bool,temp,rule},
        map = SplitBy[Transpose[snpMap[[All, 2 ;;]]], #[[2]] &];
        chrs = map[[All, 1, 2]];
        bool = snpgenoprob===None;
        errorgeno = ConstantArray[0, {Length[errorpos] + 1, 8+Boole[!bool]}];
        errorgeno[[1]] = {"Index(chr|marker|ind)","Chromosome", "SNP", "Individual", "Gender",
                                "ObservedGenotype","CalledGenotype",  "ImputedGenotype","DiplotypeProb(11|22 or 11|12|21|22)"};
        If[ bool,
            errorgeno[[1]] = Drop[errorgeno[[1]],-1]
        ];
        If[ errorpos=!={},
            (*{chr, snp,ind} instead of  {ind, chr, snp}*)
            errorgeno[[2 ;;, 1]] = errorpos[[All,{2,3,1}]];
            errorgeno[[2 ;;, 2]] = chrs[[errorpos[[All, 2]]]];
            errorgeno[[2 ;;, 3]] = Extract[map, errorpos[[All, {2,3}]]][[All, 1]];
            errorgeno[[2 ;;, 4]] = individualid[[errorpos[[All, 1]]]];
            errorgeno[[2 ;;, 5]] = individualgender[[errorpos[[All, 1]]]];
            errorgeno[[2 ;;, 6 ;; 7]] = Transpose[Extract[#, errorpos] & /@ {obsGeno, calledgeno}];
            If[ isoffspringdepth,
                errorgeno[[2 ;;, 6]] = toDelimitedString[errorgeno[[2 ;;, 6]], "|"];
            ];
            If[ bool,
                errorgeno[[2 ;;, 8]] = Extract[bestgeno, errorpos],
                temp = Extract[snpgenoprob, errorpos];
                If[ Union[Length[#]&/@temp] ==={2},
                    (*is depModel*)
                    rule = Thread[Range[2] -> {"11", "22"}],
                    (*either jointModel or indepModel*)
                    rule = Thread[Range[4] -> {"11", "12", "21", "22"}]
                ];
                errorgeno[[2 ;;, 9]] = toDelimitedString[temp, "|"];
                errorgeno[[2 ;;, 8]] = maxIndexPair[temp][[All, 2, 1]] /.rule;
                temp = 1+Flatten[Position[errorgeno[[2;;]], _?(MatchQ[#[[{2, 5}]], {"X" | "x","Male" | "male"}] &), {1}, Heads -> False]];
                errorgeno[[temp, 8]] = errorgeno[[temp, 8]] /. {"11"->"1","12"|"22" -> "2"};
            ];
            errorgeno[[2 ;;]] = SortBy[errorgeno[[2 ;;]],First];
            errorgeno[[2 ;;, 1]] = toDelimitedString[errorgeno[[2 ;;, 1]],"|"];
        ];
        errorgeno
    ]
    
getfoundererrorpos[imputedmagicSNP_, magicSNP_, minphredscore_,genothreshold_,isfounderinbred_, isfounderdepth_, isoffspringdepth_] :=
    Module[ {deltd, founderHaplo, obsGeno, snpMap, haploMap, nFounder, rawcallthreshold = 0.95,
      posA, posX, foundergender, offspringgender, founderid, sampleid, rule,ls,
      pos, bestgeno, missingpatt,indicator,obsfoundergeno, haploidIndpos, haploidchrpos, 
      errorpos, calledgeno, posnonmiss,errorgeno,errorcount,nonimputecount},
        {deltd, founderHaplo, obsGeno,snpMap,haploMap,nFounder,posA,posX,foundergender,offspringgender,
            founderid,sampleid}  = transformMagicSNP[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        pos = SplitBy[Transpose[{Range[Length[magicSNP[[3, 2 ;;]]]],magicSNP[[3, 2 ;;]]}], Last][[All, All, 1]] + 1;
        bestgeno = Transpose[imputedmagicSNP[[5 ;; 4 + nFounder, #]] & /@ pos];
        missingpatt = "NN" | "N" | "1N" | "2N" | "N1" | "N2";
        haploidIndpos = 
         If[ isfounderinbred,
             Range[nFounder],
             Flatten[Position[foundergender, "Male"]]
         ];
        haploidchrpos = If[ isfounderinbred,
                            Join[posA, posX],
                            posX
                        ];
        If[ isfounderdepth,
            obsfoundergeno = Transpose[magicSNP[[5 ;; 4 + nFounder, #]] & /@ pos];
            obsfoundergeno = ToExpression[Map[StringSplit[#, "|"] &, obsfoundergeno, {2}]];
            {errorpos, calledgeno} = geterrorposGBS[bestgeno, obsfoundergeno, minphredscore,haploidIndpos,haploidchrpos,rawcallthreshold],
            If[ isfounderinbred,
                obsfoundergeno = founderHaplo,
                obsfoundergeno = Transpose[#] & /@ Partition[founderHaplo, 2];
                obsfoundergeno = Map[StringJoin @@ # &, Map[Transpose, obsfoundergeno, {2}], {3}];
                obsfoundergeno[[haploidIndpos, haploidchrpos]] = Map[StringTake[#, 1] &,obsfoundergeno[[haploidIndpos, haploidchrpos]], {3}];
            ];
            posnonmiss = Position[obsfoundergeno, _?(! StringMatchQ[#, "N" | "NN"] &), {3},Heads -> False];
            rule = {"12"|"21"->{"12", "21"},"1N" -> {"11", "12", "21"}, "2N" -> {"12", "21", "22"}, 
               "N1" -> {"11", "12", "21"}, "N2" -> {"12", "21", "22"}, 
               "NN" -> {"11", "21", "12", "22"}, "N" -> {"1", "2"}, x_ :> {x}};
            ls = Replace[Extract[obsfoundergeno, posnonmiss], rule, {1}];
            errorpos = Pick[posnonmiss,MapThread[MemberQ[#2, #1] &, {Extract[bestgeno, posnonmiss], ls}],False];
            calledgeno = obsfoundergeno;
        ];
        errorgeno = geterrorpos[errorpos, obsfoundergeno, calledgeno, bestgeno, snpMap, founderid,foundergender,isoffspringdepth];
        indicator = Replace[calledgeno, {missingpatt -> 1, _ -> 0}, {3}];
        If[ isfounderdepth,
            errorcount = {Length[errorgeno] - 1, Length[Flatten[indicator]]};
            nonimputecount = {Count[Flatten[bestgeno],missingpatt], Length[Flatten[indicator]]},
            errorcount = {Length[errorgeno] - 1, Total[Flatten[1-indicator]]};
            nonimputecount = {Count[DeleteCases[Flatten[bestgeno indicator], 0],missingpatt], Total[Flatten[indicator]]};
        ];
        {errorgeno,errorcount,nonimputecount}
    ]
      
correcterrorgeno[errorgeno_, imputedmagicsnp_, isfounder_, isfounderinbred_] :=
    Module[ {nfounder, imputedgeno, indices, cond, missingcode,outerrorgeno, pos,outimputedmagicsnp},
        nfounder = imputedmagicsnp[[1, 2]];
        imputedgeno = SplitBy[Transpose[imputedmagicsnp[[2 ;;, 2 ;;]]], #[[2]] &];
        indices = ToExpression[StringSplit[errorgeno[[2 ;;, 1]], "|"]];
        If[ indices=!={},
            If[ isfounder,
                indices[[All, 3]] += 3,
                indices[[All, 3]] += 3 + nfounder
            ];
            (*cond = Extract[imputedgeno, indices] == errorgeno[[2 ;;, 8]];*)
            cond = Extract[imputedgeno, Join[indices[[All, ;; 2]],ConstantArray[1, {Length[indices], 1}], 2]] === errorgeno[[2 ;;, 3]];
            If[ ! cond,
                Print["correcterrorgeno: detected erroneous genotypes does not match imputedmagicsnp!"];
                Abort[]
            ];
        ];
        If[ isfounder && isfounderinbred,
            missingcode = Table["N",{Length[indices]}],
            missingcode = Table["NN",{Length[indices]}];
            pos = Flatten[Position[errorgeno[[2;;]], _?(MatchQ[#[[{2, 5}]], {"X" | "x","Male" | "male"}] &), {1}, Heads -> False]];
            missingcode[[pos]] = "N"
        ];
        imputedgeno = ReplacePart[imputedgeno, Thread[indices -> missingcode]];
        {outerrorgeno, outimputedmagicsnp} = {errorgeno, imputedmagicsnp};
        outimputedmagicsnp[[2 ;;, 2 ;;]] = Transpose[Flatten[imputedgeno, 1]];
        outerrorgeno = Transpose[Insert[Transpose[outerrorgeno], Prepend[missingcode,"ImputedGenotype"], 8]];
        {outerrorgeno, outimputedmagicsnp}
    ]
        
Options[magicImpute] = {
    founderAllelicError -> 0.005,
    offspringAllelicError -> 0.005,  
    isFounderInbred -> True,
    imputingTarget -> "All",
    imputingThreshold -> 0.9,
    detectingThreshold -> 0.9,
    sequenceDataOption ->{
        isFounderAllelicDepth -> Automatic,
        isOffspringAllelicDepth -> Automatic,
        minPhredQualScore -> 30,
        priorFounderCallThreshold -> 0.99
        },
    outputFileID -> "", 
    isPrintTimeElapsed ->True
}

magicImpute[inputmagicSNP_?(ListQ[#] || StringQ[#] &), model_?(VectorQ[#, StringQ]|| StringQ[#]&), inputpopDesign_?(ListQ[#] || StringQ[#] &),opts : OptionsPattern[]] :=
    Module[ {starttime,modells,magicSNP = inputmagicSNP, popDesign = inputpopDesign,epsF,eps,founderimputedmagicsnp,foundererrorgeno,
        foundererrorcount = {-1,-1},errorcount = {-1,-1},foundernonimputecount = {-1,-1},nonimputecount = {-1,-1},avgimputeprob = -1,outputid,isprint,imputetarget,
        foundergenocallbound,errorgenobound,imputingbound,isfounderinbred,isfounderdepth,isoffspringdepth,minphredscore,maxfoundererror = 0,
        outputfiles,imputedmagicsnp,imputedsnpgenoprob,errorgeno,truegenoposterior,labels,res,key = "magicImpute-Summary",summary},
        modells = Flatten[{model}];
        {epsF,eps} = OptionValue@{founderAllelicError,offspringAllelicError};
        {isfounderdepth,isoffspringdepth,minphredscore,foundergenocallbound} =
            OptionValue[Thread[sequenceDataOption -> {isFounderAllelicDepth,isOffspringAllelicDepth,minPhredQualScore,priorFounderCallThreshold}]];
        {isfounderinbred,imputetarget,errorgenobound,imputingbound, outputid,isprint} = 
            OptionValue@{isFounderInbred,imputingTarget,detectingThreshold,imputingThreshold,outputFileID,isPrintTimeElapsed};
        If[ outputid=!="",
            outputid = outputid<>"_"
        ];
        If[ isprint,
            starttime = SessionTime[];
            Print["magicImpute. Start date = ", DateString[]];
        ];
        imputetarget = ToLowerCase[imputetarget];
        If[ !MemberQ[{"founders","founder","offspring","all"},imputetarget],
            Print["magicImpute: imputingTarget option must take \"All\", \"Offspring\", or \"Founders\"."];
            Return[$Failed]
        ];
        If[ !(And @@ (MemberQ[{"jointModel", "indepModel", "depModel"}, #] & /@ modells)),
            Print["magicImpute: model has to take \"jointModel\", \"indepModel\", or \"depModel\"."];
            Return[$Failed]
        ];
        If[ StringQ[magicSNP],
            If[ !FileExistsQ[magicSNP],
                Print["File ", magicSNP," does not exist!"];
                Return[$Failed]
            ];
            magicSNP = Import[magicSNP,"CSV"];
        ];
        If[ StringQ[popDesign],
            If[ !FileExistsQ[popDesign],
                Print["File ", popDesign," does not exist!"];
                Return[$Failed]
            ];
            popDesign = Import[inputpopDesign,"CSV"];
        ];
        {isfounderdepth,isoffspringdepth} = checkOptionValue[magicSNP, isfounderdepth, isoffspringdepth,isfounderinbred,isprint];
        SNPValidation[magicSNP,isfounderinbred,isfounderdepth,isoffspringdepth];
        If[ isprint,
            Print["{#founders, #offspring, #markers} = ", Prepend[Dimensions[Rest[magicSNP]] - {magicSNP[[1, 2]] + 3, 1},magicSNP[[1, 2]]]]
        ];
        If[ MatchQ[imputetarget,"all"|"founders"|"founder"],
            If[ isprint,
                Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \t Start imputing founders! "];
            ];
            founderimputedmagicsnp = magicImputeFounder[magicSNP, First[modells], epsF, eps, popDesign,
                minphredscore,maxfoundererror,foundergenocallbound,isfounderinbred,isfounderdepth,isoffspringdepth,starttime,isprint,True];
            {foundererrorgeno,foundererrorcount,foundernonimputecount} = getfoundererrorpos[founderimputedmagicsnp, magicSNP, minphredscore,errorgenobound,isfounderinbred,isfounderdepth, isoffspringdepth];
            magicSNP = founderimputedmagicsnp;
            isfounderdepth = False;
        ];
        If[ !MatchQ[imputetarget,"founders"|"founder"],
            If[ isfounderdepth,
                Print["magicImpute: required to perform founder genotype imputation to obtained phased founder haplotypes."];
                Abort[]
            ];
            If[ isprint,
                Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \t Start imputing offspring! "];
            ];
            {imputedmagicsnp, imputedsnpgenoprob,truegenoposterior,errorgeno,errorcount,nonimputecount,avgimputeprob} = 
                magicImputeOffspringSplit[magicSNP, Last[modells], epsF, eps, popDesign,minphredscore,errorgenobound,imputingbound,
                isfounderinbred,isfounderdepth,isoffspringdepth,starttime,isprint];
        ];
        If[ isprint,
            Print["Time elapsed = "<>ToString[Round[SessionTime[] - starttime,0.1]]<>" Seconds. \t Start exporting! "];
        ];
        summary = Switch[imputetarget,
            "founders" | "founder", 
            Join[{{key, "Potential erroneous founder genotypes (inconsistent between estimates and called genotypes)"}}, foundererrorgeno],
            "offspring", 
            Join[{{key, "Potential erroneous offspring genotypes (inconsistent between estimates and called genotypes)"}}, errorgeno],
            "all", Join[{{key, "Erroneous founder genotypes"}}, foundererrorgeno,
                         {{key, "Erroneous offspring genotypes"}}, errorgeno],
            _, Abort[];
           ];
        outputfiles = outputid <> # & /@ {"ErroneousGenotype.csv","ImputedGenotype.csv", "PosteriorProbability.csv"};
        If[ MatchQ[imputetarget, "founders" | "founder"],
            res = MapThread[csvExport[#1,#2]&,{Most[outputfiles],{summary,founderimputedmagicsnp}}],
            res = MapThread[csvExport[#1,#2]&,{outputfiles,{summary,imputedmagicsnp,truegenoposterior}}];
        ];        
        (*Print a short summary*)
        If[ isprint,
            labels = {"Fraction of detected erroneous genotypes: ", "Fraction of imputed missing genotyes: "};
            labels[[1]] = labels[[1]] <> ToString[AccountingForm[foundererrorcount[[1]]]] <> "/" <> ToString[AccountingForm[foundererrorcount[[2]]]] <> " in founders";
            labels[[2]] = labels[[2]] <> ToString[AccountingForm[foundernonimputecount[[2]]-foundernonimputecount[[1]]]] <> "/" <> ToString[AccountingForm[foundernonimputecount[[2]]]] <> " in founders";
            If[ ! MatchQ[imputetarget, "founders" | "founder"],
                labels[[1]] = labels[[1]] <> " and " <> ToString[AccountingForm[errorcount[[1]]]] <> "/" <> ToString[AccountingForm[errorcount[[2]]]] <> " in offspring";
                labels[[2]] = labels[[2]] <> " and " <> ToString[AccountingForm[nonimputecount[[2]]-nonimputecount[[1]]]] <> "/" <> ToString[AccountingForm[nonimputecount[[2]]]] <> " in offspring";
            ];
            Print[labels[[1]],".\n",labels[[2]],"."];
            Print["Output files:\n\t", res//TableForm];
            Print["Done! Finished date =",DateString[], ". \tTime elapsed in total = ", Round[SessionTime[] - starttime,0.1], " Seconds."];
        ];
        res
    ]

(*In truemagicsnp, non-imputing/observed genotypes are set to 0*)
calImputeAccuracy[obsmagicsnp_,truemagicsnp_,estmagicsnp_, excludehalfimpute_:True] :=
    Module[ {nfounder,inputgeno = obsmagicsnp, estgeno = estmagicsnp,truegeno = truemagicsnp, estsub, frac, truesub, accu, fracfounder, 
      accufounder, fracoff, accuoff, res,index,rule,diff,inputsub,temp,geno,pos,ch},
        nfounder = obsmagicsnp[[1, 2]];
        inputgeno[[5 ;;, 2 ;;]] = Replace[obsmagicsnp[[5 ;;, 2 ;;]], {"0|0" | "N" | "NN" -> 1,"1N" | "2N" | "N1" | "N2" -> 0.5, _ -> 0}, {2}];
        truegeno[[5 ;;, 2 ;;]] *= Sign[inputgeno[[5 ;;, 2 ;;]]];
        estgeno[[5 ;;, 2 ;;]] *= Sign[inputgeno[[5 ;;, 2 ;;]]];
        (*calculate accurancy for unphased genotypes*)
        rule = {"21" -> "12", 21 -> 12, "N1" -> "1N", "N2" -> "2N"};
        estgeno[[5 ;;, 2 ;;]] = estgeno[[5 ;;, 2 ;;]] /. rule;
        truegeno[[5 ;;, 2 ;;]] = truegeno[[5 ;;, 2 ;;]] /. rule;
        If[ excludehalfimpute,
            estgeno[[5 ;;, 2 ;;]] = estgeno[[5 ;;, 2 ;;]] /. {"N1" | "N2" | "1N" | "2N" -> "NN"};
        ];
        res = Table[
            {inputsub, estsub,truesub} = Table[
                temp = #[[All, 2 ;;]] & /@SplitBy[Transpose[geno[[Join[{3}, index], 2 ;;]]], First];
                DeleteCases[Flatten[#], 0] & /@ temp, {geno, {inputgeno, estgeno,truegeno}}];
            frac = Table[
               pos = Flatten[Position[estsub[[ch]], "N" | "NN"]];
               frac = Total[inputsub[[ch, pos]]];
               pos = Flatten[Position[estsub[[ch]], "N1" | "N2" | "1N" | "2N"]];
               frac += Total[inputsub[[ch, pos]] /. {1 -> 0.5}];
               {frac, Total[inputsub[[ch]]]}, {ch, Length[estsub]}];
            frac[[All, 1]] = frac[[All, 2]] - frac[[All, 1]];
            accu = Table[
                diff = Transpose[{inputsub[[ch]], truesub[[ch]] - estsub[[ch]]}];
                pos = Flatten[Position[diff[[All, 2]], 0]];
                accu = Total[diff[[pos, 1]]];
                pos = Flatten[Position[diff[[All,2]], (11 - "1N") | (12 - "1N") | (22 - "2N") | (12 - "2N")]];
                accu += Total[diff[[pos, 1]] - 0.5];
                accu, {ch, Length[inputsub]}];
            accu = Transpose[{accu, frac[[All, 1]]}];
            {frac, accu}, {index, {Range[5, nfounder + 4],Range[nfounder + 5, Length[estgeno]]}}];
        {fracfounder, accufounder, fracoff, accuoff} = Flatten[res, 1];
        {fracfounder, accufounder, fracoff, accuoff}
    ]
    
(*calErrorDetection returns;
ntruedetecting=#true detected errors;
ntruecorrecting=#true corrected errors; 
ndetecting=#differences between observed and estimated genotypes;
nobserror=#differences between observed and true genotypes;
nobsgeno=#observed genotypes*)
calErrorDetection[obsmagicsnp_, truemagicsnp_, estmagicsnp_] :=
    Module[ {nfounder, ss, inputsnp, truesnp, estsnp, obserrorpos, 
      correctingpos, ntruedetecting = 0, ntruecorrecting = 0, rule,pos,nobserror, ndetecting, nobsgeno},
        nfounder = obsmagicsnp[[1, 2]];
        ss = nfounder + 5 ;; All;
        {inputsnp, truesnp, estsnp} = {obsmagicsnp, truemagicsnp,estmagicsnp};
        rule = {11 -> "11", 22 -> "22", 12 | 21 | "21" -> "12", "N1" -> "1N", "N2" -> "2N", 1 -> "1", 2 -> "2"};
        {truesnp[[ss, 2 ;;]], estsnp[[ss, 2 ;;]], inputsnp[[ss, 2 ;;]]} = 
          Replace[#, rule, {2}] & /@ {truesnp[[ss, 2 ;;]], estsnp[[ss, 2 ;;]],inputsnp[[ss, 2 ;;]]};
        pos = Position[Replace[inputsnp[[ss, 2 ;;]], {"N" | "NN" -> 0, _ -> 1}, {2}], 1];
        {truesnp, estsnp, inputsnp} = Extract[#, pos] & /@ {truesnp[[ss, 2 ;;]], estsnp[[ss, 2 ;;]],inputsnp[[ss, 2 ;;]]};
        rule = {"N" -> {"1", "2"}, "1N" -> {"11", "12"}, "2N" -> {"12", "22"},"NN" -> {"11", "12", "22"}, x_ :> {x}};
        inputsnp = Replace[inputsnp, rule, {1}];
        estsnp = Replace[estsnp, rule, {1}];
        nobsgeno = Length[inputsnp];
        obserrorpos = Flatten[Position[MapThread[MemberQ[#1, #2] &, {inputsnp, truesnp}],False]];
        nobserror = Length[obserrorpos];
        correctingpos = Flatten[Position[MapThread[Intersection[#1, #2] === {} &, {inputsnp, estsnp}],True]];
        ndetecting = Length[correctingpos];
        If[ ndetecting >= 1,
            {inputsnp, estsnp, truesnp} = {inputsnp, estsnp, truesnp}[[All, correctingpos]];
            ntruedetecting = Count[MapThread[MemberQ[#1, #2] &, {inputsnp, truesnp}], False];
            ntruecorrecting = Count[MapThread[MemberQ[#1, #2] &, {estsnp, truesnp}], True];
        ];
        {ntruedetecting, ntruecorrecting, ndetecting, nobserror, nobsgeno}
    ]

calErrorDetection[obsmagicsnp_, truemagicsnp_, estmagicprob_, detectingthreshold_] :=
    Module[ {estmagicsnp, nfounder, ss, bool, pos, ls},
        estmagicsnp = estmagicprob;
        nfounder = obsmagicsnp[[1, 2]];
        ss = nfounder + 5 ;; All;
        bool = Replace[obsmagicsnp[[ss,2 ;;]], {"N" | "NN" -> 0, _ -> 1}, {2}];
        estmagicsnp[[ss, 2 ;;]] *= bool;
        pos = Position[estmagicsnp[[ss, 2 ;;]], _String];
        If[ pos=!={},
            ls = ToExpression[StringSplit[Extract[estmagicsnp[[ss, 2 ;;]], pos], "|"]];
            ls = maxIndexPair[ls];
            ls = Boole[Thread[ls[[All, 1]] > detectingthreshold]] ls[[All, 2, 1]];
            ls = ls /. Thread[Range[0, 4] -> {"NN","11", "12", "21", "22"}];
            estmagicsnp[[ss, 2 ;;]] = ReplacePart[estmagicsnp[[ss, 2 ;;]], Thread[pos -> ls]];
        ];
        calErrorDetection[obsmagicsnp, truemagicsnp, estmagicsnp]
    ]
      
      
calphasingaccu[pairhaplo_] :=
    Module[ {ls, phaseerror, pair, count, switcherror, nhetero,i},
        If[ pairhaplo === {},
            {{0, 0}, {0, 0}},
            ls = {#[[1]], Length[#]} & /@ Split[pairhaplo];
            nhetero = Total[ls[[All, 2]]];
            phaseerror = {Min[Sign[Abs[ls[[All, 1, 1]] - #]].ls[[All, 2]] & /@ {ls[[All, 1,2]], IntegerReverse[ls[[All, 1, 2]]]}], nhetero-1};
            pair = ls[[All, 1]];
            count = If[ pair[[1, 1]] == pair[[1, 2]],
                        0,
                        -1
                    ];
            Do[
             If[ pair[[i, 1]] =!= pair[[i, 2]],
                 pair[[i ;;, 2]] = IntegerReverse[pair[[i ;;, 2]]];
                 count++
             ], {i, Length[pair]}];
            If[ pair[[All, 2]] != pair[[All, 1]],
                Print["Wrong in switch error calculation!"];
                Abort[]
            ];
            switcherror = {count, nhetero - 1};
            {#[[2]]-#[[1]],#[[2]]}&/@{switcherror, phaseerror}
        ]
    ]

calPhasingAccuracy[truemagicsnp_, estmagicsnp_] :=
    Module[ {nfounder, i,ch,ss, chrtruesnp, chrestsnp, bool, 
      truegeno, estgeno, pos, nchr},
        nfounder = estmagicsnp[[1, 2]];
        nchr = Length[Split[estmagicsnp[[3, 2 ;;]]]];
        Table[Sum[
           {chrtruesnp, chrestsnp} = getsubMagicSNP[#, {ch}, All] & /@ {truemagicsnp, estmagicsnp};
           bool = Replace[chrtruesnp[[ss, 2 ;;]], {12 | 21 | "12" | "21" ->1, _ -> 0}, {2}];
           bool *= Replace[chrestsnp[[ss, 2 ;;]], {12 | 21 | "12" | "21" -> 1, _ ->0}, {2}];
           {truegeno, estgeno} = ToExpression[#*bool] & /@ {chrtruesnp[[ss, 2 ;;]],chrestsnp[[ss, 2 ;;]]};
           pos = Flatten[Position[#, 1]] & /@ bool;
           Sum[calphasingaccu[Transpose[{truegeno[[i, pos[[i]]]],estgeno[[i, pos[[i]]]]}]], {i, Length[truegeno]}], {ch,nchr}], {ss, {5 ;; nfounder + 4, nfounder + 5 ;; All}}]
    ]
    
calHeteroAccuracy[truemagicsnp_, estmagicsnp_] :=
    Module[ {nfounder, truegeno, estgeno, rule, postrue, posest, 
      estheteroacc, falsenegative},
        nfounder = estmagicsnp[[1, 2]];
        {truegeno, estgeno} = {truemagicsnp, estmagicsnp};
        (*calculate accurancy for unphased heterozygote*)
        rule = {"21" | 21 | "12" -> 12};
        estgeno[[nfounder + 5 ;;, 2 ;;]] = estgeno[[nfounder + 5 ;;, 2 ;;]] /. rule;
        truegeno[[nfounder + 5 ;;, 2 ;;]] = truegeno[[nfounder + 5 ;;, 2 ;;]] /. rule;
        postrue = Position[truegeno[[nfounder + 5 ;;, 2 ;;]], 12];
        posest = Position[estgeno[[nfounder + 5 ;;, 2 ;;]], 12];
        estheteroacc = {Length[Intersection[posest, postrue]], Length[posest]};
        falsenegative = {Length[Complement[postrue, posest]], Length[postrue]};
        {estheteroacc, falsenegative}
    ]
     
genoErrorPattern[obsmagicsnp_, estmagicsnp_,truemagicsnp_] :=
    Module[ {obs, est, true, ls, set, esterror, obserror, pos, falseneg, 
      truedetect, truecorrect, falseimpute, falsepos, lab, rule, res},
     (*"TrueCorrect": genotype errors that are changed correctly;
     "TrueDetect": genotype errors that are changed wrongly;
     "FalseNegative": genotype errors that are not detected; 
     "FalsePositive": non genotype errors that are wrongly corrected;
     "FalseImpute": wrongly imputed genotypes*)
        {obs, est, true} = #[[5 ;;, 2 ;;]] & /@ {obsmagicsnp, estmagicsnp,truemagicsnp};
        {obs, est, true} = Map[ToString, {obs, est, true}, {3}] /. {"21" -> "12","N1" -> "1N", "N2" -> "2N"};
        ls = MapThread[List, {obs, est, true}, 2];
        set = Join[Thread[{"11", {"12", "22"}}],Thread[{"12", {"11", "22"}}], Thread[{"22", {"11", "12"}}], 
          Thread[{"1N", {"22"}}], Thread[{"2N", {"11"}}], {{"1", "2"}, {"2", "1"}}];
        esterror = Position[Map[MemberQ[set, #] &, ls[[All, All, 2 ;; 3]], {2}], True];
        obserror = Position[Map[MemberQ[set, #] &, ls[[All, All, {1, 3}]], {2}], True];
        pos = Intersection[obserror, esterror];
        falseneg = Pick[pos, MapThread[SameQ, {Extract[obs, pos], Extract[est, pos]}]];
        truedetect = Complement[pos, falseneg];
        truecorrect = Complement[obserror, esterror];
        pos = Complement[esterror, obserror];
        falseimpute = Pick[pos, Extract[obs, pos], "N"|"NN" | "1N" | "2N"];
        falsepos = Complement[pos, falseimpute];
        pos = {truecorrect, truedetect, falseneg, falsepos, falseimpute};
        lab = {"TrueCorrect", "TrueDetect", "FalseNegative","FalsePositive", "FalseImpute"};
        rule = Flatten[MapThread[Thread[#1 -> #2] &, {pos, lab}]];
        res = SparseArray[rule, Dimensions[obs]];
        res = Join[obsmagicsnp[[5 ;;, {1}]], res, 2];
        Join[obsmagicsnp[[;; 4]], res]
    ]

genoErrorPattern[obsmagicsnp_, estmagicsnp_] :=
    Module[ {obs, est, missgeno, posmiss, nonimpute, impute, obsgeno, 
      posobs, correction, pos, lab, rule, res},
     (*"Imputation": missing genotypes 1N/2N/NN/N are imputed;
     "NonImputation": missing genotypes 1N/2N/NN/N are not imputed;
     "Correction": obsvered genotypes 11/12/22/1/2 are changed;*)
        {obs, est} = #[[5 ;;, 2 ;;]] & /@ {obsmagicsnp, estmagicsnp};
        {obs, est} = Map[ToString, {obs, est}, {3}] /. {"21" -> "12", "N1" -> "1N","N2" -> "2N"};
        missgeno = {"1N", "2N", "NN", "N"};
        posmiss = Position[obs, Alternatives @@ missgeno];
        nonimpute = Pick[posmiss,MapThread[SameQ[#1, #2] || MatchQ[#2, "NN" | "N"] &, {Extract[obs, posmiss], Extract[est, posmiss]}]];
        impute = Complement[posmiss, nonimpute];
        obsgeno = {"11", "12", "22", "1", "2"};
        posobs = Position[obs, Alternatives @@ obsgeno];
        correction = Pick[posobs,MapThread[UnsameQ, {Extract[obs, posobs], Extract[est, posobs]}]];
        pos = {nonimpute, impute, correction};
        lab = {"NonImputed", "Imputed", "Corrected"};
        rule = Flatten[MapThread[Thread[#1 -> #2] &, {pos, lab}]];
        res = SparseArray[rule, Dimensions[obs]];
        res = Join[obsmagicsnp[[5 ;;, {1}]], res, 2];
        Join[obsmagicsnp[[;; 4]], res]
    ]

Options[plotErrorPatternGUI] = Append[Options[MatrixPlot],linkageGroupSet->All]
plotErrorPatternGUI[obsmagicsnp_?(StringQ[#]||ListQ[#]&), estmagicsnp_?(StringQ[#]||ListQ[#]&), truemagicsnp_?(StringQ[#]||ListQ[#]&), opts : OptionsPattern[]] :=
    Module[ {chrsubset,obssnp,estsnp, truesnp,tc, td, fn, fp, fi,f, statuses, nfounder, noffspring, lefttick,
       righttick, label, panel, sliders,colorrules,count, descrip,table,xylab},
        {obssnp,estsnp, truesnp} = {obsmagicsnp,estmagicsnp, truemagicsnp};
        If[ StringQ[obssnp],
            If[ !FileExistsQ[obssnp],
                Print["File ", obssnp," does not exist!"];
                Return[$Failed]
            ];
            Print[Directory[]];
            obssnp = Import[obssnp,"CSV", Path -> Directory[]];
        ];
        If[ StringQ[estsnp],
            If[ !FileExistsQ[estsnp],
                Print["File ", estsnp," does not exist!"];
                Return[$Failed]
            ];
            estsnp = Import[estsnp,"CSV", Path -> Directory[]];
        ];
        If[ StringQ[truesnp],
            If[ !FileExistsQ[truesnp],
                Print["File ", truesnp," does not exist!"];
                Return[$Failed]
            ];
            truesnp = Import[truesnp,"CSV", Path -> Directory[]];
        ];
        chrsubset = OptionValue[linkageGroupSet];
        {obssnp,estsnp, truesnp} = getsubMagicSNP[#, chrsubset,All]&/@{obssnp,estsnp, truesnp};
        statuses = genoErrorPattern[obssnp,estsnp, truesnp];
        nfounder = statuses[[1, 2]];
        noffspring = Length[statuses] - (nfounder + 4);
        {lefttick, righttick} = Table[Transpose[{Range[nfounder],f[ToString[#], IntegerLength[noffspring]] & /@ 
             Range[nfounder]}], {f, {StringPadLeft, StringPadRight}}];
        label = {"True Correct", "True Detect", "False Negative","False Positive", "False Impute"};
        xylab = {{{None, None}, {"Founder index","Founder index"}},{{"SNP index", None}, {"Offspring index","Offspring index"}}};
        panel = Panel[Column[{
            sliders = ColorSlider[Dynamic[#1],AppearanceElements -> "SwatchSpectrum"] & /@ {tc, td, fn,fp, fi};
            {tc, td, fn, fp, fi} = {Cyan,Green, Blue, Pink, Red};
            Grid[{Append[label, ""],Append[sliders,Button["Reset colors", {tc, td, fn, fp, fi} = {Cyan,Green, Blue, Pink, Red}]]}],
            Dynamic[
             colorrules = Thread[{0, "TrueCorrect", "TrueDetect", "FalseNegative","FalsePositive", "FalseImpute"} -> {White, tc, td, fn, fp,fi}];             
             Column[{MatrixPlot[statuses[[5 ;; nfounder + 4, 2 ;;]], Sequence@@FilterRules[{opts}, Options[MatrixPlot]], FrameLabel -> xylab[[1]], 
                FrameTicks -> {{lefttick, righttick}, {Automatic, Automatic}}, ColorRules -> colorrules,MaxPlotPoints -> Infinity, ImageSize -> 1000],
               MatrixPlot[statuses[[nfounder + 5 ;;, 2 ;;]], Sequence@@FilterRules[{opts}, Options[MatrixPlot]], FrameLabel -> xylab[[2]], 
                   ColorRules -> colorrules,MaxPlotPoints -> Infinity, ImageSize -> 1000]}]
             ]}]
          ];
        descrip = {"Number of genotype errors that are changed correctly", 
          "Number of genotype errors that are changed wrongly", 
          "Number of genotype errors that are not detected", 
          "Number of correctly observed genotypes that are changed", 
          "Number of wrongly imputed genotypes"};
        count = Table[Count[Flatten[statuses[[f, 2 ;;]]], #] & /@ {"TrueCorrect","TrueDetect", "FalseNegative", "FalsePositive","FalseImpute"}, {f, {5 ;; nfounder + 4, nfounder + 5 ;; All}}];
        table = TableForm[Transpose[Append[count, descrip]],TableHeadings -> {label, {"Founder", "Offspring", "Description"}}];
        TabView[{"Summary" -> table, "Pattern Plot" -> panel, "Pattern Table"->TableForm[statuses],"Observed Genotypes"->TableForm[obssnp],
             "Estimated Genotypes"->TableForm[estsnp],"True Genotypes"->TableForm[truesnp]}, 2,ImageSize->Automatic]
    ]


plotErrorPatternGUI[obsmagicsnp_?(StringQ[#]||ListQ[#]&), estmagicsnp_?(StringQ[#]||ListQ[#]&), opts : OptionsPattern[]] :=
    Module[ {chrsubset,obssnp,estsnp,nim, im, cc,f, statuses, nfounder, noffspring, lefttick,
       righttick, label, panel, sliders,colorrules,count, descrip,table,xylab},
        {obssnp,estsnp} = {obsmagicsnp,estmagicsnp};
        If[ StringQ[obssnp],
            If[ !FileExistsQ[obssnp],
                Print["File ", obssnp," does not exist!"];
                Return[$Failed]
            ];
            obssnp = Import[obssnp,"CSV", Path -> Directory[]];
        ];
        If[ StringQ[estsnp],
            If[ !FileExistsQ[estsnp],
                Print["File ", estsnp," does not exist!"];
                Return[$Failed]
            ];
            estsnp = Import[estsnp,"CSV", Path -> Directory[]];
        ];
        chrsubset = OptionValue[linkageGroupSet];
        (*Print[chrsubset];
        Abort[];*)
        {obssnp,estsnp} = getsubMagicSNP[#, chrsubset,All]&/@{obssnp,estsnp};
        statuses = genoErrorPattern[obssnp, estsnp];
        nfounder = statuses[[1, 2]];
        noffspring = Length[statuses] - (nfounder + 4);
        {lefttick, righttick} = Table[Transpose[{Range[nfounder],f[ToString[#], IntegerLength[noffspring]] & /@Range[nfounder]}], {f, {StringPadLeft, StringPadRight}}];
        label = {"NonImputed", "Imputed", "Corrected"};
        xylab = {{{None, None}, {"Founder index", "Founder index"}}, {{"SNP index", None}, {"Offspring index","Offspring index"}}};
        panel = Panel[Column[{
             sliders = ColorSlider[Dynamic[#1],AppearanceElements -> "SwatchSpectrum"] & /@ {nim, im, cc};
             {nim, im, cc} = {Red, RGBColor[0, 1, 1], Blue};
             Grid[{Append[label, ""], Append[sliders,Button["Reset colors", {nim, im, cc} = {Red, RGBColor[0, 1, 1],Blue}]]}],
             Dynamic[
              colorrules = Thread[{0, "NonImputed", "Imputed", "Corrected"} -> {White,nim, im, cc}];
              Column[{MatrixPlot[statuses[[5 ;; nfounder + 4, 2 ;;]], Sequence@@FilterRules[{opts}, Options[MatrixPlot]], 
                 FrameLabel -> xylab[[1]], FrameTicks -> {{lefttick, righttick}, {Automatic,Automatic}}, ColorRules -> colorrules, 
                 MaxPlotPoints -> Infinity, ImageSize -> 1000], 
                MatrixPlot[statuses[[nfounder + 5 ;;, 2 ;;]], Sequence@@FilterRules[{opts}, Options[MatrixPlot]], 
                	FrameLabel -> xylab[[2]], ColorRules -> colorrules,MaxPlotPoints -> Infinity, ImageSize -> 1000]}]
              ]}]
           ];
        descrip = {"Number of missing genotypes that are NOT imputed", 
           "Number of missing genotype that are imputed", 
           "Number of observed genotype that are changed"};
        count = Table[Count[Flatten[statuses[[f, 2 ;;]]], #] & /@ {"NonImputed","Imputed", "Corrected"}, {f, {5 ;; nfounder + 4,nfounder + 5 ;; All}}];
        table = TableForm[Transpose[Append[count, descrip]], TableHeadings -> {label, {"Founder", "Offspring", "Description"}}];
        TabView[{"Summary" -> table, "Pattern Plot" -> panel, "Pattern Table" -> TableForm[statuses], 
          "Observed Genotypes" -> TableForm[obssnp], "Estimated Genotypes" -> TableForm[estsnp]}, 2,ImageSize -> Automatic]
    ]
    
        
End[]

EndPackage[]


