(* Mathematica Package *)

(* Created by the Wolfram Workbench 30-Nov-2014 *)

BeginPackage["MagicReconstruct`MagicModel`",{"MagicDefinition`","MagicOrigin`","MagicOriginXY`","PedigreeOrigin`","MagicSimulate`","MagicGeneDropping`"}]
(* Exported symbols added here with SymbolName::usage *) 
    
Unprotect @@ Names["MagicReconstruct`MagicModel`*"];
ClearAll @@ Names["MagicReconstruct`MagicModel`*"];

sampleDiscretePriorProcess::usage = "sampleDiscretePriorProcess  "

relabelMarkovProcess::usage = "relabelMarkovProcess  "

sampleStartprob::usage = "sampleStartprob  "

sampleTransition::usage = "sampleTransition  "

pedPriorInfor::usage = "pedPriorInfor  "

sampleContinuedPriorProcess::usage = "sampleContinuedPriorProcess  "

toDiscreteMarkovProcess::usage = "toDiscreteMarkovProcess  "

toStartTranProb::usage = "toStartTranProb  "

magicContinuedPriorProcess::usage = "magicContinuedPriorProcess  "

magicMarkovProb::usage = "magicMarkovProb  "

magicMarkovProb2::usage = "magicMarkovProb2  "

sampleCode::usage = "sampleCode  "

(*pedContinuedPriorProcess::usage = "pedContinuedPriorProcess  "*)

Begin["`Private`"]
(* Implementation of the package *)
   
(*revision note1: transitionProbold is replaced by transitionProb, to given consistent results obtained from the rate matrix of CTMM*)
(*revision note2:
1. The magicMarkovProb and magicMarkovProbMale are calculated using  transitionProb, maleXTransitionProb, 
indepTransitionProb, and depTransitionProb, where the inter-marker transition probability are analytically dervied
from the rate matrix of CTMM.
2. For the magicMarkovProb2 and magicMarkovProbMale2, the inter-marker transition probability are calcuated numerically
from the rate matrix of CTMM. Thus the (2) functions gain generaliblity at a little cost of speed.  
3. magicMarkovProb2 and magicMarkovProbMale2 are not used. 
*)

transitionProbold = Function[{deltt, nFounder, mapR},
   Module[ {aa},
       aa = Exp[-mapR deltt] IdentityMatrix[nFounder];
       aa += (1 - IdentityMatrix[nFounder]) (1-Exp[-mapR deltt])/(nFounder - 1);
       aa
   ]
   ]

(*aa = -R IdentityMatrix[nFounder];
aa += (1 - IdentityMatrix[nFounder]) R/(nFounder - 1);
Simplify[MatrixExp[aa t] == transitionProb[t, nFounder, R]]*)
transitionProb = Function[{deltt, nFounder, mapR},
   Module[ {aa},
       aa = IdentityMatrix[nFounder] (nFounder - 1)/nFounder - (1 - IdentityMatrix[nFounder]) /nFounder;
       aa* Exp[-mapR deltt nFounder/(nFounder - 1)] + 1/nFounder
   ]
   ]
      
(*mapR =R^m for the maternally derived X chromosome of a male; mapR =(R^m+R^p)/2 for autsomes*)
maleXTransitionProb = transitionProb

depTransitionProb = transitionProb

(*effective generation= {R^m,R^p} 
for maternally and paternally derived chromosomes.
R^m = R^p for autosomes*)
indepTransitionProb = Function[{deltt, nFounder, mapR},
  Module[ {mp,aals},
      aals = Table[transitionProb[deltt, nFounder, mp], {mp, mapR}];
      KroneckerProduct @@ aals
  ]
  ]
   
(*for both autosomes and female XX chromosomes*)
magicMarkovProb[deltd_, model_, nFounder_, juncdist_] :=
    Module[ {startProb, rr, transitionRate, tranProb,Rm,Rp,inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp} = juncdist;
        Rm = 2 j1222mp + j1122mp + j1232mp;
        Rp = 2 j1211mp + j1122mp +j1213mp;
        Switch[model,
          "jointModel",
          startProb = magicStationaryProbXY[nFounder, inbred];
          rr = magicBasicRateXY[nFounder, inbred, {j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}];
          transitionRate = magicRateMatrixXY[nFounder, rr];
          tranProb = Map[MatrixExp[transitionRate #] &, deltd, {2}],                 
          "indepModel",
          startProb = Table[1./nFounder^2, {nFounder^2}];
          tranProb = Map[indepTransitionProb[#, nFounder, {Rm,Rp}]&, deltd,{2}],
          "depModel",
          startProb = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
          startProb = startProb /. {{x_, x_} :> 1./nFounder, {_, _} -> 0};
          startProb = Table[1./nFounder, {nFounder}];
          (*tranProb = Table[depTransitionProb[deltd[[i,j]], nFounder, Mean[{Rm,Rp}]], {i,Length[deltd]},{j,Length[deltd[[i]]]}];*)
          tranProb = Map[depTransitionProb[#, nFounder, Mean[{Rm,Rp}]] &, deltd,{2}];
        ];
        {startProb, tranProb}
    ]  


(*only model like "depModel" is meaningfull for the maternally derived X chromsome of a male*)
magicMarkovProbMale[deltd_,  nFounder_,juncdist_] :=
    Module[ {Rm,startProb, tranProb,inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp} = juncdist;
        Rm = 2 j1222mp + j1122mp + j1232mp;
        startProb = Table[1./nFounder, {nFounder}];
        tranProb = Map[maleXTransitionProb[#, nFounder, Rm] &, deltd,{2}];
        {startProb, tranProb}
    ]  

(*Table[Total[Flatten[magicMarkovProb[deltd, model, nFounder, juncdist[[1]]] - 
    magicMarkovProb2[deltd, model, nFounder,juncdist[[1]]]]],
     {model, {"jointModel", "depModel","indepModel"}}]*)
(*for both autosomes and female XX chromosomes*)
magicMarkovRate[model_, nFounder_, juncdist_] :=
    Module[ {inbred, j1122mp, j1211mp, j1213mp, j1222mp, j1232mp, Rm, Rp, 
      startProb, rr, transitionRate, mapR, aa},
        {inbred, j1122mp, j1211mp, j1213mp, j1222mp, j1232mp} = juncdist;
        Rm = 2 j1222mp + j1122mp + j1232mp;
        Rp = 2 j1211mp + j1122mp + j1213mp;
        Switch[model,
         "jointModel",
         startProb = magicStationaryProbXY[nFounder, inbred];
         rr = magicBasicRateXY[nFounder, inbred, {j1122mp, j1211mp, j1213mp, j1222mp, j1232mp}];
         transitionRate = magicRateMatrixXY[nFounder, rr],
         "indepModel",
         startProb = Table[1./nFounder^2, {nFounder^2}];
         transitionRate = Table[
             aa = -mapR IdentityMatrix[nFounder];
             aa += (1 - IdentityMatrix[nFounder]) mapR/(nFounder - 1);
             aa,{mapR,{Rm,Rp}}];
         (*Kroneckersum[a,b] =!= ArrayFlatten[outer[Plus, a,b]]*)
         transitionRate = kroneckerMatrixSum[transitionRate[[1]],transitionRate[[2]]],
         "depModel",
         (*startProb = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1];
         startProb = startProb /. {{x_, x_} :> 1./nFounder, {_, _} -> 0};*)
         startProb = Table[1./nFounder, {nFounder}];
         mapR = Mean[{Rm, Rp}];
         transitionRate = -mapR IdentityMatrix[nFounder];
         transitionRate += (1 - IdentityMatrix[nFounder]) mapR/(nFounder - 1);
         ];
        {startProb, transitionRate}
    ] 

(*only model like "depModel" is meaningfull for the maternally derived X chromsome of a male*)
magicMarkovRateMale[nFounder_,juncdist_] :=
    Module[ {Rm,startProb, transitionRate,inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp},
        {inbred,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp} = juncdist;
        Rm = 2 j1222mp + j1122mp + j1232mp;
        startProb = Table[1./nFounder, {nFounder}];
        transitionRate = -Rm IdentityMatrix[nFounder];
        transitionRate += (1 - IdentityMatrix[nFounder]) Rm/(nFounder - 1);
        {startProb, transitionRate}
    ]  
    
magicMarkovProb2[deltd_, model_, nFounder_, juncdist_] :=
    Module[ {startProb, transitionRate,tranProb},
        {startProb, transitionRate} = magicMarkovRate[model, nFounder, juncdist];
        tranProb = Map[MatrixExp[transitionRate #] &, deltd, {2}];
        {startProb, tranProb}
    ]
    
magicMarkovProbMale2[deltd_, nFounder_, juncdist_] :=
    Module[ {startProb, transitionRate,tranProb},
        {startProb, transitionRate} = magicMarkovRateMale[nFounder, juncdist];
        tranProb = Map[MatrixExp[transitionRate #] &, deltd, {2}];
        {startProb, tranProb}
    ]

(*
1) Random mating mapping population with non-overlappign generations, by default, 
    popDeisgn = {magicOrigPrior[nFounder,popDesign],magicOrigPriorXY[nFounder,popDesign]}
2) Diversity Outbred (DO) populations:
    popDesign = {magicOrigPriorDO[nPower, preCCfreq, crossPopSize, gCross, crossScheme],
                 magicOrigPriorDOXY[nPower, preCCfreq, crossPopSize, gCross, crossScheme]};    
3) Fixed pedigree with non-overlapping generations: 
    popDeisgn = pedOrigPrior[popPed, chrLength,interferStrength,isObligate,isOogamy,sampleSize];    
*)
magicJuncdist[nFounder_,popDesign_,isfounderinbred_,isincludeXX_] :=
    Module[ {juncdist},
        (*juncdist={{f,j1122,j1211,j1213,j1222,j1232},{fmp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}}*)
        If[ VectorQ[popDesign, StringQ],
            juncdist = {magicOrigPrior[nFounder,popDesign,isfounderinbred]};
            If[ isincludeXX,
                AppendTo[juncdist,magicOrigPriorXY[nFounder,popDesign,isfounderinbred]]
            ],
            juncdist = popDesign;
        ];
        juncdist
    ]
 
 magicContinuedPriorProcess[nFounder_,popDesign_,isfounderinbred_,model_,posA_,posX_] :=
     Module[ {nFGL,juncdist,startProbAutosome, tranRateAutosome,startProbFemale, tranRateFemale,startProbMale, tranRateMale},
         juncdist = magicJuncdist[nFounder,popDesign,isfounderinbred,posX=!={}];
         nFGL = If[ isfounderinbred,
                    nFounder,
                    2 nFounder
                ];
         If[ posA=!={},
             If[ !(VectorQ[First[juncdist], NonNegative] && Length[First[juncdist]] == 6&&juncdist[[1,1]]<=1),
                 Print["magicContinuedPriorProcess: wrong popDesign for autosomes!","popDesign = ",popDesign];
                 Abort[]
             ];
             {startProbAutosome, tranRateAutosome} = magicMarkovRate[model, nFGL, First[juncdist]];
         ];
         If[ posX=!={},
             If[ !(VectorQ[Last[juncdist], NonNegative] && Length[Last[juncdist]] == 6&&juncdist[[-1,1]]<=1),
                 Print["Wrong popDesign for X chromosomes!","popDesign = ",popDesign];
                 Abort[]
             ];
             {startProbFemale, tranRateFemale} = magicMarkovRate[model, nFGL, Last[juncdist]];
             {startProbMale, tranRateMale} = magicMarkovRateMale[nFGL, Last[juncdist]];
         ];
         {juncdist,startProbAutosome, tranRateAutosome,startProbFemale, tranRateFemale,startProbMale, tranRateMale}
     ] 
    
(*magicDiscretePriorProcess[nFounder_,popDesign_,isfounderinbred_,model_,posA_, posX_, deltd_] :=
    Module[ {juncdist,startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale,
        tranRateAutosome, tranRateFemale, tranRateMale},
        {juncdist,startProbAutosome, tranRateAutosome,startProbFemale, tranRateFemale,startProbMale, tranRateMale} = 
            magicContinuedPriorProcess[nFounder,popDesign,isfounderinbred,model,posA,posX];
        tranProbAutosome=Map[MatrixExp[tranRateAutosome #] &, deltd[[posA]], {2}];
        tranProbFemale=Map[MatrixExp[tranRateFemale #] &, deltd[[posX]], {2}];
        tranProbMale=Map[MatrixExp[tranRateMale #] &, deltd[[posX]], {2}];
        {juncdist,startProbAutosome, tranProbAutosome,startProbFemale, tranProbFemale,startProbMale, tranProbMale}
    ]  *)
        
pedPriorInfor[popDesign_, isfounderinbred_,isJointModel_] :=
    Module[ {pedigree, sampleInfor, indlist, nFounder, founderFGL, pedPriorList},
        {pedigree, sampleInfor} = Rest[splitPedigreeInfor[popDesign]];
        indlist = Union[sampleInfor[[2 ;;, 2]]];
        nFounder = Count[pedigree[[2 ;;, -1]], {0, 0}];
        If[ Complement[indlist, pedigree[[nFounder + 2 ;;, 2]]] =!= {},
            Print["pedPriorInfor: some sample IDs are found in the breeding pedigree!"];
            Abort[]
        ];    
        (*inbred founder parents*)
        founderFGL = If[ isfounderinbred,
                         Transpose[{Range[nFounder], Range[nFounder]}],
                         Partition[Range[2 nFounder], 2]
                     ];
        pedPriorList = pedAncestryPrior[pedigree, founderFGL, indlist, isJointModel];
        pedPriorList = Thread[pedPriorList[[2 ;;, 2]] -> N[pedPriorList[[2 ;;, 5;;]]]];
        {pedPriorList, sampleInfor}
    ]
  
pedpriorjointModel[prior_, posA_, posX_] :=
    Module[ {startProb, tranRate, iniAA, rateAA, iniSex, rateSex},
        startProb = tranRate = Table[0, {Length[posA] + Length[posX]}];
        If[ posA =!= {},
            {iniAA, rateAA} = First[prior];
            startProb[[posA]] = Table[iniAA, {Length[posA]}];
            tranRate[[posA]] = Table[rateAA,{Length[posA]}];
        ];
        If[ posX =!= {},
            (*For female, the markov process refers to the joint process along two XX chromosomes*)
            (*For male, the markov process refers to the process along maternally derived X chromosome.*)
            {iniSex, rateSex} = Last[prior];
            startProb[[posX]] = Table[iniSex, {Length[posX]}];
            tranRate[[posX]] = Table[rateSex,{Length[posX]}];
        ];
        {startProb, tranRate}
    ]

      
pedpriorindepModel[prior_, posA_, posX_,gender_] :=
    Module[ {startProb, tranRate, iniAm, rateAm, iniAp, rateAp, iniSm, 
      rateSm, iniSp, rateSp},
        startProb = tranRate = Table[0, {Length[posA] + Length[posX]}];
        If[ posA =!= {},
            {{iniAm, rateAm}, {iniAp, rateAp}} = First[prior];
            startProb[[posA]] = Table[Flatten[KroneckerProduct[iniAm, iniAp]], {Length[posA]}];
            tranRate[[posA]] = Table[kroneckerMatrixSum[rateAm,rateAp],{Length[posA]}];
        ];
        If[ posX =!= {},
            {{iniSm, rateSm}, {iniSp, rateSp}} = First[prior];
            If[ ToLowerCase[gender] === "male",
            	startProb[[posX]] = Table[iniSm, {Length[posX]}];
                tranRate[[posX]] = Table[rateSm,{Length[posX]}],
                startProb[[posX]] = Table[Flatten[KroneckerProduct[iniSm, iniSp]], {Length[posX]}];
                tranRate[[posX]] = Table[kroneckerMatrixSum[rateSm,rateSp],{Length[posX]}];                
            ];
        ];
        {startProb, tranRate}
    ]
  
pedpriordepModel[prior_, posA_, posX_, gender_] :=
    Module[ {startProb, tranRate, iniAm, rateAm, iniAp, rateAp, 
      iniSm, rateSm, iniSp, rateSp},
        startProb = tranRate = Table[0, {Length[posA] + Length[posX]}];
        If[ posA =!= {},
            {{iniAm, rateAm}, {iniAp, rateAp}} = First[prior];
            startProb[[posA]] = Table[(iniAm + iniAp)/2, {Length[posA]}];
            tranRate[[posA]] = Table[(rateAm + rateAp)/2,{Length[posA]}];
        ];
        If[ posX =!= {},
            {{iniSm, rateSm}, {iniSp, rateSp}} = Last[prior];
            If[ ToLowerCase[gender] === "male",
            	startProb[[posX]] = Table[iniSm, {Length[posX]}];
                tranRate[[posX]] = Table[rateSm,{Length[posX]}],
                startProb[[posX]] = Table[(iniSm + iniSp)/2, {Length[posX]}];
                tranRate[[posX]] = Table[(rateSm + rateSp)/2,{Length[posX]}];                
            ];
        ];
        {startProb, tranRate}
    ]
         
pedContinuedPriorProcess[prior_, model_, gender_, posA_, posX_] :=
    Switch[model,
         "jointModel",
         pedpriorjointModel[prior,posA, posX],
         "indepModel",
         pedpriorindepModel[prior, posA, posX,gender],
         "depModel",
         pedpriordepModel[prior, posA, posX,gender]
     ]
     
(*pedDiscretePriorProcess[prior_, model_, gender_, posA_, posX_, deltd_] :=
    Module[ {startProb, tranRate,tranProb},
        {startProb, tranRate} = pedContinuedPriorProcess[prior, model, gender, posA, posX];
        tranProb = Map[MatrixExp[tranRate #] &, deltd, {2}];
        {startProb, tranProb}
    ]*)       
     
(*sampleMarkovProcess is a list of rules for each group (according memeberID or gender): 
  MemberID/gender ->{startProb0, tranProb0}, not accounting for the ordering of funnelcode*)
(*startProb0 = a list of initial distribution for each linkage group; 
tranRate0 = a list of transition rate matrices for each linkage group*)
sampleContinuedPriorProcess[nFounder_, popDesign_, isfounderinbred_,model_, posA_, posX_, samplegender_, sampleid_] :=
    Module[ {pedigreekey = "Pedigree-Information", isPedigree, pedPriorList, juncdist, startprobls, tranratels,startprobAutosome, 
      tranrateAutosome, startprobFemale, tranrateFemale, startprobMale, tranrateMale, sampleLabel,
      sampleMarkovProcess, gender, genders,funnelcode,fglcode,haplocode,diplocode,startprob, tranrate,i},
        isPedigree = VectorQ[First[popDesign]] && popDesign[[1, 1]] == pedigreekey;
        If[ isPedigree,
            (*pedPriorList is a list of rule: MemberID -> junctionprior*)
            (*sampleLabel, columns ={sampleID, MemberID in breeding pedigree, funnelcode}*)
            {pedPriorList, sampleLabel} = pedPriorInfor[popDesign, isfounderinbred, model === "jointModel"];
            fglcode = sampleLabel[[2 ;;, 3]];
            If[ !isfounderinbred,
                fglcode = Transpose[sampleLabel[[2 ;;, 3]]];
                fglcode = Transpose[Riffle[2 fglcode - 1, 2 fglcode]];
            ];
            haplocode = toHaplocode[#]&/@fglcode;
            diplocode = toDiplocode[#]&/@fglcode;
            sampleLabel = Join[sampleLabel, Join[{{ "Haplocode", "Diplocode"}},Transpose[{haplocode, diplocode}]], 2];
            (*If[ sampleLabel[[2;;,1]]!=sampleid,
                Print["sampleContinuedPriorProcess: Inconsistent sampleid between popDesign and magicSNP!"];                
                Abort[]
            ];*)
            sampleMarkovProcess = Table[
              gender = samplegender[[Flatten[Position[sampleLabel[[2 ;;, 2]], pedPriorList[[i, 1]]]]]];
              If[ ! (SameQ @@ gender),
                  Print["sampleContinuedPriorProcess: Inconsistent genders!"];
                  Abort[]
              ];
              gender = First[gender];
              {startprobls, tranratels} = pedContinuedPriorProcess[pedPriorList[[i, 2]], model, gender, posA, posX];
              pedPriorList[[i, 1]] -> {startprobls, tranratels}, {i, Length[pedPriorList]}],            
            {juncdist, startprobAutosome, tranrateAutosome, startprobFemale, 
              tranrateFemale, startprobMale, tranrateMale} = magicContinuedPriorProcess[nFounder, popDesign, isfounderinbred, model, posA, posX];
            funnelcode = Range[nFounder];
            fglcode = 
              If[ isfounderinbred,
                  funnelcode,
                  Riffle[2 funnelcode - 1, 2 funnelcode]
              ];
            haplocode = toHaplocode[fglcode];
            diplocode = toDiplocode[fglcode];
            sampleLabel = Join[{{"SampleID", "PriorLabel(gender)", "Funnelcode", "Haplocode", "Diplocode"}}, 
                               Transpose[{sampleid, samplegender, Table[funnelcode, {Length[sampleid]}], 
                             Table[haplocode, {Length[sampleid]}],Table[diplocode, {Length[sampleid]}]}]];
            (*For random mating popDesign, priorlabels among samples are given by gender*)
            genders = Union[sampleLabel[[2 ;;, 2]]];
            sampleMarkovProcess = Table[
              startprob = tranrate = Table[0, {Length[posA] + Length[posX]}];
              If[ posA =!= {},
                  {startprob[[posA]],  tranrate[[posA]]} = Transpose[Table[{startprobAutosome, tranrateAutosome},{Length[posA]}]];
              ];
              Switch[ToLowerCase[gender],
                   "hermaphrodite"|"notapplicable",
                    0,
                   "female",
                   {startprob[[posX]], tranrate[[posX]]} = Transpose[Table[{startprobFemale, tranrateFemale},{Length[posX]}]],
                   "male",
                   {startprob[[posX]], tranrate[[posX]]} = Transpose[Table[{startprobMale, tranrateMale},{Length[posX]}]],
                   _, 
                   Print["sampleContinuedPriorProcess: wrong gender!"];
                   Abort[]
               ];
              gender -> {startprob, tranrate}, {gender, genders}];
        ];
        {sampleLabel, sampleMarkovProcess}
    ]
  
toDiscreteMarkovProcess[continuedMarkovProcess_, deltd_] :=
    Module[ {markovprocess = continuedMarkovProcess,ch,i},
        Do[markovprocess[[i, 2, 2]] = 
          Table[MatrixExp[markovprocess[[i, 2, 2, ch]] #] & /@deltd[[ch]], {ch, Length[deltd]}], {i, Length[markovprocess]}];
        markovprocess
    ]

toStartTranProb[samplelabel_, discretemarkovprocess_] := 
 Module[{res, startProb, tranProb, haplocode, diplocode,ind},
  res = Table[
    {haplocode, diplocode} = samplelabel[[ind + 1, 4 ;; 5]];
    {startProb, tranProb} = samplelabel[[ind + 1, 2]] /. discretemarkovprocess;
    If[! OrderedQ[haplocode],
     {startProb, tranProb} = relabelMarkovProcess[{startProb, tranProb}, haplocode,diplocode];
     ];
    {startProb, tranProb}, {ind, Length[samplelabel] - 1}];
  Transpose[res]
  ]    
    

sampleDiscretePriorProcess[nFounder_, popDesign_, isfounderinbred_,model_, posA_, posX_, samplegender_, sampleid_,deltd_] :=
    Module[ {sampleLabel, continuedMarkovProcess,discreteMarkovProcess},
        {sampleLabel, continuedMarkovProcess} = sampleContinuedPriorProcess[nFounder, popDesign, isfounderinbred, 
                       model, posA, posX, samplegender, sampleid];
        discreteMarkovProcess = toDiscreteMarkovProcess[continuedMarkovProcess, deltd];
        {sampleLabel, discreteMarkovProcess}
    ]
    
(*fglcode is a permutation of founder genome labels, natural integers starting from 1. toDiplocode returns reordering of nFgl^2 states*)
toDiplocode[fglcode_] :=
    Module[ {nfgl, states, newstates,pos},
        nfgl = Length[fglcode];
        (*origGenotype[nfgl][[1, 2]]*)
        states =  Flatten[Outer[List, Range[nfgl], Range[nfgl]], 1];
        newstates = states /. Thread[Range[nfgl] -> fglcode];
        (*pos=Flatten[FirstPosition[states,#]&/@newstates];
        states[[pos]]==newstates*)
        pos = SortBy[Transpose[{Range[Length[states]], newstates}], Last][[All, 1]];
        If[ newstates[[pos]]!= states,
            Print["toDiplocode: wrong relabeling by haplocode!"];
            Abort[]
        ];
        pos
    ]
    
toHaplocode[fglcode_] :=
    Module[ {nfgl,pos},
        nfgl = Length[fglcode];
        pos = SortBy[Transpose[{Range[nfgl], fglcode}], Last][[All,1]];
        If[ fglcode[[pos]]!= Range[nfgl],
            Print["toHaplocode: wrong relabeling by haplocode!"];
            Abort[]
        ];
        pos
    ]    
    
relabelMarkovProcess[markovprocess_, haplocode_, diplocode_] :=
    Module[ {startProb, transition,code,ishaplo,ch},
        {startProb, transition} = markovprocess;
        Transpose[Table[
            ishaplo = Length[startProb[[ch]]] === Length[haplocode];
            code = If[ ishaplo,
                       haplocode,
                       diplocode
                   ];
            {startProb[[ch,code]], Map[#[[code, code]] &, transition[[ch]], {Depth[transition[[ch]]] - 3}]},{ch,Length[startProb]}]]
    ]    

(*samplemarkovprocess refers a given chromosome*)
(*samplemarkovprocess: List of rule id->{startprob,tranrate}*)    
sampleStartprob[samplelabel_, samplemarkovprocess_] :=
    Module[ {haplocode, diplocode, code,prob, ind,ismaleX},
        Table[{haplocode, diplocode} = samplelabel[[ind + 1, 4;;5]];
              prob = First[samplelabel[[ind + 1, 2]] /. samplemarkovprocess];
              If[ ! OrderedQ[haplocode],
                  ismaleX = Length[prob] === Length[haplocode];
                  code = If[ ismaleX,
                             haplocode,
                             diplocode
                         ];
                  prob = prob[[code]]
              ];
              prob, {ind, Length[samplelabel] - 1}]
    ]    

(*samplemarkovprocess and ismaleX refer a given chromosome*)
(*samplemarkovprocess: List of rule id->{startprob,tranrate or tranprob for all intervals between neighbor markers}*)
sampleTransition[samplelabel_, samplemarkovprocess_, isdepModel_,ismaleX_,indices_] :=
    Module[ {process,code, tran},
    	code = sampleCode[samplelabel, isdepModel, ismaleX];
    	process = samplemarkovprocess;
    	process[[All, 2]] = process[[All, 2, 2, indices]];    	
    	tran = samplelabel[[2 ;;, 2]] /. process; 
    	If[Head[indices]===Integer, 
    		MapThread[#1[[#2, #2]] &, {tran, code}],
    		MapThread[#1[[All,#2, #2]] &, {tran, code}]
    	]
    ]
    
sampleCode[samplelabel_, isdepModel_,ismaleX_] :=
    Module[ {haplocode, diplocode, code, ind},
        Table[
         {haplocode, diplocode} = samplelabel[[ind + 1, 4;;5]];
         If[ ! OrderedQ[haplocode],
             code = If[ ismaleX[[ind]]||isdepModel,
                        haplocode,
                        diplocode
                    ],
             code =All
         ], {ind, Length[samplelabel] - 1}]
    ]    
    
sampleTransition[samplelabel_, samplemarkovprocess_,samplecode_, indices_] :=
    Module[ {process, tran},
    	process = samplemarkovprocess;
    	process[[All, 2]] = process[[All, 2, 2, indices]];
    	tran = samplelabel[[2 ;;, 2]] /. process; 
    	If[Head[indices]===Integer, 
    		MapThread[#1[[#2, #2]] &, {tran, samplecode}],
    		MapThread[#1[[All,#2, #2]] &, {tran, samplecode}]
    	]
    ]    

End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicReconstruct`MagicModel`*"];

EndPackage[]

