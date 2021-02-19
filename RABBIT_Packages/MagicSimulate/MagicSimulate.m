(* Mathematica Package *)

(* Created by the Wolfram Workbench Oct 2, 2013 *)

BeginPackage["MagicSimulate`",{"MagicDefinition`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicSimulate`*"];
ClearAll @@ Names["MagicSimulate`*"];

setfounderFGL::usage = "setfounderFGL  "

setIniPop::usage = "setIniPop  "

setIniPopFGL::usage = "setIniPopFGL  "

simIniPopFGLDO::usage = "simIniPopFGLDO  "

simPedigree::usage = "simPedigree  "

simPedigreeDO::usage = "simPedigreeDO  "

appendInbreeding::usage = "appendInbreeding  "

isCenterAlignment::usage = "An option of pedigreePlot to specifiy if plot is centered or not"

xShiftRight::usage = "An option of xShiftRight to specifiy the amounts of shift right in each generation"

pedigreeMemberList::usage = "pedigreeMemberList is an option"

obligateRescale::usage = "obligateRescale is an option"

pedigreePlot::usage = "pedigreePlot  "

simMagicFgl::usage = "simMagicFgl[mateScheme,founderFGL, chrLength,interferStrength, isObligate,isOogamy] simulates the breeding pedigree, and then drops founder genome lables (FGL)  on the pedigree to obtain FGL list of each individual at all generations. chrLength is a list of  chromosome length for each pair of chromosomes in centiMorgan. founderFGL speices the founder haploid genome labels, a matrix with dimension nfounder x 2. mateScheme is a list of mating scheme. isOogammy denotes if the population is Oogammy. isObligate denote where \nThe format of the return result. The first element of result is the input parameters, the result[[1]]={chrLength, founderFGL, mateScheme, isOogamy, isObligate, interferStrength}. The result[[i+2]] saves the F(i) population, a list of individuals in order of their IDs from 1 to the size of the population. Each individual is reprsented as {generation id,family id,individual id, individual gender, {father id,mother id},fgllist}, where family id is same as individual id in current version. The gender is represented as 1= female, 2=male, and 0=hermaphrodite. The individual's genome is represented as a FGL list. The FGL list is a list of pairs of maternally derived chromosomes and paternally derived chromosomes. The dimensions of the FGL list nnPairChr x 2 because of diploid populations. The FGL of each chromesomes is a list of breakpoints, for example, {{0,10},{47.0, 15},{67.7,2},{100,-1}} denotes the FGL of the segment [0,47.0) is 10, the FGL of the segment [47.0,67.7) is 15, the FGL of the segment [67.7, 100] is 2. The last element {Length of the chromosome (cM), -1} represents the end of the chromosome."

toFglDiplo::usage = "toFglDiplo  "

toFglDiploGrid::usage = "toFglDiploGrid  "

calFglDiplo::usage = "calFglDiplo  "

calFglDiploGrid::usage = "calFglDiploGrid  "

calIdentityJunc::usage = "calIdentityJunc  "

calFglSummary::usage = "calFglSummary  "

simOrigPrior::usage = "simOrigPrior  "

simIdentityPrior::usage = "simIdentityPrior  "

simAncestryPrior::usage = "simAncestryPrior  "

compilePopdesign::usage = "compilePopdesign  "

compilePopPedigree::usage = "compilePopPedigree  "

nextGeneration::usage = "nextGeneration  "

crossover::usage = "crossover  "

Begin["`Private`"]
(* Implementation of the package *)

compilePopdesign[popid_String] :=
    Module[ {nfounder,matescheme, pop, inbredscheme, nintercross, ninbred},
        pop = ToLowerCase[popid];
        Switch[pop,
         "f2", pop = "2wayril-self1",
         "cp" | "f1", pop = "2wayril-self0";
        ];
        Which[
         StringMatchQ[pop, DigitCharacter .. ~~ "wayril-self" ~~ DigitCharacter ..],
         {nfounder, ninbred} = ToExpression[StringSplit[pop, "wayril-self"]];
         inbredscheme = "Selfing";
         nintercross = Log[2, nfounder];
         matescheme = Join[Table["Pairing", {nintercross}],Table[inbredscheme, {ninbred}]],
         StringMatchQ[pop, DigitCharacter .. ~~ "wayril-sib" ~~ DigitCharacter ..],
         {nfounder, ninbred} = ToExpression[StringSplit[pop, "wayril-sib"]];
         inbredscheme = "Sibling";
         nintercross = Log[2, nfounder] - 1;
         matescheme = Join[Table["Pairing", {nintercross}],Table[inbredscheme, {ninbred}]];,
         True,
         Print["compilePopdesign: unknow population " popid,"!"];
         Abort[];
         ];
        {nfounder, matescheme}
    ]   
    
compilePopPedigree[popid_String, isoogamy_: False] :=
    Module[ {pop, k, s, temp, ped, nfounder, matescheme,j},
        pop = ToLowerCase[popid];
        If[ StringMatchQ[pop, "bc" ~~ DigitCharacter ..],
            pop = pop <> "-self0"
        ];
        If[ StringMatchQ[pop, "bc" ~~ DigitCharacter .. ~~ "-self" ~~ DigitCharacter ..],
            {k, s} = ToExpression[StringDelete[StringSplit[pop, "-"], "bc" | "self"]];
            ped = {{"Generation", "MemberID","Female=1/Male=2/Hermaphrodite=0", {"MotherID","FatherID"}}, 
                    {0, 1, 0, {0, 0}}, {0, 2, 0, {0, 0}}, {1, 3,0, {1, 2}}};
            temp = Join[Table[Range[ped[[-1, j]] + 1, k + ped[[-1, j]]], {j, 2}], {Table[ped[[-1, 3]], {k}], 
               Thread[{1, Range[ped[[-1, 2]], k + ped[[-1, 2]] - 1]}]}];
            ped = Join[ped, Transpose[temp]];
            temp = Join[Table[Range[ped[[-1, j]] + 1, s + ped[[-1, j]]], {j, 2}], {Table[ped[[-1, 3]], {s}], 
               Transpose[Table[Range[ped[[-1, 2]], s + ped[[-1, 2]] - 1], {2}]]}];
            ped = Join[ped, Transpose[temp]],
            {nfounder, matescheme} = compilePopdesign[pop];
            ped = simPedigree[nfounder, isoogamy, matescheme]
        ];
        ped
    ]
      
RandomDerangement[n_] :=
    Module[ {set = Range[n], res},
        res = RandomSample[set];
        While[Count[res - set, 0] != 0,
         res = RandomSample[set]
        ];
        res
    ]

matingParents[matingScheme_, n_Integer] :=
    Module[ {i, j,temp,nrep},
        Which[
            StringMatchQ[matingScheme, "WF-NE-" ~~ DigitCharacter ..],
            temp = ToExpression[Last[StringSplit[matingScheme, "-"]]];
            If[ temp>0,
                RandomChoice[Range[n], {temp, 2}],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],
            StringMatchQ[matingScheme, "RM1-NE-" ~~ DigitCharacter ..],
            temp = ToExpression[Last[StringSplit[matingScheme, "-"]]];
            If[ temp>=1,
                Table[RandomSample[Range[n], 2], {temp}],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],
            StringMatchQ[matingScheme, "RM2-NE-" ~~ DigitCharacter ..],
            temp = ToExpression[Last[StringSplit[matingScheme, "-"]]];
            If[ temp>=2,
                temp = Transpose[{RandomChoice[Range[1, n, 2], Round[temp/2]], RandomChoice[Range[2, n, 2], Round[temp/2]]}];
                (*temp = Table[RandomSample[Range[n], 2], {Round[temp/2]}];*)
                RandomSample[Join[temp, temp]],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],
            (*StringMatchQ[matingScheme, "WF-E-" ~~ DigitCharacter ..],
            temp = ToExpression[Last[StringSplit[matingScheme, "-"]]];
            If[ temp>0,
                RandomSample[Transpose[{Range[temp], RandomSample[Range[temp]]}]],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],            
            StringMatchQ[matingScheme, "RM1-E-" ~~ DigitCharacter ..],
            temp = ToExpression[Last[StringSplit[matingScheme, "-"]]];
            If[ temp>=1,
                RandomSample[Transpose[{Range[temp], RandomDerangement[temp]}]],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],*)
            StringMatchQ[matingScheme, "Selfing" ~~ DigitCharacter .. ~~ "-" ~~ DigitCharacter ..],
            temp = First[StringCases[matingScheme,"Selfing" ~~ temp2 : DigitCharacter ... ~~ "-" ~~ temp : DigitCharacter .. :> {temp2, temp}]];
            {nrep, temp} = ToExpression[#] & /@ temp;
            If[ !(IntegerQ[nrep] && nrep >= 1),
                Print["matingParents: the  replication number ", nrep," must be a postive integer in ",matingScheme,"!"];
                Abort[];
            ];
            If[ !(IntegerQ[temp/nrep] && temp/nrep>= 1),
                Print["matingParents: the population size ", temp,
                    " must be a postive integer times the replication number ", nrep, 
                    " in the mating scheme " matingScheme,"!"];
                Abort[];
            ];
            If[ temp/nrep>n,
                Print["matingParents: the population size ", temp,
                    " must be no greater than  the replication number ", nrep, 
                    " times the previous population size ", n, " in the mating scheme ", matingScheme,"!"];
                Abort[];
            ];
            temp = Sort[RandomSample[Range[n], temp/nrep]];
            temp = Flatten[Transpose[Table[temp, {nrep}]]];
            Transpose[{temp, temp}],
            True,            
            (*"WF-NE": fixed populatin size, allow selfing, random number of offspring *)
            (*"WF-E": fixed populatin size, allow selfing, each parent contributes eactly 2 gametes*)
            (*"RM11-NE": fixed populatin size, not allow selfing, random number of offspring *)        
            (*"RM11-E": fixed populatin size, not allow selfing,  each parent contributes eactly 2 gametes*)
            Switch[matingScheme, 
             "Selfing", Table[{i, i}, {i, n}],
             "Sibling", Flatten[Table[{2 i - 1, 2 i}, {i, Floor[n/2]}, {2}], 1],         
             "Pairing", Partition[Range[n], 2],
             "CrossPairing", Partition[Flatten[Partition[Range[n], 4][[All, {1, 3, 2, 4}]]], 2],
             "HalfDiallel", RandomSample[Flatten[Table[{i, j}, {i, n - 1}, {j, i + 1, n}], 1]],
             "FullDiallel", RandomSample[DeleteCases[Flatten[Array[List, {n, n}], 1], _?(#[[1]] == #[[2]] &)]],    
             "CM-E", Partition[Append[Range[n], 1], 2, 1],
             "CPM-E",RotateRight[Riffle[Partition[Range[n], 2], Partition[Range[n], 2]]],
             "MAI-E",Join[Partition[Range[n], 2], Partition[Range[n], 2]],
             "WF-NE", RandomChoice[Range[n], {n, 2}],
             "WF-E",RandomSample[Transpose[{Range[n], RandomSample[Range[n]]}]],
             "RM1-NE",Table[RandomSample[Range[n], 2], {n}],         
             "RM1-E",RandomSample[Transpose[{Range[n], RandomDerangement[n]}]],         
             (*"RM11-NE",Transpose[{RandomChoice[Range[1, n, 2], n],RandomChoice[Range[2, n, 2], n]}],        
             "RM11-E",
             temp = Transpose[{RandomSample[Range[1, n, 2]], RandomSample[Range[2, n, 2]]}];
             temp2 = Transpose[{RandomSample[Range[1, n, 2]],RandomSample[Range[2, n, 2]]}];
             RandomSample[Join[temp, temp2]],*)
             "RM2-NE",    
             temp = Transpose[{RandomChoice[Range[1, n, 2], n/2], RandomChoice[Range[2, n, 2], n/2]}];
             (*temp = Table[RandomSample[Range[n], 2], {Round[n/2]}];*)
             RandomSample[Join[temp, temp]],
             "RM2-E",
             temp = Transpose[{RandomSample[Range[1, n, 2]], RandomSample[Range[2, n, 2]]}];
             RandomSample[Join[temp, temp]],
              _, Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                 Abort[];
             ]
        ]
    ]

matingParents[matingScheme_, gender_List] :=
    Module[ {gendergroup, nm,nf,temp,i,rep,pairs,newgender},
        (*1= female, 2=male*)
        gendergroup = Table[Flatten[Position[gender, i]], {i, 1, 2}];
        (*assuming nf=nm in the current version*)
        {nf,nm} = Length[#]&/@gendergroup;
        (* ordered mating pair: (mother=1, father=2)*)
        Which[
            StringMatchQ[matingScheme, "RM1-NE"~~ "" | "-" ~~DigitCharacter ...],
            temp = StringSplit[matingScheme, "-"];
            If[ Length[temp]>2,
                nf = nm = ToExpression[Last[temp]]/2
            ];
            If[ IntegerQ[nf],
                pairs = Transpose[{RandomChoice[gendergroup[[1]], nf+nm], RandomChoice[gendergroup[[2]],nf+nm]}],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],            
            StringMatchQ[matingScheme, "RM2-NE"~~ "" | "-" ~~DigitCharacter ...],
            temp = StringSplit[matingScheme, "-"];
            If[ Length[temp]>2,
                nf = nm = ToExpression[Last[temp]]/2
            ];
            If[ IntegerQ[nf],
                temp = Transpose[{RandomChoice[gendergroup[[1]], nf], RandomChoice[gendergroup[[2]],nm]}];
                pairs = RandomSample[Join[temp, temp]],
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],
            (*StringMatchQ[matingScheme, "RM11-NE"~~ "" | "-" ~~DigitCharacter ...],
            temp = StringSplit[matingScheme, "-"];
            If[Length[temp]>2,nf=nm=ToExpression[Last[temp]]/2];        
            If[ IntegerQ[nf],
                temp = Transpose[{RandomChoice[gendergroup[[1]], nf], RandomChoice[gendergroup[[2]],nm]}];
                (*Riffle and alternative male and female ==> each sample of two parents produces one son, and one daughter*)
                 pairs = Riffle[temp, temp],                
                Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];
            ],*)            
            StringMatchQ[matingScheme, "RM1-E"~~ "" | "-" ~~DigitCharacter ...],
            temp = StringSplit[matingScheme, "-"];
            If[ Length[temp]>2,
                temp = ToExpression[Last[temp]]/2;
                If[ (temp!=nf)||(temp!=nm),
                    Print["Warning: the population size is fixed to "<>ToString[Length[Flatten[gendergroup]]]<>", and not equal to the size specified in ", matingScheme, "!"]
                ];
            ];
            temp = Table[Transpose[{RandomSample[gendergroup[[1]]], RandomSample[gendergroup[[2]]]}],{2}];
            pairs = RandomSample[Join@@temp],  
            StringMatchQ[matingScheme, "RM2-E"~~ "" | "-" ~~DigitCharacter ...],
            temp = StringSplit[matingScheme, "-"];
            If[ Length[temp]>2,
                temp = ToExpression[Last[temp]]/2;
                If[ (temp!=nf)||(temp!=nm),
                    Print["Warning: the population size is fixed to "<>ToString[Length[Flatten[gendergroup]]]<>", and not equal to the size specified in ", matingScheme,"!"]
                ];
            ];
            temp = Transpose[{RandomSample[gendergroup[[1]]], RandomSample[gendergroup[[2]]]}];
             (*Genders for the two offspring produced in each sampling are randomly assigned*)
            pairs = RandomSample[Join[temp, temp]],            
             (*"RM11-E",
             temp = Transpose[{RandomSample[gendergroup[[1]]], RandomSample[gendergroup[[2]]]}];
             (*Riffle and alternative male and female ==> each sample of two parents produces one son, and one daughter*)
             pairs = Riffle[temp, temp],*)   
            StringMatchQ[matingScheme, "Sibling" ~~ DigitCharacter ..],
            rep = ToExpression[StringDrop[matingScheme, 7]];
            temp = Transpose[gendergroup[[All, ;; Min[nm, nf]]]];
            pairs = Flatten[Transpose[Table[temp, {rep}]], 1],
            True,  
            Switch[matingScheme,
             "Sibling",
             temp = Transpose[gendergroup[[All, ;; Min[nm,nf]]]];
             pairs = Riffle[temp,temp],         
             "Pairing", 
             pairs = Transpose[gendergroup[[All, ;; Min[nm,nf]]]],
             "CrossPairing",
             temp = Flatten[Reverse[#] & /@ Partition[gendergroup[[2]], 2]];
             temp = Join[temp, Complement[gendergroup[[2]], temp]];
             pairs = Transpose[{gendergroup[[1]], temp}[[All, ;; Min[nm, nf]]]],             
             "HalfDiallel", 
             pairs = RandomSample[Flatten[Outer @@ Prepend[gendergroup, List], 1]],         
             "FullDiallel", 
             temp = Flatten[Outer @@ Prepend[gendergroup, List], 1];
             pairs = RandomSample[Riffle[temp, temp]],
             "CM-E", pairs = Partition[Append[Range[nm+nf], 1], 2, 1],
             "CPM-E",pairs = RotateRight[Riffle[Partition[Range[nm+nf], 2], Partition[Range[nm+nf], 2]]],
             "MAI-E",pairs = Join[Partition[Range[nm+nf], 2], Partition[Range[nm+nf], 2]],
             _, Print["matingParents: Wrong mating Scheme ",matingScheme,"!"];
                Abort[];             
             ];
        ];
        newgender = Table[1,{Length[pairs]}];
        newgender[[2;;;;2]] = 2;
        (*Print["{matingScheme,pairs,newgender}=",{matingScheme,pairs,newgender}];*)
        {pairs,newgender}
    ]
    
(*any pedigree member = {generation, individual id, gender, {mother id,father id}}*)
(*gender:1=female, 2=male,0=hermaphrodite*)
nextGeneration[pop_, matingScheme_, size_: 1] :=
    Module[ {genders,newgenders,pairs, next},
        genders = pop[[All,3]];
        Which[ Union[genders] == {0},
            pairs = matingParents[matingScheme, Length[pop]];
            next = Thread[{-1, -1, 0, pairs}],
            Union[genders] == {1,2},
            {pairs,newgenders} = matingParents[matingScheme, genders];
            next = Thread[{-1, -1, newgenders, pairs}],
            True,
            Print[genders];
            Print["nextGeneration: there are no male (or female) indiviudals !"];
            Abort[];             
        ];
        next[[All, 2]] = If[ Length[next]==Length[pop],
                             pop[[All, 2]],
                             Range[Length[next]]
                         ];
        next = Flatten[Transpose[Table[next, {size}]], 1];
        next
    ]

crossover[chrpair_,isObligate_,interferStrength_] :=
    Module[ {len, breaks, temp,complex, rec, map, ls, res, ind, ii, jj,ch},
        len = chrpair[[1, -1, 1]];
        (*len is in unit of centiMorgan*)
        (*breaks = RandomReal[{0, len}, Max[1,RandomInteger[PoissonDistribution[len/100]]]];*)
        temp = RandomInteger[PoissonDistribution[len/100]] (1 + interferStrength);
        If[ isObligate,
            temp = Max[1 + interferStrength, temp]
        ];
        breaks = Sort[RandomReal[{0, len}, temp]];
        If[ Length[breaks]>=3,
            breaks = breaks[[RandomInteger[{1, 1 + interferStrength}] ;; ;;1 + interferStrength]];
        ];
        breaks = Union[{0}, breaks, {len}];
        complex = Table[
          rec = Transpose[{breaks, ch[[intervalIndexFunction[ch[[All, 1]]][breaks], 2]]}];
          map = intervalIndexFunction[rec[[All, 1]]];
          ls = Join[rec[[;; -2]], ch[[2 ;; -2]]];
          Gather[Transpose[{ls, map[ls[[All, 1]]]}], Last[#1] == Last[#2] &][[All, All, 1]], {ch, chrpair}];
        res = ind = Range[Length[breaks] - 1];
        ii = ind[[RandomChoice[{1, 2}] ;; ;; 2]];
        jj = Complement[ind, ii];
        res[[ii]] = complex[[1, ii]];
        res[[jj]] = complex[[2, jj]];
        Append[SplitBy[Flatten[res, 1], Last][[All, 1]], {len, -1}]
    ]

setfounderFGL[nFounder_, isInbredFounder_] :=
    Module[ {i},
        If[ isInbredFounder,
            Table[{i, i}, {i, nFounder}],
            Table[{2 i - 1, 2 i}, {i, nFounder}]
        ]
    ]
      
(*any pedigree member = {generation, individual id, gender, {mother id, father id}};
generation: non-overlapping, starting from 0 for the founder populations; 
            set generation =-1 for overlapping generations.
gender:1=female,2=male, 0=hermaphrodite;
individual id: natural number starting from 1;
The parents of founders are set to {0,0};founders are always in the beginning. 
All the parents precede their offspring.
*)
(*If isOogammy=True, the founder population consists of 1=female, 2=male, female, male, and so on;
If is Oogammy =False, the gender =0 for all founders*)

setIniPop[nFounder_,isOogamy_] :=
    Module[ {iniPop,i},
        iniPop = Table[{0, -1, 0, {0, 0}},{nFounder}];
        iniPop[[All, 2]] = Range[nFounder];
        If[ isOogamy,
            Do[iniPop[[i ;; ;; 2, 3]] = i, {i, 1, 2}]
        ];
        (*Add column heads to the list of {generation, member id, gender, mother id, father id}*)
        Join[{{"Generation", "MemberID","Female=1/Male=2/Hermaphrodite=0",{"MotherID","FatherID"}}},iniPop]
    ]  
    
setIniPopFGL[iniPop_,founderFGL_, chrLength_] :=
    Module[ {inipopFGL,nFounder,ls,ls2,i},
        nFounder = Length[iniPop]-1;
        inipopFGL = Join[Rest[iniPop],Table[{0},{nFounder}],2];
        inipopFGL[[All, -1]] = Table[
          ls = {{0, founderFGL[[i, 1]]}, {#, -1}} & /@ chrLength;
          ls2 = {{0, founderFGL[[i, 2]]}, {#, -1}} & /@ chrLength;
          Transpose[{ls, ls2}], {i, nFounder}];
        Join[{Join[First[iniPop],{"Genome"}]},inipopFGL]
    ]    

simIniPopFGLDO[nPower_, preCCfreq_, popSize_, chrLength_,interferStrength_, isObligate_, isOogamy_] :=
    Module[ {nFounder, founderFGL, inipop, inipopFGL, mateScheme,i},
        nFounder = 2^nPower;
        inipop = Table[
          founderFGL = RandomSample[Range[nFounder]];
          founderFGL = Transpose[{founderFGL, founderFGL}];
          inipop = setIniPop[nFounder, isOogamy];
          inipopFGL = setIniPopFGL[inipop, founderFGL, chrLength];
          mateScheme = Join[Table["Pairing", {nPower - 1}],Table["Sibling", {RandomChoice[Rule @@ Reverse[Transpose[preCCfreq]]]}]];
          simMagicFgl[mateScheme, inipopFGL, interferStrength, isObligate, isOogamy][[If[ OddQ[i],
                                                                                          -2,
                                                                                          -1
                                                                                      ]]], {i, popSize}];
        inipop[[All, 1]] = 0;
        inipop[[All, 2]] = Range[popSize];
        inipop[[All, 4]] = {0, 0};
        Join[{{"Generation", "MemberID","Female=1/Male=2/hermaphrodite=0",{"MotherID","FatherID"}, "Genome"}}, inipop]
    ]
 
    
simPedigree[iniPop_,mateScheme_] :=
    Module[ {ini,popped,i},
        ini = SplitBy[Rest[iniPop], First];
        popped = Join[Most[ini],FoldList[nextGeneration[#1, #2] &, Last[ini], mateScheme]];
        Do[
          popped[[i, All, 4]] = popped[[i - 1, All, 2]][[#]] & /@ popped[[i, All, 4]];
          popped[[i, All, 2]] +=Max[popped[[i - 1, All, 2]]] - popped[[i, 1, 2]] + 1, {i, 
           Length[ini] + 1, Length[popped]}];
        popped[[All, All, 1]] = Range[0, Length[popped] - 1];
        popped = Flatten[popped, 1];
        Join[{First[iniPop]},popped]
    ]
    
simPedigree[nFounder_, isOogamy_,mateScheme_] :=
    simPedigree[setIniPop[nFounder,isOogamy],mateScheme]    

appendInbreeding[pedigree_, inbredScheme_, ginbred_] :=
    Module[ {isOogamy, pop0, newt0, newind0, cycles, res, genderset,i},
        isOogamy = pedigree[[-1, 3]] != 0;
        pop0 = Last[SplitBy[pedigree[[2 ;;]], First]];
        If[ ! (Length[pop0] == 2),
            Print["appendInbreeding: the last generation of the input pedigre must have exactly two individuals!"];
            Abort[]
        ];
        If[ ! MemberQ[{"AlternatingBackcross", "FatherDaughter"}, inbredScheme],
            Print["appendInbreeding: inbredScheme must be AlternatingBackcross (=alternating backcross swaps between mother-son and father-daugher) or FatherDaughter (alternating backcross swaps between father-daugher and sibling)!"];
            Abort[]
        ];
        newt0 = pop0[[-1, 1]];
        newind0 = 
         If[ NumberQ[pop0[[-1, 2]]],
             pop0[[-1, 2]],
             Length[pedigree] - 1
         ];
        cycles = Ceiling[ginbred/2];
        Switch[inbredScheme,
         "AlternatingBackcross",
         (*Alternating backcross swaps between mother-son and father-daugher*)
         res = Table[0, {cycles 2}];
         genderset = If[ isOogamy,
                         {2, 1},
                         {0, 0}
                     ];
         res[[1]] = {newt0 + 1, newind0 + 1, 0, pop0[[All, 2]]};
         res[[2]] = {newt0 + 2, newind0 + 2, 0, {pop0[[1, 2]], newind0 + 1}};
         res[[;; 2, 3]] = genderset;
         Do[
          res[[2 i - 1]] = {newt0 + 2 i - 1, newind0 + 2 i - 1, 0, res[[{2 (i - 1), 2 (i - 1) - 1}, 2]]};
          res[[2 i]] = {newt0 + 2 i, newind0 + 2 i, 0, res[[{2 (i - 1), 2 i - 1}, 2]]};
          res[[{2 i - 1, 2 i}, 3]] = genderset, {i, 2, cycles}],
         "FatherDaughter",
         (*Alternating backcross swaps between father-daugher and sibling*)
         res = Table[0, {cycles 3}];
         genderset = If[ isOogamy,
                         {1, 1, 2},
                         {0, 0, 0}
                     ];
         res[[1]] = {newt0 + 1, newind0 + 1, 0, pop0[[All, 2]]};
         res[[2]] = {newt0 + 2, newind0 + 2, 0, {newind0 + 1, pop0[[2, 2]]}};
         res[[3]] = {newt0 + 2, newind0 + 3, 0, {newind0 + 1, pop0[[2, 2]]}};
         res[[;; 3, 3]] = genderset;
         Do[
          res[[3 i - 2]] = {newt0 + 2 i - 1, newind0 + 3 i - 2, 0, res[[{3 (i - 1) - 1, 3 (i - 1)}, 2]]};
          res[[3 i - 1]] = {newt0 + 2 i, newind0 + 3 i - 1, 0, res[[{3 i - 2, 3 (i - 1)}, 2]]};
          res[[3 i]] = {newt0 + 2 i, newind0 + 3 i, 0, res[[{3 i - 2, 3 (i - 1)}, 2]]};
          res[[3 i - 2 ;; 3 i, 3]] = genderset, {i, 2, cycles}],
         _, 0
         ];
        Join[pedigree, Flatten[Take[SplitBy[res, First], ginbred], 1]]
    ]

simPedigreeDO[nPower_, preCCfreq_, crossPopSize_, gCross_, crossScheme_] :=
    Module[ {isOogamy = True, nFounder = 2^nPower, inipop, inipop2, pedls,
       nsibling, mateScheme, i0, i,ped, founders, crossped, nf, temp,res},
        inipop = setIniPop[nFounder, isOogamy];
        pedls = Table[
          inipop2 = inipop;
          temp = RandomSample[Range[nFounder]];
          temp[[2;;;;2]]+=nFounder;
          inipop2[[2 ;;, 2]] = temp;
          nsibling = RandomChoice[Rule @@ Reverse[Transpose[preCCfreq]]];
          mateScheme = Join[{"Pairing", "Pairing"}, Table["Sibling", {nsibling}]];
          ped = simPedigree[inipop2, mateScheme];
          ped, {crossPopSize}];
        i0 = nFounder + 1;
        Do[
         temp = Range[Length[pedls[[i]]] - i0] + 
           If[ i == 1,
               2 nFounder,
               Max[pedls[[i - 1, i0 + 1 ;;, 2]]]
           ];
         pedls[[i, i0 + 1 ;;, {2, 4}]] = pedls[[i, i0 + 1 ;;, {2, 4}]] /. Thread[pedls[[i, i0 + 1 ;;, 2]] -> temp], {i, Length[pedls]}];
        (*Print[pedigreePlot[#,ImageSize\[Rule]300]&/@pedls];*)
        pedls = MapIndexed[Drop[#1, {-(1 + Boole[EvenQ[First[#2]]])}] &, pedls];
        (*founders are female, male, female,.. alternatively*)
        founders = pedls[[All, -1, 2]];
        inipop = setIniPop[Length[founders], isOogamy];
        crossped = simPedigree[inipop, Table[crossScheme, {gCross}]];
        nf = Length[founders] - Boole[OddQ[Length[founders]]];
        temp = Max[pedls[[All, 2 ;;, 2]]];
        temp = Thread[crossped[[2 ;;, 2]] -> Join[Take[founders, nf],Range[temp + 1, temp + Length[crossped] - 1 - nf]]];
        crossped[[2 ;;, {2, 4}]] = crossped[[2 ;;, {2, 4}]] /. temp;
        crossped[[2 ;; nf + 1]] = pedls[[;; nf, -1]];
        crossped[[nf + 2 ;;, 1]] += Max[crossped[[2 ;; nf + 1, 1]]];
        pedls = Most[#] & /@ pedls;
        res = Join[pedls[[1]], Sequence @@ pedls[[2 ;;, 2 ;;]], Rest[crossped]];
        res[[2 ;;, 1]] -= (res[[-1, 1]] - gCross);
        res[[2 ;;]] = SortBy[res[[2 ;;]], 1];
        res = DeleteDuplicates[res];
        res
    ]
         
Options[pedigreePlot] = Join[{isCenterAlignment->True, xShiftRight-> None},Options[LayeredGraphPlot]];
pedigreePlot[pedigree_, opts : OptionsPattern[]] :=
    Module[ {ped = pedigree,data,yls,xyls,alignls,i,y,bool,shift,xy},
        (*assume that the pedigree members are sorted by generation in the firsted column, so that parents are always before children*)
        data = DeleteCases[ped[[2 ;;, {2, 4}]], {_, {0, 0}}];
        data = Flatten[Thread[#[[2]] -> #[[1]]] & /@ data];
        yls = SplitBy[ped[[2 ;;, {1, 3}]], First];
        alignls = TrueQ[#]&/@Flatten[{OptionValue[isCenterAlignment]}];
        alignls = PadRight[alignls, Length[yls], Last[alignls]];
        bool = Boole[And @@ alignls];
        shift = Flatten[{OptionValue[xShiftRight]}];
        If[ First[shift]===None,
            shift = Table[0,{Length[yls]}],
            shift = PadRight[shift,Length[yls],0]
        ];
        xyls = Table[
           y = yls[[i]];
           xy = If[ alignls[[i]],
                    Transpose[{Range[-Ceiling[Length[y]/2] + 1, Floor[Length[y]/2]]- 
                       If[ EvenQ[Length[y]],
                           1/2,
                           0
                       ] bool, Max[yls[[All, 1]]] - y[[All, 1]]}],
                    Transpose[{Range[1, Length[y]] + y[[1, 2]] - 2, Max[yls[[All, 1]]] - y[[All, 1]]}]
                ];
           xy[[All,1]]+=shift[[i]];
           xy, {i, Length[yls]}];
        xyls = Flatten[xyls, 1];
        xyls = Thread[ped[[2 ;;, 2]] -> xyls];
        LayeredGraphPlot[data, Join[Evaluate[FilterRules[{opts}, Options[LayeredGraphPlot]]],
            {VertexLabeling->False,DirectedEdges -> True,VertexCoordinateRules -> xyls,
             VertexRenderingFunction -> ({LightYellow,Opacity[0.5], EdgeForm[Red],
             If[ ped[[Position[ped[[2 ;;, 2]], #2][[1, 1]] + 1, 3]]==2,
                 Rectangle[# - .2, # + .2],
                 Disk[#, .2]
             ], 
             Black,Text[Style[#2, 12], #1]} &)}]]
    ]    
  
Options[simMagicFgl] = Options[calFglSummary] = {pedigreeMemberList->Automatic,obligateRescale->True}
(*If isOogamy=True,the last pair of chromosomes are sex chromosomes.*)(*XY for male and XX for female.*)
(*for each linkage group = {maternally derived chromosome, paternally derived chromosome}*)
(*If ith individual is female (gender =1), the paternal dervied chromosome is the X chromosome*)
(*popfgl[[1]] ={popPed,iniPopFGL, interferStrength,isObligate, isOogamy}*)
simMagicFgl[inputpopPed_,iniPopFGL_, interferStrength_,isObligate_,isOogamy_,opts : OptionsPattern[]] :=
    Module[ {popPed=inputpopPed,rule,popfgl, nFounder,gender, chrLength,memberlist,pos,factor,i},    	
        nFounder = Length[iniPopFGL]-1;
        If[popPed[[2 ;; nFounder + 1,2]]=!= iniPopFGL[[2 ;;,2]],
        	rule = Thread[popPed[[2 ;; nFounder + 1,2]]->iniPopFGL[[2 ;;nFounder + 1,2]]];    
        	popPed[[2 ;;, {2,4}]]= popPed[[2 ;;, {2,4}]]/.rule;
        ];        
        If[ popPed[[2 ;; nFounder + 1,;;4]] =!= iniPopFGL[[2 ;;, ;;4]],        	
        	(*Put[popPed,iniPopFGL,"temptest.txt"];
            Print["here1"];*)               
            Print["Inconsistent founders between popPed and iniPopFGL!"];
            Abort[]
        ];
        (*By defaults, startpoint of each chromosome is 0, and the endpoint is the length of the chromosome*)
        chrLength = iniPopFGL[[-1, -1, All, 1, -1, 1]];
        popfgl = Join[Rest[popPed], Table[{0}, {Length[popPed]-1}], 2];
        popfgl[[;; nFounder, -1]] = iniPopFGL[[2;;,-1]];
        Do[
         popfgl[[i, -1]] = Transpose[Map[crossover[#, isObligate, interferStrength] &, popfgl[[popfgl[[i, 4]], 5]], {2}]];
         If[ isOogamy&&Length[chrLength]>1,
             gender = popfgl[[i, 3]];
             (*popfgl[[i,4,2]] is the father id*)
             popfgl[[i, -1, -1, 2]] = popfgl[[popfgl[[i, 4, 2]], -1, -1, gender]];
         ], {i, nFounder + 1, Length[popfgl]}];
        memberlist = OptionValue[pedigreeMemberList];
        Which[
            memberlist === All,
            0,
            memberlist === Automatic,
            popfgl = Select[popfgl, First[#] == popfgl[[-1, 1]] &],            
            VectorQ[memberlist],
            pos = Flatten[Position[popfgl[[All, 2]], #, {1}, Heads -> True] & /@ memberlist];
            If[ pos==={},
                Print["pedigreeMemberList contains no individuals in the input peidgree!"];
                Abort[]
            ];
            popfgl = popfgl[[pos]];            
        ];     
        (*correcting for obligate crossover*)
        If[ isObligate&& TrueQ[OptionValue[obligateRescale]],
            (*transform from centiMorgan into Morgan*)
            chrLength/=100;
            factor = 1. + Exp[-chrLength]/(chrLength);
            Do[popfgl[[All, -1, i, All, -1, 1]]*= factor[[i]], {i, Length[factor]}];
        ];
        Join[{{popPed,iniPopFGL,interferStrength,isObligate, isOogamy,memberlist}},popfgl]
    ]

(*any pedigree member = {generation, individual id, gender, mother id, father id}*)
(*popfgl[[1]] ={popPed,iniPopFGL, interferStrength,isObligate, isOogamy}*)
simMagicFgl[mateScheme_?(VectorQ[#, StringQ]&),iniPopFGL_,interferStrength_,isObligate_,isOogamy_,opts : OptionsPattern[]] :=
    Module[ {popped,iniPop,popfgl},
        iniPop = iniPopFGL[[All,;;-2]];
        popped = simPedigree[iniPop,mateScheme];        
        popfgl = simMagicFgl[popped, iniPopFGL, interferStrength,isObligate, isOogamy,opts];
        popfgl[[1]] = Join[popfgl[[1]],{mateScheme}];
        popfgl
    ] 
                      
toFglDiplo[fglhap_] :=
    Module[ {fun},
        fun[z_] :=
            Module[ {junc, temp,l},
                junc = Union[Flatten[z[[All, All, 1]]]];
                temp = Transpose[Table[z[[l, intervalIndexFunction[z[[l, All, 1]]][junc], 2]], {l,Length[z]}]];
                Transpose[{junc, temp}]
            ];
        Map[fun, fglhap, {Depth[fglhap] - 4}]
    ]
(*popfgl[[1]] ={popPed, iniPopFGL,interferStrength,isObligate,isOogamy}*)
(*or popfgl[[1]] ={popPed, iniPopFGL,interferStrength,isObligate, isOogamy, mateScheme}*)
calFglDiplo[popfgl_] :=
    Module[ {popfgldiplo = popfgl,chrLength, isObligate},
        isObligate = popfgldiplo[[1, 4]];
        chrLength = popfgldiplo[[-1, -1, All, 1, -1, 1]];
        popfgldiplo[[2 ;;, -1]] = toFglDiplo[popfgldiplo[[2 ;;, -1]]];
        popfgldiplo
    ] 
    

toFglDiploGrid[fgldiplo_, grid_] :=
    Module[ {fun},
        fun[z_] :=
            Module[ {temp,ch},
                Table[
                 temp = 
                  z[[ch, intervalIndexFunction[z[[ch, All, 1]]][grid[[ch]]], 2]];
                 Transpose[{grid[[ch]], temp}], {ch, Length[grid]}]
            ];
        Map[fun, fgldiplo, {Depth[fgldiplo] - 5}]
    ]  
    
calFglDiploGrid[popfgldiplo_, grid_] :=
    Module[ {popgrid = popfgldiplo,chrLength},
        chrLength = popgrid[[-1, -1, All, -1, 1]];
        If[ Length[grid] != Length[chrLength],
            Print["The #chromosomes are different between grid and popfgldiplo!"];
            Abort[];
        ];
        If[ ! ((And @@ Thread[grid[[All, -1]] < chrLength]) && (And @@Thread[grid[[All, 1]] >= 0])),
            Print["The locations of some grid points are out of range!"];
            Abort[];
        ];
        popgrid[[2 ;;, -1]] = toFglDiploGrid[popgrid[[2 ;;, -1]], grid];
        popgrid
    ]
        
toIdentityJunction[diplo_] :=
    Module[ {label, rule, freq, f},
        label = {"J1112", "J1121", "J1122", "J1211", "J1213", "J1222", "J1232"};
        rule = {{{x_, x_}, {y_, y_}} :> "J1122", 
            {{x_, x_}, {x_, y_}} :>"J1112", {{x_, x_}, {y_, x_}} :>"J1121", 
            {{x_, y_}, {x_, x_}} :>"J1211", {{x_, y_}, {y_, y_}} :>"J1222", 
            {{x_, y_}, {x_, z_}} :> "J1213", {{x_, y_}, {z_, y_}} :>"J1232", 
            {{_, _}, {_, _}} :> "Wrong!"};
        freq = Replace[diplo, rule, {Depth[diplo] - 3}];
        f = Function[{x}, Count[x, #] & /@ label];
        Map[f, freq, {Depth[freq] - 2}]
    ]       

calIdentityJunc[popfgldiplo_] :=
    Module[ {pop = popfgldiplo,res, isOogamy, chrLength, geno, ibd, diplo, junc, ls,phi11},
        isOogamy = pop[[1, 5]];
        chrLength = pop[[-1, -1, All, -1, 1]];
        chrLength = chrLength/100.;
        geno = pop[[2 ;;, -1]];
        geno = geno[[All, All, All, 2]];
        ibd = geno[[All, All, 1]];
        ibd = Boole[Map[SameQ @@ # &, ibd, {2}]];
        diplo = Map[Partition[Most[#], 2, 1] &, geno, {2}];
        (*{"J1112", "J1121", "J1122", "J1211", "J1213", "J1222", "J1232"}*)
        junc = toIdentityJunction[diplo];
        junc = Map[Join[{Total[#[[{2, 3, 6, 7}]]], Total[#[[{1, 3, 4, 5}]]], Total[#]}, #] &, junc, {2}];
        (*{"Rm","Rp","rho","J1112", "J1121", "J1122", "J1211", "J1213", "J1222", "J1232"}*)
        If[ isOogamy && Length[chrLength] > 1,
            phi11 = Mean[Transpose[ibd[[All, ;; -2]]]];
            ls = Join[Transpose[{phi11}], 
              Total[junc[[All, ;; -2]], {2}]/Total[Most[chrLength]], 2];
            ls = Transpose[{ls, Join[Transpose[{ibd[[All, -1]]}], junc[[All, -1]]/Last[chrLength],2]}],
            phi11 = Mean[Transpose[ibd]];
            ls = Partition[Join[Transpose[{phi11}], Total[junc, {2}]/Total[chrLength], 2], 1];
        ];
        res = pop;
        res[[2 ;;, -1]] = ls;
        res
    ]    
    
(*popfgl[[1]] ={popPed, iniPopFGL,interferStrength,isObligate, isOogamy}*)
(*or popfgl[[1]] ={popPed,iniPopFGL,interferStrength,isObligate, isOogamy,mateScheme}*)
(*gender:1=female,2=male, 0=hermaphrodite;*)
calFglSummary[popfgljunc_] :=
    Module[ {pop = popfgljunc,resAA, tls,isOogamy,res, ls, temp, pos,gender},
        resAA = (Total[#]/Length[#]) & /@ SplitBy[Transpose[{pop[[2 ;;, 1]],pop[[2 ;;, -1, 1]]}], First];
        tls = resAA[[All, 1]];
        isOogamy = pop[[1, 5]];
        If[ isOogamy && Length[pop[[-1, -1]]]>1,
            res = ConstantArray[0, {Length[tls], 4}];
            res[[All, ;; 2]] = resAA;
            ls = Transpose[{pop[[2 ;;, {1, 3}]], pop[[2 ;;, -1, -1]]}];
            Do[
             temp = Total[#]/Length[#] & /@ SplitBy[Select[ls, (#[[1, 2]] == gender &)], First];
             temp[[All, 1]] = temp[[All, 1, 1]];
             pos = Flatten[Position[tls, #, {1}, 1, Heads -> False] & /@temp[[All, 1]]];
             res[[pos, gender + 2]] = temp[[All, 2]], {gender, 2}],
            res = resAA
        ];
        res[[All, 2 ;;, 1]] = N[res[[All, 2 ;;, 1]]];
        (*{t,{{finb,j1122,j1222,j1232}}*)
        (*{t,{finb,j1122,j1222,j1232},female={finb_mp,j1122mp,j1211mp,j1213mp, j1222mp,j1232mp},male={finb_mp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}*)
        Join[{pop[[1]]}, res]
    ]   

(*pop =popped or matescheme*)
calFglSummary[pop_, inipopFGL_,interferStrength_,isObligate_,isOogamy_,sampleSize_, opts : OptionsPattern[]] :=
    Module[ {popfgl,popfgldiplo,popfgljunc,summary,res,i},
        Monitor[
            res = Sum[
                popfgl = simMagicFgl[pop, inipopFGL, interferStrength, isObligate, isOogamy,opts];
                popfgldiplo = calFglDiplo[popfgl];
                popfgljunc = calIdentityJunc[popfgldiplo];
                summary = calFglSummary[popfgljunc];
                Rest[summary], {i, sampleSize}], ProgressIndicator[i, {1, sampleSize}]];
        Join[{Append[First[summary],sampleSize]}, res/sampleSize]
    ]
           
(*pop =popped or matescheme*)       
simOrigPrior[pop_, inipopFGL_,interferStrength_,isObligate_,isOogamy_,sampleSize_] :=
    Module[ {summary,i},
        summary = calFglSummary[pop, inipopFGL, interferStrength, isObligate, isOogamy, sampleSize, pedigreeMemberList -> Automatic,obligateRescale->True];
        summary = Rest[summary[[-1]]];
        (*{{"finb","J1112","J1121","J1122","J1211","J1213","J1222","J1232"},..} to 
          {{finbred,j1122,j1211,j1213,j1222,j1232},{finbred_mp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}}*)
        Table[{summary[[i, 1]], summary[[i, 4]],Mean[summary[[i, {2, 5}]]], summary[[i, 6]], 
                                                 Mean[summary[[i, {3, 7}]]], summary[[i, 8]]},{i,Min[2,Length[summary]]}]
    ]
    
(*pop =popped or matescheme*)       
simOrigPrior[pop_, founderFGL_,chrLength_,interferStrength_,isObligate_,isOogamy_,sampleSize_] :=
    Module[ {inipop,inipopFGL},
        inipop = setIniPop[Length[founderFGL],isOogamy];
        inipopFGL = setIniPopFGL[inipop,founderFGL, chrLength];
        simOrigPrior[pop, inipopFGL,interferStrength,isObligate,isOogamy,sampleSize]
    ]    

simIdentityPrior[pop_, inipopFGL_,interferStrength_,isObligate_,isOogamy_,sampleSize_, opts : OptionsPattern[]] :=
    Module[ {popfgl,popfgldiplo,popfgljunc,res,lab,i},
        Monitor[
            res = Sum[
                popfgl = simMagicFgl[pop, inipopFGL, interferStrength, isObligate, isOogamy,opts];
                popfgldiplo = calFglDiplo[popfgl];
                popfgljunc = calIdentityJunc[popfgldiplo];
                Rest[popfgljunc], {i, sampleSize}], ProgressIndicator[i, {1, sampleSize}]];
        res = res/sampleSize;
        res = Join[inipopFGL[[{1}]], res];
        lab = {"phi(11)^mp","R^m","R^p","rho^mp","J1112^mp", "J1121^mp", "J1122^mp", "J1211^mp", "J1213^mp", "J1222", "J1232^mp"};
        res[[1, -1]] = 
          If[ isOogamy && Length[inipopFGL[[2, -1]]] > 1,
              {StringJoin[#, "_AA"] & /@ lab, StringJoin[#, "_XX/XY"] & /@ lab},
              {StringJoin[#, "_AA"] & /@ lab}
          ];
        res
    ]
    
mergef[list_] :=
    {#[[1, 1]], Total[#[[All, 2]]]} & /@ SplitBy[SortBy[list, First], First]

toAncestryJunc[ancestrydiplo_, nfgl_] :=
    Module[ {fmi, fpi, Rmij, Rpij, rhoij,phiij, Jiiij, Jiiji, Jiijj, Jijii, Jijik, Jijjj, Jijkj, 
      label, rule, phi, changepoint, ls, pos,lab},
        label = {"J1112", "J1121", "J1122", "J1211", "J1213", "J1222", "J1232"};
        rule = {{{x_, x_}, {y_, y_}} :> "J1122", 
            {{x_, x_}, {x_, y_}} :>"J1112", {{x_, x_}, {y_, x_}} :>"J1121", 
            {{x_, y_}, {x_, x_}} :>"J1211", {{x_, y_}, {y_, y_}} :>"J1222", 
            {{x_, y_}, {x_, z_}} :> "J1213", {{x_, y_}, {z_, y_}} :>"J1232", 
            {{_, _}, {_, _}} :> "Wrong!"};
        {phi, changepoint} = ancestrydiplo;
        phiij = SparseArray[Thread[phi[[All, 1]] -> phi[[All, 2]]], {nfgl, nfgl}];
        changepoint = Transpose[{changepoint[[All, 1]] /. rule, changepoint[[All, 1]], changepoint[[All, 2]]}];
        Do[
         ls = Pick[changepoint, changepoint[[All, 1]], lab];
         Switch[lab,
          "J1112",
          pos = ls[[All, 2, 2]];
          Jiiij = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl}],
          "J1121",
          pos = Reverse[#] & /@ ls[[All, 2, 2]];
          Jiiji = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl}],
          "J1122",
          pos = ls[[All, 2, All, 1]];
          Jiijj = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl}],
          "J1211",
          pos = ls[[All, 2, 1]];
          Jijii = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl}],
          "J1213",
          pos = {#[[1, 1]], #[[1, 2]], #[[2, 2]]} & /@ ls[[All, 2]];
          Jijik = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl, nfgl}],
          "J1222",
          pos = ls[[All, 2, 1]];
          Jijjj = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl}],
          "J1232",
          pos = {#[[1, 1]], #[[1, 2]], #[[2, 1]]} & /@ ls[[All, 2]];
          Jijkj = SparseArray[Thread[pos -> ls[[All, 3]]], {nfgl, nfgl, nfgl}]
          ], {lab, label}];
        fmi = Total[phiij, {2}];
        fpi = Total[phiij];
        Rmij = Jiijj + Jiiji + Jijjj + Map[Total, Jijkj, {1}];
        Rpij = Jiijj + Jiiij + Transpose[Jijii] + Total[Jijik];
        rhoij = Rmij + Rpij - Jiijj;
        {fmi, fpi, Rmij, Rpij, rhoij} = SparseArray[#] & /@ {fmi, fpi, Rmij, Rpij, rhoij};
        {fmi, fpi, phiij, Rmij, Rpij, rhoij, Jiiij, Jiiji, Jiijj, Jijii, Jijik, Jijjj, Jijkj}
    ]

Options[simAncestryPrior] = Options[simMagicFgl]
    
simAncestryPrior[pop_, inipopFGL_, interferStrength_, isObligate_,isOogamy_, sampleSize_, opts : OptionsPattern[]] :=
    Module[ {fgl, nfgl, chrLength, runsize = 1000, sizelist, reslist, ls, 
      popfgl, popfgldiplo, geno, res, results, lab,run,ind,chr,i},
        fgl = Union[Flatten[inipopFGL[[2 ;;, -1, 1, All, 1, 2]]]];
        nfgl = Length[fgl];
        If[ fgl != Range[nfgl],
            Print["FGLs are expected to be natural numbers starting from 1!"];
            Abort[]
        ];
        chrLength = inipopFGL[[-1, -1, All, -1, -1, 1]];
        chrLength = chrLength/100.;
        sizelist = Append[Table[runsize, {Quotient[sampleSize, runsize]}],Mod[sampleSize, runsize]];
        If[ Last[sizelist] == 0,
            sizelist = Most[sizelist]
        ];
        reslist = Table[0, {Length[sizelist]}];
        Do[
         PrintTemporary["Iterations =" <> ToString[runsize (run-1)+1]<>" ~ "<>ToString[Min[runsize (run),sampleSize]] <> ". " <> DateString[]];
         Monitor[
         ls = Table[
           popfgl = simMagicFgl[pop, inipopFGL, interferStrength, isObligate, isOogamy, opts];
           popfgldiplo = calFglDiplo[popfgl];
           geno = popfgldiplo[[2 ;;, -1]];
           geno[[All, All, All, 2]], {i,sizelist[[run]]}], ProgressIndicator[i, {1, sizelist[[run]]}]];
         ls = Transpose[ls];
         res = ConstantArray[0, {Length[ls], 2}];
         Do[
          res[[ind, 1]] = Table[Tally[ls[[ind, All, chr, 1]]], {chr, Length[chrLength]}];
          res[[ind, 2]] = Table[Tally[Flatten[Map[Partition[Most[#], 2, 1] &, ls[[ind, All, chr]]],1]], {chr, Length[chrLength]}];
          0, {ind, Length[res]}];
         reslist[[run]] = Transpose[#] & /@ res;
         0, {run, Length[reslist]}];
        results = popfgldiplo[[2 ;;]];
        If[ isOogamy && Length[chrLength] > 1,
            results[[All, -1]] = Table[{Table[mergef[Flatten[reslist[[All, ind, ;; -2, i]], 2]], {i, 2}],
                     Table[mergef[Flatten[reslist[[All, ind, {-1}, i]], 2]], {i,2}]}, {ind, Length[results]}];
            results[[All, -1, 1, All, All, 2]] /= (sampleSize* Total[chrLength[[;; -2]]]);
            results[[All, -1, 2, All, All, 2]] /= (sampleSize* Total[chrLength[[{-1}]]]),
            results[[All, -1]] = Table[{Table[mergef[Flatten[reslist[[All, ind, All, i]], 2]], {i,2}]}, {ind, Length[results]}];
            results[[All, -1, 1, All, All, 2]] /= (sampleSize*Total[chrLength]);
        ];
        results[[All, -1]] = Map[toAncestryJunc[#, nfgl] &, results[[All, -1]], {2}];
        results = Join[inipopFGL[[{1}]], results];
        lab = {"fi^m","fi^p","phiij^mp","Rij^m","Rij^p","rhoij^mp", "Jiiij^mp", "Jiiji^mp", "Jiijj^mp", "Jijii^mp", "Jijik^mp", "Jijjj^mp", "Jijkj^mp"};
        results[[1, -1]] = 
         If[ isOogamy && Length[chrLength] > 1,
             {StringJoin[#, "_AA"] & /@ lab, StringJoin[#, "_XX/XY"] & /@ lab},
             {StringJoin[#, "_AA"] & /@ lab}
         ];
        results
    ]    
    
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicSimulate`*"];

EndPackage[]

(*
PrependTo[$Path, "D:/Chaozhi/WorkspaceCommon/MagicSimulateApp"];
ParallelEvaluate[
  PrependTo[$Path, "D:/Chaozhi/WorkspaceCommon/MagicSimulateApp"]];
Needs["MagicSimulate`"]
ParallelNeeds["MagicSimulate`"]
*)
