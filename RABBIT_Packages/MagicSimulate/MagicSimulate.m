(* Mathematica Package *)

(* Created by the Wolfram Workbench Oct 2, 2013 *)

BeginPackage["MagicSimulate`"]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicSimulate`*"];
ClearAll @@ Names["MagicSimulate`*"];

setfounderFGL::usage = "setfounderFGL  "

pedigreePlot::usage = "pedigreePlot  "

simPedigree::usage = "simPedigree  "

simMagicFgl::usage = "simMagicFgl[founderFGL, mateScheme,chrLength,interferStrength, isObligate,isOogamy] simulates the breeding pedigree, and then drops founder genome lables (FGL)  on the pedigree to obtain FGL list of each individual at all generations. chrLength is a list of  chromosome length for each pair of chromosomes in centiMorgan. founderFGL speices the founder haploid genome labels, a matrix with dimension nfounder x 2. mateScheme is a list of mating scheme. isOogammy denotes if the population is Oogammy. isObligate denote where \nThe format of the return result. The first element of result is the input parameters, the result[[1]]={chrLength, founderFGL, mateScheme, isOogamy, isObligate, interferStrength}. The result[[i+2]] saves the F(i) population, a list of individuals in order of their IDs from 1 to the size of the population. Each individual is reprsented as {generation id,family id,individual id, individual gender, {father id,mother id},fgllist}, where family id is same as individual id in current version. The gender is represented as 1= female, 2=male, and 0=hermaphrodite. The individual's genome is represented as a FGL list. The FGL list is a list of pairs of maternally derived chromosomes and paternally derived chromosomes. The dimensions of the FGL list nnPairChr x 2 because of diploid populations. The FGL of each chromesomes is a list of breakpoints, for example, {{0,10},{47.0, 15},{67.7,2},{100,-1}} denotes the FGL of the segment [0,47.0) is 10, the FGL of the segment [47.0,67.7) is 15, the FGL of the segment [67.7, 100] is 2. The last element {Length of the chromosome (cM), -1} represents the end of the chromosome."

simMagicGeno::usage = "simMagicGeno  "

simMagicErrorGeno::usage = "simMagicErrorGeno  "

toFglGeno::usage = "toFglGeno  "

toFglGenoGrid::usage = "toFglGenoGrid  "

isLastGeneration::usage = "isLastGeneration is an option"

calFglGeno::usage = "calFglGeno  "

calFglGenoGrid::usage = "calFglGenoGrid  "

calFglJunc::usage = "calFglJunc  "

calFglSummary::usage = "calFglSummary  "

pedOrigPrior::usage = "pedOrigPrior  "

Begin["`Private`"]
(* Implementation of the package *)

indexByInterpolation[data_?(OrderedQ[#] && VectorQ[#, NumericQ] &)] :=
    Module[ {ls, f},
        ls = Transpose[{data, Range[Length[data]] - 1}];
        f = Interpolation[ls, InterpolationOrder -> 0];
        Function[x,
          Round[f[x]] + Switch[Depth[x],
               1, Total[Boole[Thread[Rest[data] == x]]],
               2, Total[Boole[Outer[Equal, x, Rest[data]]], {2}]
              ]
         ]
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
    Module[ {gendergroup, nm,nf,temp,i,pairs,newgender},
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
            True,  
            Switch[matingScheme,
             "Sibling",
             temp = Transpose[gendergroup[[All, ;; Min[nm,nf]]]];
             pairs = Riffle[temp,temp],         
             "Pairing", 
             pairs = Transpose[gendergroup[[All, ;; Min[nm,nf]]]],
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
            next = Thread[{-1, -1, 0, pairs,"NA"}],
            Union[genders] == {1,2},
            {pairs,newgenders} = matingParents[matingScheme, genders];
            next = Thread[{-1, -1, newgenders, pairs,"NA"}],
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
          rec = Transpose[{breaks, ch[[indexByInterpolation[ch[[All, 1]]][breaks], 2]]}];
          map = indexByInterpolation[rec[[All, 1]]];
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
    
(*any pedigree member = {generation, individual id, gender, {mother id, father id}, genomes};
generation: non-overlapping, starting from 0 for the founder populations; 
            set generation =-1 for overlapping generations.
gender:1=female,2=male, 0=hermaphrodite;
individual id: natural number starting from 1;
The parents of founders are set to {0,0};founders are always in the beginning. 
All the parents precede their offspring.
*)
(*If isOogammy=True, the founder population consists of 1=female, 2=male, female, male, and so on;
If is Oogammy =False, the gender =0 for all founders*)
simPedigree[founderFGL_, mateScheme_, isOogamy_] :=
    Module[ {founderpop, popped,i},
        founderpop = {-1, -1, 0, {0, 0}, #} & /@ founderFGL;
        founderpop[[All, 2]] = Range[Length[founderFGL]];
        If[ isOogamy,
            Do[founderpop[[i ;; ;; 2, 3]] = i, {i, 1, 2}]
        ];
        popped = FoldList[nextGeneration[#1, #2] &, founderpop, mateScheme];
        Do[{popped[[i, All, 2]], popped[[i + 1, All, 4]]} += 
          popped[[i - 1, -1, 2]], {i, 2, Length[popped] - 1}];
        i = Length[popped];
        popped[[i, All, 2]] += popped[[i - 1, -1, 2]];
        popped[[All, All, 1]] = Range[0, Length[popped] - 1];
        Flatten[popped, 1]
    ]

Options[pedigreePlot] = Options[LayeredGraphPlot];
pedigreePlot[popped_, opts : OptionsPattern[]] :=
    Module[ {data},
        data = DeleteCases[popped[[All, {4, 2}]], {{0, 0}, _}];
        data = Flatten[Thread[#[[1]] -> #[[2]]] & /@ data];
        LayeredGraphPlot[data, Join[Evaluate[FilterRules[{opts}, Options[LayeredGraphPlot]]],
            {VertexLabeling -> True, DirectedEdges -> True}]]
    ]    


Options[simMagicFgl] = Options[calFglSummary] = {isLastGeneration->True}
    (*If isOogamy=True,the last pair of chromosomes are sex chromosomes.*)(*XY for male and XX for female.*)
(*for each linkage group = {maternally derived chromosome, paternally derived chromosome}*)
(*If ith individual is female (gender =1), the paternal dervied chromosome is the X chromosome*)
(*popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy}*)
simMagicFgl[popPed_, chrLength_, interferStrength_,isObligate_,isOogamy_,opts : OptionsPattern[]] :=
    Module[ {popfgl = popPed, nFounder, founderfgl, ls, ls2, gender, i},
        nFounder = Count[popfgl[[All, 4]], {0, 0}];
        founderfgl = popfgl[[;; nFounder, 5]];
        popfgl[[;; nFounder, 5]] = Table[
          ls = {{0, founderfgl[[i, 1]]}, {#, -1}} & /@ chrLength;
          ls2 = {{0, founderfgl[[i, 2]]}, {#, -1}} & /@ chrLength;
          Transpose[{ls, ls2}], {i, nFounder}];
        Do[
         popfgl[[i, 5]] = Transpose[Map[crossover[#, isObligate, interferStrength] &, popfgl[[popfgl[[i, 4]], 5]], {2}]];
         If[ isOogamy&&Length[chrLength]>1,
             gender = popfgl[[i, 3]];
             (*popfgl[[i,4,2]] is the father id*)
             popfgl[[i, -1, -1, 2]] = popfgl[[popfgl[[i, 4, 2]], -1, -1, gender]];
         ], {i, nFounder + 1, Length[popfgl]}];
        If[ TrueQ[OptionValue[isLastGeneration]],
            popfgl = Select[popfgl, First[#] == popfgl[[-1, 1]] &];
        ];
        Join[{{popPed, chrLength, interferStrength,isObligate, isOogamy}},popfgl]
    ]
  
(*any pedigree member = {generation, individual id, gender, {mother id, father id}, genomes}*)
(*popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy,founderFGL, mateScheme}*)
simMagicFgl[founderFGL_, mateScheme_,chrLength_,interferStrength_,isObligate_,isOogamy_,opts : OptionsPattern[]] :=
    Module[ {popped,popfgl},
        popped = simPedigree[founderFGL, mateScheme, isOogamy];
        popfgl = simMagicFgl[popped, chrLength, interferStrength,isObligate, isOogamy,opts];
        popfgl[[1]] = Join[popfgl[[1]],{founderFGL, mateScheme}];
        popfgl
    ] 
    
checkfounderGeno[simSetup_, founderHaplo_] :=
    Module[ {popped,nFounder, isInbredFounder, lengthChr,temp},
        {popped, lengthChr} = simSetup[[;; 2]];
        nFounder = Count[popped[[All, 4]], {0, 0}];
        isInbredFounder = 
          popped[[;; nFounder, 5, 1]] === popped[[;; nFounder, 5, 2]];
        If[ TrueQ[nFounder != Length[founderHaplo] - 3],
            Print["checkfounderGeno: The values for the number of founders between popFgl and founderGeno are not consistent!"];
            Abort[];
        ];
        temp = #[[{1, -1}, 2]] & /@ SplitBy[Transpose[founderHaplo[[2 ;; 3, 2 ;;]]], First];
        If[ Length[temp] != Length[lengthChr],
            Print["checkfounderGeno: The values for the number of founders between popFgl and founderGeno are not consistent!"];
            Abort[];
        ];
        If[ ! (And @@ Thread[temp[[All, 2]] <= lengthChr]) && 
          And @@ Thread[temp[[All, 1]] > 0],
            Print["checkfounderGeno: SNP marker locations are out of the range!"];
            Abort[];
        ];
        True
    ]
fgltodiplo[fgl_, fhaplo_] :=
    Module[ {ls, zz, temp, i,chr,orig},
        Table[
         ls = Table[zz = fgl[[chr, orig]];
                    temp = Table[Select[fhaplo[[chr, All, {1, zz[[i, 2]] + 1}]], zz[[i, 1]] <= #[[1]] < zz[[i + 1, 1]] &], {i, Length[zz] - 1}];
                    Flatten[temp, 1][[All, 2]], {orig, Dimensions[fgl][[2]]}];
         StringJoin @@ # & /@ Transpose[ls], {chr, Length[fgl]}]
    ]
(*to set each founder allele as a string/charactor*)    
simMagicGeno[popfgl_, founderHaplo_] :=
    Module[ {popgeno = popfgl,fhaplo,pos},
        checkfounderGeno[popfgl[[1]], founderHaplo];
        fhaplo = #[[All, 2 ;;]] & /@ SplitBy[Transpose[founderHaplo[[2 ;;, 2 ;;]]], First];
        popgeno[[2 ;;, -1]] = fgltodiplo[#, fhaplo] & /@ popgeno[[2 ;;, -1]];
        (*gender\[Equal]2===male,only X haplotype*)
        pos = Flatten[Position[popgeno[[2 ;;, 3]], 2]] + 1;
        popgeno[[pos, -1]] = Map[StringTake[#, 1] &, popgeno[[pos, -1]], {2}];
        popgeno
    ]
    
simMagicErrorGeno[popgeno_, eps_] :=
    Module[ {rules, obsgeno, geno, err, pos, obs,ii,ch},
        (*haplotypes "1","2","N" reference to X of male XY sex chromosomes*)
        rules = {{"1", "2", "N", "11", "12", "1N", "21", "22", "2N", "N1", "N2", "NN"},
                   {"2", "1", "N", "21", "22", "2N", "11", "12", "1N", "N1", "N2", "NN"},
                   {"1", "2", "N", "12", "11", "1N", "22", "21", "2N", "N2", "N1", "NN"},
                   {"2", "1", "N", "22", "21", "2N", "12", "11", "1N", "N2", "N1", "NN"}};
        (*The rules[[ii]] for ii =1,2,3,4, refer to 0, 1, 1, 2 allelic errors respectively.*)
        rules = Thread[rules[[1]] -> #] & /@ rules;
        obsgeno = popgeno;
        Do[
            geno = obsgeno[[2 ;;, -1, ch]];
            err = RandomChoice[{(1 - eps)^2, eps (1 - eps), (1 - eps) eps, eps^2} -> Range[4], Dimensions[geno]];
            Do[
                pos = Position[err, ii];
                obs = Extract[geno, pos] /. rules[[ii]];
                geno = ReplacePart[geno, Thread[pos -> obs]], {ii, 2, 4}];
            obsgeno[[2 ;;, -1, ch]] = geno, {ch, Length[obsgeno[[2, -1]]]}];
        obsgeno
    ]    
          
toFglGeno[fglhap_] :=
    Module[ {fun},
        fun[z_] :=
            Module[ {junc, temp,l},
                junc = Union[Flatten[z[[All, All, 1]]]];
                temp = Transpose[Table[z[[l, indexByInterpolation[z[[l, All, 1]]][junc], 2]], {l,Length[z]}]];
                Transpose[{junc, temp}]
            ];
        Map[fun, fglhap, {Depth[fglhap] - 4}]
    ]

toFglGenoGrid[fglgeno_, grid_] :=
    Module[ {fun},
        fun[z_] :=
            Module[ {temp,ch},
                Table[
                 temp = 
                  z[[ch, indexByInterpolation[z[[ch, All, 1]]][grid[[ch]]], 2]];
                 Transpose[{grid[[ch]], temp}], {ch, Length[grid]}]
            ];
        Map[fun, fglgeno, {Depth[fglgeno] - 5}]
    ]    

toJunction[diplo_] :=
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

(*popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy}*)
(*or popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy,founderFGL, mateScheme}*)
calFglGeno[popfgl_,obligateRescale_:True] :=
    Module[ {popfglgeno = popfgl,geno,chrlen,isObligate,factor,ch},
        {chrlen, isObligate} = popfglgeno[[1, {2, 4}]];
        geno = toFglGeno[popfglgeno[[2 ;;, -1]]];
        (*correcting for obligate crossover*)
        If[ isObligate&&obligateRescale,
            (*transform from centiMorgan into Morgan*)
            chrlen/=100;
            factor = 1. + Exp[-chrlen]/(chrlen);
            Do[geno[[All, ch, All, 1]] *= factor[[ch]], {ch, Length[factor]}];
        ];
        popfglgeno[[2 ;;, -1]] = geno;
        popfglgeno
    ] 

calFglGenoGrid[popfglgeno_, grid_] :=
    Module[ {popgrid = popfglgeno, chrlen},
        chrlen = popgrid[[1, 2]];
        If[ Length[grid] != Length[chrlen],
            Print["The #chromosomes are different between grid and popfglgeno!"];
            Abort[];
        ];
        If[ ! ((And @@ Thread[grid[[All, -1]] < chrlen]) && (And @@Thread[grid[[All, 1]] >= 0])),
            Print["The locations of some grid points are out of range!"];
            Abort[];
        ];
        popgrid[[2 ;;, -1]] = toFglGenoGrid[popgrid[[2 ;;, -1]], grid];
        popgrid
    ]
    
calFglJunc[popfglgeno_,opts : OptionsPattern[]] :=
    Module[ {pop = popfglgeno,res, isOogamy, chrlen, geno, ibd, diplo, junc, ls},
        isOogamy = pop[[1, 5]];
        chrlen = pop[[2, -1, All, -1, 1]]/100.;
        geno = pop[[2 ;;, -1]];
        geno = geno[[All, All, All, 2]];
        ibd = geno[[All, All, 1]];
        ibd = Boole[Map[SameQ @@ # &, ibd, {2}]];
        diplo = Map[Partition[Most[#], 2, 1] &, geno, {2}];
        (*{"J1112", "J1121", "J1122", "J1211", "J1213", "J1222", "J1232"}*)
        junc = toJunction[diplo];
        If[ isOogamy && Length[chrlen] > 1,
            ls = Join[Transpose[{Mean[Transpose[ibd[[All, ;; -2]]]]}], 
              Total[junc[[All, ;; -2]], {2}]/Total[Most[chrlen]], 2];
            ls = Transpose[{ls, Join[Transpose[{ibd[[All, -1]]}], junc[[All, -1]]/Last[chrlen],2]}],
            ls = Partition[Join[Transpose[{Mean[Transpose[ibd]]}], Total[junc, {2}]/Total[chrlen], 2], 1];
        ];
        res = pop;
        res[[2 ;;, -1]] = ls;
        res
    ]    
    
(*popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy}*)
(*or popfgl[[1]] ={popPed, chrLength, interferStrength,isObligate, isOogamy,founderFGL, mateScheme}*)
(*gender:1=female,2=male, 0=hermaphrodite;*)
calFglSummary[popfgljunc_] :=
    Module[ {pop = popfgljunc,resAA, tls, res, ls, temp, pos,gender},
        resAA = (Total[#]/Length[#]) & /@ SplitBy[Transpose[{pop[[2 ;;, 1]],pop[[2 ;;, -1, 1]]}], First];
        tls = resAA[[All, 1]];
        If[ pop[[1, 5]] && Length[pop[[2, -1]]] == 2,
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

calFglSummary[popPed_, chrLength_,interferStrength_,isObligate_,isOogamy_,sampleSize_, opts : OptionsPattern[]] :=
    Module[ {popfgl,popfglgeno,popfgljunc,summary,res,i},
        Monitor[
            res = Sum[
                popfgl = simMagicFgl[popPed, chrLength, interferStrength, isObligate, isOogamy,opts];
                popfglgeno = calFglGeno[popfgl,True];
                popfgljunc = calFglJunc[popfglgeno];
                summary = calFglSummary[popfgljunc];
                Rest[summary], {i, sampleSize}], ProgressIndicator[i, {1, sampleSize}]];
        Join[{Append[First[summary],sampleSize]}, res/sampleSize]
    ]
    
calFglSummary[founderFGL_, mateScheme_, chrLength_,interferStrength_,isObligate_,isOogamy_,sampleSize_, opts : OptionsPattern[]] :=
    Module[ {popfgl,popfglgeno,popfgljunc,summary,res,i},
        Monitor[
            res = Sum[
                popfgl = simMagicFgl[founderFGL, mateScheme, chrLength, interferStrength, isObligate, isOogamy,opts];
                popfglgeno = calFglGeno[popfgl,True];
                popfgljunc = calFglJunc[popfglgeno];
                summary = calFglSummary[popfgljunc];
                Rest[summary], {i, sampleSize}], ProgressIndicator[i, {1, sampleSize}]];
        Join[{Append[First[summary],sampleSize]}, res/sampleSize]
    ]    
       
pedOrigPrior[popPed_, chrLength_,interferStrength_,isObligate_,isOogamy_,sampleSize_] :=
    Module[ {summary,i},
        summary = calFglSummary[popPed, chrLength, interferStrength, isObligate, isOogamy, sampleSize, isLastGeneration -> True];
        summary = Rest[summary[[-1]]];
        (*{{"finb","J1112","J1121","J1122","J1211","J1213","J1222","J1232"},..} to 
          {{finbred,j1122,j1211,j1213,j1222,j1232},{finbred_mp,j1122mp,j1211mp,j1213mp,j1222mp,j1232mp}}*)
        Table[{summary[[i, 1]], summary[[i, 4]],Mean[summary[[i, {2, 5}]]], summary[[i, 6]], 
        	                                     Mean[summary[[i, {3, 7}]]], summary[[i, 8]]},{i,Min[2,Length[summary]]}]
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
