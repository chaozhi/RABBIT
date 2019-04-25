(* Mathematica Package *)

(* Created by the Wolfram Workbench Feb 22, 2016 *)

BeginPackage["MagicDefinition`"]
(* Exported symbols added here with SymbolName::usage *) 

linkageGroupSet::usage = "linkageGroupSet is an option to specify the list of linkage groups. "

founderAllelicError::usage = "founderAllelicError is an option to specify the allelic error probability in founders. "

offspringAllelicError::usage = "offspringAllelicError is an option to specify the allelic error probabilityin offspring. "

isFounderInbred::usage = "isFounderInbred is an option to specify whether the founders are completely inbred."

isFounderAllelicDepth::usage = "isFounderAllelicDepth is an option to specify whether founder data are allelic depths of genotyping-by-sequencing (GBS) (True) or called genotypes from GBS or SNP Array (default False)."

sequenceDataOption::usage = "sequenceDataOption is an option to specify option rules for allelic depths of sequencing data. "

isOffspringAllelicDepth::usage = "isOffspringAllelicDepth is an option to specify whether offspring data are allelic depths of genotyping-by-sequencing (GBS) (True) or called genotypes from GBS or SNP Array (default False)."

priorFounderCallThreshold::usage = "priorFounderCallThreshold is an option to specify the threshould for prior calling of missing founder genotypes. Before founder genotype imputation, single locus calling of founder genotypes is performed if the posterior probability of the true genotype is greater than the threshould. "

sampleSize::usage = "sampleSize is an option to specify the number of samples. "

outputFileID::usage = "outputFileID is an option to specify the stem of output filename."

isPrintTimeElapsed::usage = "isPrintTimeElapsed is an option to secify whether to print and some details and time elapsed."

isRunInParallel::usage = "isRunInParallel is an option to specify whether to compute in parallel."

detectingThreshold::usage = "detectingThreshold is an option to specify a correction threshold. An observed genotype is corrected only if the posterior probability of the true genotype is greater than the threshold and is greater than the posterior probability of the observed genotype."

imputingThreshold::usage = "imputingThreshold is an option to specify an imputing threshold. A missing genotype is imputed only if its posterior probabaility is greater than the threshold. "

minPhredQualScore::usage = "minPhredQualScore is an option to specify base-calling error probability P=10^(-minPhredQualScore/10)."

origGenotype::usage = "origGenotype[nFgl] gives {{genotypes, diplotypes}, {geno2diplo, diplo2geno}}.  The genotypes and diplotypes are the definitions where the inbred founders are labeled by natural numbers starting from 1,geno2diplo gives the mapping from genotypes to dipotypes, and diplo2geno gives the mapping from diplotypes to genotypes."

intervalIndexFunction::usage = "intervalIndexFunction[list] is a function for searching the grid of a one-dimensional ordered list; intervalIndexFunction[list][x] returns index i so that list[[i-1]]<=x<list[[i]]; intervalIndexFunction[list][{x1,x2..}] return a list of indices."

maxIndexPair::usage = "maxIndexPair[list] returns {max,indices of max} for each element at the second-to-last level of the list, where indices is a list since max may occur multiple times."

toDelimitedString::usage = "toPipeDelimitedString[list,delimiter] first transforms into string for each element at the last level, then tranforms into delimited-string for eacb element at the second-to-last level of the list. By default, delimiter=\"|\"."

csvExport::usage = "csvExport[file,table] exports table into comma-delimited file, where the depth of table must be 3. Note loss of precision for numeric number."

csvExport2::usage = "csvExport2 is similar to csvExport but for a large table"

(*readTable::usage = "readTable[file,format,chunkSize:1000] performs memory-efficient import of tabular data in file in format and default chuncksize of 1000."*)

readTable::usage = "readTable  "

readDictTable::usage = "readDictTable[file,format] read multiple tables that are seperated by one line. The line consists of three elements: table separator, asssoication key, and description. It return an association between keys and tables."

kroneckerMatrixSum::usage = "kroneckerMatrixSum  "

removeDownValues::usage = "removeDownValues  "

myParallelNeeds::usage = "myParallelNeeds  "

listErrorPlot::usage = "listErrorPlot  "

crossTabulate::usage = "crossTabulate returns a cross table from two-column input data."

independenceGTest::usage = "independenceGTest  "

randomIntegerPartition::usage = "randomIntegerPartition  "

renameDuplicates::usage = "renameDuplicates  "

notebookEvaluate::usage = "notebookEvaluate  "




Begin["`Private`"]
(* Implementation of the package *)
    
(*diplotypes refer to a single locus, phased genotypes*)
(*geno2diplo[[i]] refer to the index of the ith genotype, in the diplotypes*)
origGenotype[nFgl_] :=
    Module[ {genotypes,diplotypes,geno2diplo,diplo2geno, i,j},
        diplotypes = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        (*lower triangular matrix*)
        genotypes = Flatten[Table[{i, j}, {i, nFgl}, {j, i}], 1];
        diplo2geno = Flatten[Position[genotypes, #, {1}, 1] & /@ (Sort[#, Greater] & /@diplotypes)];
        geno2diplo = {#[[All, 1]], #[[1, 2]]} & /@ 
            SplitBy[SortBy[Transpose[{Range[Length[diplotypes]], diplo2geno}], Last],Last];
        geno2diplo = geno2diplo[[All, 1]];
        {{genotypes,diplotypes},{geno2diplo,diplo2geno}}
    ]

intervalIndexFunction[list_?(VectorQ[#, NumericQ]&&OrderedQ[#,Less] &)] :=
    Module[ {ls, f},
        ls = Transpose[{list, Range[Length[list]] - 1}];
        f = Interpolation[ls, InterpolationOrder -> 0];
        Function[x,
          Round[f[x]] + Switch[Depth[x],
               1, Total[Boole[Thread[Rest[list] == x]]],
               2, Total[Boole[Outer[Equal, x, Rest[list]]], {2}]
              ]
         ]
    ]

maxIndexPair[list_List] :=
    Map[{Max[#], Flatten[Position[#, Max[#], {1}, Heads -> False]]} &, list, {-2}]
    
toDelimitedString[list_List,delimiter_String:"|"] :=
    Map[StringJoin @@ Riffle[#, delimiter] &, Map[ToString, list, {-1}], {-2}]    

csvExport2[file_String, table_] :=
    Module[ {ls,i},
        ls = Replace[table,x_?NumericQ :> ToString[x, FormatType -> InputForm], {2}];
        ls = Table[StringReplace[TextString[ls[[i]]], {"{" | "}" -> "", ", " -> ","}],{i,Length[ls]}];
        Export[file, ls, "List"]
    ]
    
csvExport[file_String, table_,nrow_Integer:100] :=
    Module[ {ls,ii,iils},
        ls = Table[0,{Length[table]}];
        iils = Partition[Range[Length[ls]], UpTo[nrow]];
        Do[
            ls[[ii]] = Replace[table[[ii]],x_?NumericQ :> ToString[x, FormatType -> InputForm], {2}];
            ls[[ii]] = StringReplace[TextString[#], {"{" | "}" -> "", ", " -> ","}]&/@ls[[ii]],{ii,iils}];
        Export[file, ls, "List"]
    ]    
        
(*http://stackoverflow.com/questions/7525782/import-big-files-arrays-with-mathematica*)
(*does not work in v11.3, works in v11.0.0*)
(*readTable[file_String?FileExistsQ, format_String, chunkSize_Integer:1000] :=
    Module[ {stream, dataChunk, result, linkedList, add},
        SetAttributes[linkedList, HoldAllComplete];
        add[ll_, value_] :=
            linkedList[ll, value];
        stream = StringToStream[Import[file, "String"]];
        Internal`WithLocalSettings[Null,(*main code*)
            result = linkedList[];
            While[dataChunk =!= {}, 
             dataChunk = ImportString[StringJoin[Riffle[ReadList[stream, "String", chunkSize], "\n"]], format];
             result = add[result, dataChunk];
            ];
            result = Flatten[result, Infinity, linkedList],(*clean-up*)
         Close[stream]
         ];
        Join @@ result
    ]*)
    
readTable[file_String?FileExistsQ, format_String]:=Import[file,format]    
    
readDictTable[file_String?FileExistsQ, format_String] :=
    Module[ {data, sep},
        data = Import[file, format];
        sep = data[[1, 1]];
        data = Partition[Split[data, #1[[1]] != sep && #2[[1]] != sep &], 2];
        Association[#[[1, 1, 2]] -> #[[2]] & /@ data]
    ] 

(*kroneckerMatrixSum[m1_?MatrixQ, m2_?MatrixQ] :=
    KroneckerProduct[m1, IdentityMatrix[Length[m2]]] + 
     KroneckerProduct[IdentityMatrix[Length[m1]], m2]*)
     
kroneckerMatrixSum[a_, b_] /; MatrixQ[a] && MatrixQ[b] :=
    Catch@Module[ {n, p, m, q},
              {n, p} = Dimensions[a];
              {m, q} = Dimensions[b];
              If[ n != p || m != q,
                  Throw[$Failed]
              ];
              KroneckerProduct[a, IdentityMatrix[m]] +
               KroneckerProduct[IdentityMatrix[n], b]
          ]     

(*http://mathematica.stackexchange.com/questions/19536/how-to-clear-parts-of-a-memoized-function*)
SetAttributes[removeDownValues, HoldAllComplete];
removeDownValues[p : f_[___]] :=
    DownValues[f] = DeleteCases[DownValues[f, Sort -> False], 
      HoldPattern[Verbatim[HoldPattern][p] :> _]];

Options[myParallelNeeds] = {
  Path:>$Path
  }
  
myParallelNeeds[context_String, opts : OptionsPattern[]] :=
    myParallelNeeds[{context}, opts]
myParallelNeeds[contextlist_List, OptionsPattern[]] :=
    Module[ {path},
        SetSharedVariable[path];
        {path} = OptionValue@{Path};
        path = DeleteDuplicates[Join[Flatten[{path}],{Directory[]}]];
        $Path = DeleteDuplicates[Join[$Path, path]];
        ParallelEvaluate[$Path = DeleteDuplicates[Join[$Path, path]]];
        ParallelNeeds[#] & /@ contextlist;
    ]    


(*data ={{x1,ylow1,ymid1,ytop1},{x2,ylow2,ymid2,ytop2},...}*)
listErrorPlot[data_, xlogscale_, ylogscale_, opts : OptionsPattern[]] :=
    Module[ {ls, fun,i},
        ls = SortBy[data, First];
        ls = Table[ls[[All, {1, i}]], {i, 2, Dimensions[data][[2]]}];
        fun = Switch[{xlogscale, ylogscale},
          {False, False}, ListPlot,
          {True, False}, ListLogLinearPlot,
          {False, True}, ListLogPlot,
          {True, True}, ListLogLogPlot,
          _,
          Print["Wrong {xlogscale,ylogscale} =", {xlogscale, ylogscale}, 
           ". xlogscale (ylogscale) must be True or False."];
          Abort[]
          ];
        fun[ls, opts, Filling -> {1 -> {3}}]
    ]
    
crossTabulate[data_?MatrixQ] :=
    Module[ {ls, rule},
        If[ Dimensions[data][[2]] != 2,
            Print["crossTabulate: data matrix with two columns are expected!"];
            Abort[]
        ];
        ls = Tally[data];
        rule = Union[#] & /@ Transpose[ls[[All, 1]]];
        rule = Thread[# -> Range[Length[#]]] & /@ rule;
        ls[[All, 1]] = Transpose[MapThread[#1 /. #2 &, {Transpose[ls[[All, 1]]], rule}]];
        ls = Rule @@ # & /@ ls;
        Association["CrossTab" -> SparseArray[ls], "RowNames" -> rule[[1, All, 1]], "ColumnNames" -> rule[[2, All, 1]]]
    ]    

(*data: list of cross tables, return: list of G-test for each table; 
each G-test ={outcome of G^2, degree of freedom, -Log[10,P-value], Lodscore} *)
independenceGTest[data_?MatrixQ] :=
    First[independenceGTest[{data}]]
independenceGTest[data_?ListQ] :=
    Module[ {count = data, rowsum, colsum, nn, expect, temp, g2, df, 
      pvalue, neglogpvalue, lod, pos, res},
        res = ConstantArray[0, {Length[data], 4}];
        rowsum = Total[count, {3}];
        colsum = Total[count, {2}];
        df = ((Length[#] & /@ rowsum) - 1) ((Length[#] & /@ colsum) - 1);
        pos = Intersection @@ (Flatten[Position[(Length[#] & /@ #) - 1, _?Positive]] & /@ {rowsum,colsum});
        If[ pos=!={},
            {rowsum, colsum, df, count} = {rowsum, colsum, df, count}[[All, pos]];
            nn = Total[rowsum, {2}];
            expect = MapThread[KroneckerProduct[#1, #2] &, {rowsum, colsum}]/N[nn];
            count = Flatten[#] & /@ count;
            temp = Log[count /. {0 -> 1}] - Log[Flatten[#] & /@ expect];
            g2 = 2 MapThread[#1.#2 &, {count, temp}];
            pvalue = MapThread[SurvivalFunction[ChiSquareDistribution[#1], #2] &, {df, g2}];
            neglogpvalue = -Log[10, pvalue];
            lod = InverseSurvivalFunction[ChiSquareDistribution[1], #1] & /@pvalue;
            lod *= Log[10, E]/2;
            res[[pos]] = Transpose[{g2, df, neglogpvalue, lod}];
        ];
        res
    ]
        
randomIntegerPartition[n_, dist_] :=
    Module[ {size, res, sum, pos},
        size = 2 Ceiling[n/Mean[dist]];
        res = {};
        While[True,
         res = Join[res, RandomInteger[dist, size]];
         If[ Total[res] > n,
             Break[]
         ]
         ];
        sum = Accumulate[res];
        pos = Position[Sign[sum - n], 0 | 1, {1}, 1, Heads -> False][[1, 1]];
        res = Take[res, pos];
        res[[-1]] -= sum[[pos]] - n;
        res
    ]
        
renameDuplicates[list_?(VectorQ[#, StringQ] &)] :=
    Module[ {dup, res = list, pos},
        If[ DuplicateFreeQ[list],
            list,
            dup = Select[Tally[list], #[[2]] > 1 &];
            Do[
             pos = Flatten[Position[list, dup[[i, 1]]]];
             res[[pos]] = 
              dup[[i, 1]] <> "_" <> ToString[#] & /@ Range[dup[[i, 2]]];
             0, {i, Length[dup]}];
            res
        ]
    ]
            
notebookEvaluate[nbfilename_String] :=
    Module[ {nb},
        nb = NotebookOpen[FindFile[nbfilename]];
        NotebookEvaluate[nb];
        NotebookClose[nb]
    ]            
            
End[]

EndPackage[]

