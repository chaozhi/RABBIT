(* Mathematica Package *)

(* Created by the Wolfram Workbench Aug 24, 2013 *)

BeginPackage["ContinuousTimeHmm`"]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["ContinuousTimeHmm`*"];
ClearAll @@ Names["ContinuousTimeHmm`*"];

CtViterbi::usage = "CtViterbi "

CtForward::usage = "CtForward "

CtBackward::usage = "CtBackward "

CtPosteriorProb::usage = "CtPosteriorProb "

CtLogLiklihood::usage = "CtLogLiklihood "

CtPathLogProb::usage = "CtPathLogProb "

CtPosteriorDecoding::usage = "CtPosteriorDecoding "

CtViterbiScale::usage = "CtViterbiScale  "

CtPathSampling::usage = "CtPathSampling  "

CtDiscreteFormat::usage = "CtDiscreteFormat  "

CtJumpFormat::usage = "CtJumpFormat  "

CtStringFormat::usage = "CtStringFormat  "

CtProbPlot::usage = "CtProbPlot  "

CtPathPlot::usage = "CtPathPlot  "

CtSeqPlot::usage = "CtSeqPlot  "


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
      

CtViterbiScale[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {logdataProb = Log[dataProbSeq/.{0.->0}], nSeq = Length[dataProbSeq], logtransProb, vlogProb, vIndex, 
      temp, ii, vPathLogProb, vPath, t},
        vlogProb = vIndex = Table[0, {nSeq}, {Length[startProb]}];
        vlogProb[[1]] = Log[startProb/.{0.->0}] + logdataProb[[1]];
        Do[
         logtransProb = Log[tranProbSeq[[t-1]]/.{0.->0}];
         temp = ((vlogProb[[t - 1]] + #) & /@ Transpose[logtransProb]) + logdataProb[[t]];
         vlogProb[[t]] = Max[#] & /@ temp;
         vIndex[[t]] = Flatten[MapThread[Position[#1, #2, {1}, 1, Heads -> False] &, {temp, vlogProb[[t]]}]], {t, 2, nSeq}];
        vPathLogProb = Max[vlogProb[[-1]]];
        ii = Position[vlogProb[[-1]], vPathLogProb, {1}, 1, Heads -> False][[1, 1]];
        vPath = Reverse[NestList[{#[[1]] - 1, vIndex[[#[[1]], #[[2]]]]} &, {nSeq, ii}, nSeq - 1][[All, 2]]];
        {vPathLogProb, vPath}
    ]  


CtViterbi[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {nSeq = Length[dataProbSeq],vProb, vIndex, 
      temp, ii,transProb,vPathProb, vPath, t},
        vProb = vIndex = Table[0, {nSeq}, {Length[startProb]}];
        vProb[[1]] = startProb dataProbSeq[[1]];
        Do[
         transProb = tranProbSeq[[t-1]];
         temp = ((vProb[[t - 1]] #) & /@ Transpose[transProb]) dataProbSeq[[t]];
         vProb[[t]] = Max[#] & /@ temp;
         vIndex[[t]] = Flatten[MapThread[Position[#1, #2, {1}, 1, Heads -> False] &, {temp, vProb[[t]]}]], {t, 2, nSeq}];
        vPathProb = Max[vProb[[-1]]];
        ii = Position[vProb[[-1]], vPathProb, {1}, 1, Heads -> False][[1, 1]];
        vPath = Reverse[NestList[{#[[1]] - 1, vIndex[[#[[1]], #[[2]]]]} &, {nSeq, ii}, nSeq - 1][[All, 2]]];
        {Log[vPathProb], vPath}
    ]    

CtViterbi2[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {nSeq = Length[dataProbSeq],vProb, vIndex, 
      temp, ii,transProb,vPathProb, vPath, t},
        vProb = vIndex = Table[0, {nSeq}, {Length[startProb]}];
        vProb[[1]] = startProb dataProbSeq[[1]];
        Do[
         transProb = tranProbSeq[[t-1]];
         temp = ((vProb[[t - 1]] #) & /@ Transpose[transProb]) dataProbSeq[[t]];
         vProb[[t]] = Max[#] & /@ temp;
         vIndex[[t]] = Flatten[MapThread[RandomChoice[Position[Thread[#1 == #2], True, {1},Heads -> False]]&,{temp,vProb[[t]]}]], {t, 2, nSeq}];
        vPathProb = Max[vProb[[-1]]];
        (*ii = Position[vProb[[-1]], vPathProb, {1}, 1, Heads -> False][[1, 1]];*)
        ii = Position[Thread[vProb[[-1]] == vPathProb], True, {1},Heads -> False][[1,1]];
        vPath = Reverse[NestList[{#[[1]] - 1, vIndex[[#[[1]], #[[2]]]]} &, {nSeq, ii}, nSeq - 1][[All, 2]]];
        {Log[vPathProb], vPath}
    ]     

CtForward[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {nSeq = Length[dataProbSeq],forwardProb, forwardScale, t},
        forwardProb = forwardScale = Table[0, {nSeq}];
        forwardProb[[1]] = startProb dataProbSeq[[1]];
        forwardScale[[1]] = Total[forwardProb[[1]]];
        forwardProb[[1]] /= forwardScale[[1]];
        Do[
         forwardProb[[t]] = (forwardProb[[t - 1]].tranProbSeq[[t-1]]) dataProbSeq[[t]];
         forwardScale[[t]] = Total[forwardProb[[t]]];
         forwardProb[[t]] /= forwardScale[[t]], {t, 2, nSeq}];
        {forwardProb, forwardScale}
    ]
    
(*backwardProb may be \[GreaterEqual]1 because it is scaled by the fowardScale*)
CtBackward[tranProbSeq_, dataProbSeq_, forwardScale_] :=
    Module[ {nSeq = Length[dataProbSeq],backwardProb,t},
        backwardProb = Table[0, {nSeq}];
        backwardProb[[-1]] = Table[1, {Length[tranProbSeq[[-1]]]}];
        Do[
         backwardProb[[t]] = tranProbSeq[[t]].(dataProbSeq[[t+1]] backwardProb[[t + 1]]);
         backwardProb[[t]] /= forwardScale[[t + 1]], {t, nSeq - 1, 1, -1}];
        backwardProb
    ]
        
(*
complileCtForward = 
  Compile[{{startProb, _Real, 1}, {tranProbSeq, _Real,3}, {dataProbSeq, _Real, 2}},
   Module[ {nSeq = Length[dataProbSeq], t},
       forwardScale = ConstantArray[0, nSeq];
       forwardProb = ConstantArray[0, {nSeq, Length[startProb]}];
       forwardProb[[1]] = startProb dataProbSeq[[1]];
       forwardScale[[1]] = Total[forwardProb[[1]]];
       forwardProb[[1]] = forwardProb[[1]]/forwardScale[[1]];
       Do[
          forwardProb[[t]] = (forwardProb[[t-1]].tranProbSeq[[t - 1]]) dataProbSeq[[t]];
          forwardScale[[t]] = Total[forwardProb[[t]]];
          forwardProb[[t]] = forwardProb[[t]]/forwardScale[[t]], {t, 2, nSeq}];
       MapThread[Join[{#1}, #2] &, {forwardScale, forwardProb}]
   ], {{forwardProb, _Real, 2}, {forwardScale, _Real, 1}}
   ]
CtForward2[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {forwardScaleProb},
        forwardScaleProb = 
         complileCtForward[startProb, tranProbSeq, dataProbSeq];
        {forwardScaleProb[[All, 2 ;;]], forwardScaleProb[[All, 1]]}
    ]
      
CtBackward2 = 
 Compile[{{tranProbSeq, _Real, 3}, {dataProbSeq, _Real, 2}, {forwardScale, _Real, 1}},
  Module[ {nSeq = Length[dataProbSeq], transitionProb, t},
      backwardProb = ConstantArray[0, {nSeq, Length[tranProbSeq[[-1]]]}];
      backwardProb[[-1]] = Table[1, {Length[tranProbSeq[[-1]]]}];
      Do[transitionProb = tranProbSeq[[t]];
         backwardProb[[t]] = transitionProb.(dataProbSeq[[t + 1]] backwardProb[[t + 1]]);
         backwardProb[[t]] = backwardProb[[t]]/forwardScale[[t + 1]], {t,nSeq - 1, 1, -1}];
      backwardProb
  ], {{backwardProb, _Real, 2}}
  ]    
*)  
    
CtPosteriorProb[forwardProb_, backwardProb_] :=
    forwardProb backwardProb

CtLogLiklihood[forwardScale_] :=
    Total[Log[forwardScale]]

CtLogLiklihood[startProb_, tranProbSeq_, dataProbSeq_] :=
    Total[Log[Last[CtForward[startProb, tranProbSeq, dataProbSeq]]]]
    
CtPathLogProb[startProb_, tranProbSeq_, dataProbSeq_,  path_] :=
    Module[ {logl, transitionProb,t},
        logl = Log[startProb[[First[path]]]] + Log[dataProbSeq[[1,path[[1]]]]];
        Do[
         transitionProb = tranProbSeq[[t-1]];
         logl += Log[transitionProb[[path[[t - 1]], path[[t]]]]] + 
           Log[dataProbSeq[[t,path[[t]]]]], {t, 2, Length[dataProbSeq]}];
        logl
    ] 

(*CtPosteriorPath[posteriorProb_] :=
    (Position[#, Max[#], {1}, 1, Heads -> False] & /@posteriorProb)[[All, 1, 1]]*)  
        
CtPosteriorDecoding[startProb_, tranProbSeq_, dataProbSeq_] :=
    Module[ {forwardProb, forwardScale, backwardProb, logLikelihood, posteriorProb},
        {forwardProb, forwardScale} = CtForward[startProb, tranProbSeq, dataProbSeq];
        backwardProb = CtBackward[tranProbSeq, dataProbSeq, forwardScale];
        logLikelihood = CtLogLiklihood[forwardScale];
        posteriorProb = CtPosteriorProb[forwardProb, backwardProb];
        {logLikelihood, posteriorProb}
    ]    
    
    
CtPathSampling[startProb_, tranProbSeq_, dataProbSeq_, size_:1] :=
    Module[ {nSeq = Length[dataProbSeq], forwardProb, forwardScale, states = Range[Length[startProb]], transitionProb,weights,res,t},
        {forwardProb, forwardScale} = CtForward[startProb, tranProbSeq, dataProbSeq];
        res = ConstantArray[0, {size, nSeq}];
        res[[All, -1]] = RandomChoice[forwardProb[[-1]] -> states, size];
        Do[
            transitionProb = tranProbSeq[[t]];
            weights = Abs[Transpose[forwardProb[[t]] transitionProb[[All, res[[All, t + 1]]]]]];
            res[[All, t]] = RandomChoice[# -> states] & /@ weights, {t,nSeq - 1, 1, -1}];
        {CtLogLiklihood[forwardScale], res}
    ]
      
CtDiscreteFormat[x_List] :=
    Flatten[ConstantArray[#[[2]], #[[1]]] & /@ 
      Transpose[{Differences[x[[All, 1]]], x[[;; -2, 2]]}], 1]      
          
CtJumpFormat[x_List] :=
    Module[ {ls},
        ls = {Length[#], #[[1]]} & /@ Split[x];
        ls[[All, 1]] = Prepend[Most[Accumulate[ls[[All, 1]]]] + 1, 1];
        Append[ls, {Length[x] + 1, Null}]
    ]
      
CtJumpFormat[x_String, delimiter_String: "-"] :=
    Partition[Append[ToExpression[StringSplit[x, delimiter]], Null], 2]
 
CtStringFormat[jumpPath_List, delimiter_String: "-"] :=
    StringJoin[Riffle[ToString[#] & /@ Most[Flatten[jumpPath]], delimiter]]      
        
CtProbPlot[timeseq_?(VectorQ[#, NumericQ] &),array_?(ArrayQ[#, 2, NumericQ] &), rectcolor_] :=
    Module[ {xxls = timeseq, yyls, ydim = Dimensions[array][[2]], matrix = array, rect, ls,i,j},
        matrix = Append[matrix, matrix[[-1]]];
        xxls = Append[xxls, Last[xxls]];
        yyls = Range[ydim];
        rect = Table[
          ls = SplitBy[Transpose[{xxls, matrix[[All, j]]}], Last][[All, 1]];
          Table[{Lighter[rectcolor, 1 - ls[[i, 2]]], 
              Rectangle[{ls[[i, 1]], yyls[[j]] - 0.5}, {ls[[i + 1, 1]],yyls[[j]] + 0.5}]}, {i, Length[ls] - 1}], {j, ydim}];
        Graphics[rect]
    ] /; Length[timeseq] == Length[array]

CtProbPlot[timeseqlist_?(VectorQ[#, ListQ] &), arraylist_?(VectorQ[#, ArrayQ] &), colorlist_?VectorQ, boundarystyle___] :=
    Module[ {xxls = timeseqlist, len, ii, g, g2},
        len = xxls[[;; -2, -1]];
        xxls[[2 ;;]] += Accumulate[len];
        g = Table[CtProbPlot[xxls[[ii]], arraylist[[ii]], colorlist[[ii]]], {ii, Length[xxls]}];
        g2 = Graphics[{boundarystyle, 
             Line[Thread[{#, {0, Dimensions[First[arraylist]][[2]]}}]]} & /@xxls[[;; -2, -1]]];
        Show[g, g2]
    ] /; Length[timeseqlist] == Length[arraylist] == Length[colorlist]     
    
CtPathPlot[timeseq_?(VectorQ[#, NumericQ] &),discretepath_?(VectorQ[#, IntegerQ] &), linestyle___] :=
    Module[ {ylist, xx, yy, l},
        ylist = CtJumpFormat[discretepath];
        ylist[[-1, 1]] -= 1;
        xx = Riffle[ylist[[;; -2, 1]], ylist[[2 ;; -1, 1]]];
        yy = Riffle[ylist[[;; -2, 2]], ylist[[;; -2, 2]]];
        l = Line[Transpose[{timeseq[[xx]], yy}]];
        Graphics[{linestyle, l}]
    ] /; Length[timeseq] == Length[discretepath]     
    
CtPathPlot[timeseqlist_?(VectorQ[#, ListQ] &), pathlist_?(VectorQ[#, ListQ] &), stylelist_?VectorQ] :=
    Module[ {xxls = timeseqlist, len, ii, g},
        len = xxls[[;; -2, -1]];
        xxls[[2 ;;]] += Accumulate[len];
        g = Table[CtPathPlot[xxls[[ii]], pathlist[[ii]], stylelist[[ii]]], {ii,Length[xxls]}];
        Show[g]
    ] /; Length[timeseqlist] == Length[pathlist] == Length[stylelist]    
    
Options[CtSeqPlot] = Options[ListPlot];

CtSeqPlot[timeseq_?(VectorQ[#, NumericQ] &), yconst_?NumericQ, 
  opts : OptionsPattern[]] :=
    ListPlot[Thread[{timeseq, yconst}], opts]

CtSeqPlot[timeseqlist_?(VectorQ[#, ListQ] &), yconst_?NumericQ, 
  opts : OptionsPattern[]] :=
    Module[ {xxls = timeseqlist, len, g, ii},
        len = xxls[[;; -2, -1]];
        xxls[[2 ;;]] += Accumulate[len];
        g = Table[
          ListPlot[Thread[{xxls[[ii]], yconst}], opts], {ii, Length[xxls]}];
        Show[g]
    ]    
    
End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["ContinuousTimeHmm`*"];

EndPackage[]

