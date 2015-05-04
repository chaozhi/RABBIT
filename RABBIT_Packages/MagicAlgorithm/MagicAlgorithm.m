(* Mathematica Package *)

(* Created by the Wolfram Workbench 30-Nov-2014 *)

BeginPackage["MagicAlgorithm`",{"ContinuousTimeHmm`"}]
(* Exported symbols added here with SymbolName::usage *) 

Unprotect @@ Names["MagicAlgorithm`*"];
ClearAll @@ Names["MagicAlgorithm`*"];

origLogLiklihood::usage = "origLogLiklihood[startProb, tranProb, dataProb]  "

origPosteriorDecoding::usage = "origPosteriorDecoding[startProb, tranProb, dataProb]   "

origViterbiDecoding::usage = "origViterbiDecoding[startProb, tranProb, dataProb]  "

origPathSampling::usage = "origPathSampling[startProb, tranProb, dataProb, samplesize]  "

Begin["`Private`"]
(* Implementation of the package *)
      
origLogLiklihood[startProb_, tranProb_, dataProb_] :=
    Module[ {chr},
        Sum[CtLogLiklihood[startProb[[chr]], tranProb[[chr]], dataProb[[chr]]],{chr,Length[dataProb]}]
    ]
        
origPosteriorDecoding[startProb_, tranProb_, dataProb_] :=
    Module[ {logl, prob, chr},
        Table[
         {logl, prob} = CtPosteriorDecoding[startProb[[chr]], tranProb[[chr]], dataProb[[chr]]];
         prob = N[Round[prob, 10^(-5)]];
         {logl, prob}, {chr, Length[dataProb]}]
    ]        
    
origViterbiDecoding[startProb_, tranProb_, dataProb_] :=
    Module[ {logpathprob,path,chr},
        Table[
          {logpathprob, path} = CtViterbi[startProb[[chr]], tranProb[[chr]], dataProb[[chr]]];
          {logpathprob, CtJumpFormat[path]}, {chr, Length[dataProb]}]
    ]
    
origPathSampling[startProb_, tranProb_, dataProb_,size_:1] :=
    Module[ {logl,paths,chr},
        Table[
          {logl,paths} = CtPathSampling[startProb[[chr]], tranProb[[chr]], dataProb[[chr]],size];
          {logl,CtJumpFormat[#]&/@paths}, {chr, Length[dataProb]}]
    ]        
    

End[]

SetAttributes[#, {Protected,ReadProtected}]&/@ Names["MagicAlgorithm`*"];

EndPackage[]

