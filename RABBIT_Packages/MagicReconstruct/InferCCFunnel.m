(* Mathematica Package *)

BeginPackage["MagicReconstruct`InferCCFunnel`"]

inferCCFunnel::usage = "inferCCFunnel  "
(* Exported symbols added here with SymbolName::usage *)  

Begin["`Private`"] (* Begin Private Context *)

getlinepedinfor[gen_, lineid_, funnel_] :=
    Module[ {matingScheme, ped, nfounder = 8, saminfor},
        matingScheme =Join[Table["Pairing", {2}], Table["Sibling", {gen - 3}]];
        ped = simPedigree[nfounder, True, matingScheme];
        saminfor = Join[{{"ProgenyLine", "MemberID", "Funnelcode"}}, 
          Thread[{lineid, ped[[-1, 2]], {funnel}}]];
        mergePedigreeInfor[ped, saminfor]
    ]

getproppedinfor[gen_, lineid_, funnel_] :=
    Module[ {matingScheme, ped, nfounder = 8, saminfor, ls, ls2, posls, 
      funnells, prop,pos},
        ls = {{1, 3}, {1, 4}, {1, 5}, {1, 6}, {1, 7}, {1, 8}, {2, 3}, {2, 
           4}, {2, 5}, {2, 6}, {2, 7}, {2, 8}, {3, 4}, {3, 5}, {3, 6}, {3, 
           7}, {3, 8}, {4, 5}, {4, 6}, {4, 7}, {4, 8}, {5, 7}, {5, 8}, {6, 
           7}, {6, 8}};
        ls2 = {{1, 2, 3, 4}, {1, 2, 5, 6}, {1, 2, 7, 8}, {3, 4, 5, 6}, {3, 
           4, 7, 8}, {5, 6, 7, 8}};
        posls = Join[ls, ls2];
        funnells = Table[
          prop = funnel;
          prop[[pos]] = Flatten[Reverse[Partition[prop[[pos]], Length[pos]/2]]];
          prop, {pos, posls}];
        matingScheme = Join[Table["Pairing", {2}], Table["Sibling", {gen - 3}]];
        ped = simPedigree[nfounder, True, matingScheme];
        saminfor = Join[{{"ProgenyLine", "MemberID", "Funnelcode"}},Thread[{lineid, ped[[-1, 2]], funnells}]];
        mergePedigreeInfor[ped, saminfor]
    ]

estlinefunnel[magicsnp_List, model_String, epsF_?NonNegative, 
  eps_?NonNegative, ngeneration_Integer, initfunnel_List, maxit_: 20, 
  isprint_: True] :=
    Module[ {lineid = magicsnp[[-1, 1]], funnel = initfunnel, pedinfor, 
      logl, his, proppedinfor, propinfor, nprop, propmagicsnp, proplogl, 
      prop, pos,it},
        pedinfor = getlinepedinfor[ngeneration, lineid, funnel];
        logl = First[calOrigLogl[magicsnp, model, epsF, eps, pedinfor]];
        his = ConstantArray[0, {maxit + 1, 2}];
        his[[1]] = {funnel, logl};
        Do[
         proppedinfor = getproppedinfor[ngeneration, lineid, funnel];
         propinfor = splitPedigreeInfor[proppedinfor][[-1]];
         nprop = Length[propinfor] - 1;
         propmagicsnp = Join[Most[magicsnp], Table[magicsnp[[-1]], {nprop}]];
         proplogl = calOrigLogl[propmagicsnp, model, epsF, eps, proppedinfor];
         proplogl = Join[Rest[propinfor], List /@ proplogl, 2];
         pos = Ordering[proplogl[[All, -1]], -1][[1]];
         {prop, proplogl} = proplogl[[pos, {3, 4}]];
         If[ isprint,
             PrintTemporary["{iteration,lineid,funnel,logl,prop,proplogl}=", {it, lineid, 
                funnel, logl, prop, proplogl}];
         ];
         If[ proplogl >= logl,
             {funnel, logl} = {prop, proplogl};
         ];
         his[[it + 1]] = {funnel, logl};
         If[ his[[it + 1, 2]] == his[[it, 2]],
             Break[]
         ], {it, maxit + 1}];
        his = DeleteCases[his, {0, 0}];
        {lineid, his[[-1, 1]], his}
    ]

Options[inferCCFunnel] = {
	founderAllelicError -> 0.005,
	offspringAllelicError -> 0.005, 
	outputFileID -> "",
	isPrintTimeElapsed -> True
	}
	
inferCCFunnel[magicSNP_List, model_String, ngeneration_List, initfunnel_List: Automatic, opts : OptionsPattern[]] :=
    Module[ {starttime, resls, nfounder, lineid, maxit = 20, count, bool, 
      linemagicsnp, epsF, eps, outputid, outputfile, isprint, ngen, 
      funnel, res, ii},
        {epsF, eps, outputid, isprint} = 
         OptionValue@{founderAllelicError, offspringAllelicError, 
           outputFileID, isPrintTimeElapsed};
        If[ outputid =!= "",
            outputfile = outputid <> "_inferCCFunnel.txt",
            outputfile = "inferCCFunnel.txt"
        ];
        nfounder = magicSNP[[1, 2]];
        lineid = magicSNP[[5 + nfounder ;;, 1]];
        If[ ngeneration[[2 ;;, 1]] != lineid,
            Print["Inconsistent line ID between magicSNP and ngeneration!"];
            Abort[]
        ];
        If[ initfunnel[[2 ;;, 1]] != lineid,
            Print["Inconsistent line ID between magicSNP and initfunnel!"];
            Abort[]
        ];
        bool = (Union[Flatten[#]] =!= Range[nfounder]) & /@ 
          initfunnel[[2 ;;, 2]];
        If[ MemberQ[bool, True],
            Print[Join[{{"LineID", "WrongInitialFunnel"}}, 
               Pick[Rest[initfunnel], bool]] // MatrixForm];
            Abort[]
        ];
        If[ isprint,
            starttime = SessionTime[];
            Print["magicReconstruct. Start date = ", DateString[], 
             "\toutputfile = ", outputfile];
        ];
        resls = Table[0, {Length[initfunnel] - 1}];
        count = 0;
        SetSharedVariable[resls, count];
        Monitor[ParallelDo[
          funnel = initfunnel[[ii + 1, 2]];
          ngen = ngeneration[[ii + 1, 2]];
          linemagicsnp = magicSNP[[Join[Range[nfounder + 4], nfounder + 4 + {ii}]]];
          resls[[ii]] = 
           res = estlinefunnel[linemagicsnp, model, epsF, eps, ngen, funnel,maxit, False];
          If[ isprint,
              PrintTemporary["TimeElapsed = ",Round[SessionTime[] - starttime, 0.1], 
                " seconds. {No., Lineid,Funnel,History}=", Join[{ii}, res]];
          ];
          count++, {ii, Length[initfunnel] - 1}], 
         ProgressIndicator[count, {0, Length[initfunnel] - 1}]];
        resls = Join[{{"LineID", "Funnel", "History"}}, resls];
        If[ isprint,
            Print["Done! Finished date =", DateString[], ". \tTime elapsed = ",
               Round[SessionTime[] - starttime, 0.1], " Seconds."];
        ];
        Put[Sequence @@ resls, outputfile]
    ] 

End[] (* End Private Context *)

EndPackage[]