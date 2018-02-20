(* Mathematica Package *)

(* Created by the Wolfram Workbench Jul 10, 2015 *)

(*removeDownValues from MagicDefinition` is used in this package*)

BeginPackage["PedigreeOrigin`",{"MagicDefinition`"}]

(*pedigree is arranged so that all the founders on the top, and parents are before children*)

pedIdentitySummary::usage = "pedIdentitySummary[pedigree,founderFGL,isAutosome,pairlist] "

pedIdentityConciseSummary::usage = "pedIdentityConciseSummary  "

pedIdentityPrior::usage = "pedIdentityPrior[pedigree,founderFGL,isAutosome,indlist]  "

pedIdentityConcisePrior::usage = "pedIdentityConcisePrior  "

pedAncestrySummary::usage = "pedAncestrySummary[pedigree,founderFGL,isAutosome,pairlist] "

pedAncestryConciseSummary::usage = "pedAncestryConciseSummary  "

pedAncestryPrior::usage = "pedAncestryPrior[pedigree,founderFGL,isAutosome,indlist]  "

pedAncestryConcisePrior::usage = "pedAncestryConcisePrior  "

pedAncestryCompatibleQ::usage = "pedAncestryCompatibleQ  "

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package *)

$RecursionLimit = Infinity;

ClearAll[x,x1,x2,x3,x4,x5,x6,i,j,k]

getinitial2[{{x1_,x2_},{x3_,x4_}}] = Boole[Unequal@@#&/@Tuples[{{x1,x2},{x3,x4}}]]

getinitial3[{{x1_,x2_},{x3_,x4_},{x5_,x6_}}] = Boole[Unequal@@#&/@Tuples[{{x1,x2},{x3,x4},{x5,x6}}]]
   
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*phi12=two-gene non-ibd probability; phi12[a] ={mm,mp,pm,pp}*)
(*Boole[Unequal @@ # & /@ Tuples[$founderfgl[[{a,b}]]]]*)
phi12[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    phi12[a, b] = N[getinitial2[$founderfgl[[{a,b}]]]]    
phi12[a_, b_] :=
    phi12[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a >= b,         
         mm = Boole[a =!= b] Mean[phi12[$mothers[[a]], b][[{1, 3}]]];
         mp = Mean[phi12[$mothers[[a]], b][[{2, 4}]]];
         If[ $isautosome,
             pm = Mean[phi12[$fathers[[a]], b][[{1, 3}]]];
             pp = Boole[a =!= b] Mean[phi12[$fathers[[a]], b][[{2, 4}]]],
             If[ $genders[[a]]==1,
                 (*female a, X chromsome*)
                 pm = phi12[$fathers[[a]], b][[1]];
                 pp = Boole[a =!= b] phi12[$fathers[[a]], b][[2]],
                 (*male a, Y chromsome*)
                 pm = phi12[$fathers[[a]], b][[3]];
                 pp = Boole[a =!= b] phi12[$fathers[[a]], b][[4]]
             ]
         ];
         {mm, mp, pm, pp},
         a < b,
         phi12[b, a][[{1, 3, 2, 4}]]
         ]
    ];
        
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"},{"m","p"}}]*)
(*phi123=three-gene non-ibd probability; phi123[a,b,c]={mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}*)
(*Boole[Unequal @@ # & /@ Tuples[$founderfgl[[{a,b,c}]]]]*)
phi123[a_?(#<=$nfounder&), b_?(#<=$nfounder&),c_?(#<=$nfounder&)] :=
    phi123[a, b] = N[getinitial3[$founderfgl[[{a,b,c}]]]]
phi123[a_, b_, c_] :=
    phi123[a, b, c] =
    Module[ {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp,indicator},
        Which[
         a == b == c,
         Table[N[0], {8}],
         a >= b >= c,          
         {mmm, mmp, mpm, mpp} = (Mean[phi123[$mothers[[a]], b, c][[{0, 4} + #]]] & /@Range[4]);
         If[ $isautosome,
             {pmm, pmp, ppm, ppp} = (Mean[phi123[$fathers[[a]], b, c][[{0, 4} + #]]] & /@Range[4]),
             If[ $genders[[a]]==1,
                 {pmm, pmp, ppm, ppp} = phi123[$fathers[[a]], b, c][[;; 4]],
                 {pmm, pmp, ppm, ppp} = phi123[$fathers[[a]], b, c][[5;;]]
             ]
         ];
         indicator = {Boole[a=!=b],Boole[a=!=b],1,1} {Boole[b=!=c],1, 1, Boole[b=!=c]};
         {mmm, mmp, mpm, mpp} = indicator {mmm, mmp, mpm, mpp};
         indicator = {1,1,Boole[a=!=b],Boole[a=!=b]} {Boole[b=!=c],1, 1, Boole[b=!=c]};
         {pmm, pmp, ppm, ppp} = indicator {pmm, pmp, ppm, ppp};
         {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp},
         (*{{a,b,c},{a,c,b},{b,a,c},{b,c,a},{c,a,b},{c,b,a}}*)
         a >= c >= b,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
         {mmm,mpm,mmp,mpp,pmm,ppm,pmp,ppp}*)
         phi123[a, c, b][[{1, 3, 2, 4, 5, 7, 6, 8}]],
         b >= a >= c,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
         {mmm,mmp,pmm,pmp,mpm,mpp,ppm,ppp}*)
         phi123[b, a, c][[{1, 2, 5, 6, 3, 4, 7, 8}]],
         b >= c >= a,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
         {mmm,mpm,pmm,ppm,mmp,mpp,pmp,ppp}*)
         phi123[b, c, a][[{1, 3, 5, 7, 2, 4, 6, 8}]],
         c >= a >= b,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
         {mmm,pmm,mmp,pmp,mpm,ppm,mpp,ppp}*)
         phi123[c, a, b][[{1, 5, 2, 6, 3, 7, 4, 8}]],
         c >= b >= a,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
         {mmm,pmm,mpm,ppm,mmp,pmp,mpp,ppp}*)
         phi123[c, b, a][[{1, 5, 3, 7, 2, 6, 4, 8}]],
         True,
         Print["Missing scenario: {a,b,c} = ", {a, b, c}];
         Abort[];
         ]
    ];
                    
(*R12[a]= {Rm,Rp} map expansion for maternally and patenrally derived chromosomes*)
R12[a_?(#<=$nfounder&)] :=
    R12[a] = N[{0,0}]
R12[a_] :=
    R12[a] = 
    If[ $isautosome,
        {Mean[R12[$mothers[[a]]]] + phi12[$mothers[[a]], $mothers[[a]]][[2]], Mean[R12[$fathers[[a]]]] + phi12[$fathers[[a]], $fathers[[a]]][[2]]},
        {Mean[R12[$mothers[[a]]]] + phi12[$mothers[[a]], $mothers[[a]]][[2]], R12[$fathers[[a]]][[$genders[[a]]]]}
    ];
        
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*junc1122=two-locus junction density; junc1122[a,b]={mm,mp,pm,pp}*)
junc1122[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    junc1122[a, b] = N[{0,0,0,0}]
junc1122[a_, b_] :=
    junc1122[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a == b,
         {mm, pp} = R12[a];
         mp = pm = Mean[junc1122[a, $mothers[[a]]][[{3, 4}]]];
         {mm, mp, pm, pp},
         a > b,
         {mm, mp} = (Mean[junc1122[$mothers[[a]], b][[{0, 2} + #]]] & /@Range[2]);
         If[ $isautosome,
             {pm, pp} = (Mean[junc1122[$fathers[[a]], b][[{0, 2} + #]]] & /@Range[2]),
             If[ $genders[[a]]==1,
                 {pm, pp} = junc1122[$fathers[[a]], b][[;;2]],
                 {pm, pp} = junc1122[$fathers[[a]], b][[3;;]]
             ]
         ];
         {mm, mp, pm, pp},
         a < b,
         junc1122[b, a][[{1, 3, 2, 4}]]
         ]
    ];
    
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*junc1232=two-locus junction density; junc1232[a,b]={mm,mp,pm,pp}*)
(*phi123=three-gene non-ibd probability; phi123[a,b,c]={mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}*)
junc1232[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    junc1232[a, b] = N[{0,0,0,0}]
junc1213[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    junc1213[a, b] = N[{0,0,0,0}]    
junc1232[a_, b_] :=
    junc1232[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a >= b,
         mm = Boole[a=!=b] (Mean[junc1232[$mothers[[a]], b][[{1,3}]]] + phi123[$mothers[[a]], $mothers[[a]], b][[3]]);
         mp = (Mean[junc1232[$mothers[[a]], b][[{2,4}]]] + phi123[$mothers[[a]], $mothers[[a]], b][[4]]);
         If[ $isautosome,
             pm = (Mean[junc1232[$fathers[[a]], b][[{1, 3}]]] + phi123[$fathers[[a]], $fathers[[a]], b][[3]]);
             pp = Boole[a=!=b] (Mean[junc1232[$fathers[[a]], b][[{2, 4}]]] + phi123[$fathers[[a]], $fathers[[a]], b][[4]]),
             If[ $genders[[a]]==1,
                 {pm, pp} = {1,  Boole[a=!=b]} junc1232[$fathers[[a]], b][[;;2]],
                 {pm, pp} = {1,  Boole[a=!=b]} junc1232[$fathers[[a]], b][[3;;]]
             ]
         ];
         {mm, mp, pm, pp},
         a < b,
         junc1213[b, a][[{1, 3, 2, 4}]]
         ]
    ];
junc1213[a_, b_] :=
    junc1213[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a >= b,
         mm = Boole[a=!=b] Mean[junc1213[$mothers[[a]], b][[{1,3}]]];
         mp = Mean[junc1213[$mothers[[a]], b][[{2,4}]]];
         If[ $isautosome,
             pm = Mean[junc1213[$fathers[[a]], b][[{1, 3}]]];
             pp = Boole[a=!=b] Mean[junc1213[$fathers[[a]], b][[{2, 4}]]],
             If[ $genders[[a]]==1,
                 {pm, pp} = {1, Boole[a=!=b]} junc1213[$fathers[[a]], b][[;;2]],
                 {pm, pp} = {1, Boole[a=!=b]} junc1213[$fathers[[a]], b][[3;;]]
             ]
         ];
         {mm, mp, pm, pp},
         a < b,
         junc1232[b, a][[{1, 3, 2, 4}]]
         ]
    ];
    
(*Internal package variables: $isautosome,$mothers,$fathers,$genders,$founderfgl,$nfounder*)    
setstage[pedigree_,founderFGL_,isAutosome_,pairlist_] :=
    Module[ {rule,pairlist2,wronggender},
        $isautosome = isAutosome;
        $founderfgl = founderFGL;
        $nfounder = Length[founderFGL];
        $fgllabset = Union[Flatten[founderFGL]];
        $nfgllab = Length[$fgllabset];
        If[ $nfounder!= Count[pedigree[[2 ;;, 4]], {0, 0}],
            Print["setstage: inconsistent $nfounder between pedigree and founderFGL!"];
            Abort[]
        ];
        rule = Thread[pedigree[[2 ;;, 2]] -> Range[Length[pedigree[[2 ;;, 2]]]]];
        {$mothers, $fathers} = Transpose[pedigree[[2 ;;, 4]] /. rule];
        $genders = pedigree[[2 ;;, 3]];
        If[ Complement[Union[$genders], {0, 1, 2}] =!= {},
            Print["The genders must be Female=1/Male=2/hermaphrodite=0!"];
            Abort[]
        ];
        wronggender = Select[Transpose[{pedigree[[2 ;;, 2]], $genders, $genders[[$mothers]]}], !MemberQ[{List, 1, 0}, #[[3]]] &];
        If[ wronggender =!= {},
            Print["The geneders of Mothers must be Female=1/hermaphrodite=0!"];
            Print[Join[{{pedigree[[1]],pedigree[[2]],"motherer gender"}},wronggender]//MatrixForm];
            Abort[]
        ];
        wronggender = Select[Transpose[{pedigree[[2 ;;, 2]], $genders, $genders[[$fathers]]}], !MemberQ[{List, 2, 0}, #[[3]]] &];
        If[ wronggender =!= {},
            Print["The geneders of Mothers must be Male=2/hermaphrodite=0!"];
            Print[Join[{{pedigree[[1]],pedigree[[2]],"motherer gender"}},wronggender]//MatrixForm];
            Abort[]
        ];
        pairlist2 = pairlist /. rule;
        pairlist2
    ]  
        
clearIdentityMemory[] :=
    Module[ {},
        removeDownValues[phi12[_?IntegerQ, _?IntegerQ]];
        removeDownValues[phi123[_?IntegerQ, _?IntegerQ,_?IntegerQ]];
        removeDownValues[R12[_?IntegerQ]];
        removeDownValues[junc1122[_?IntegerQ, _?IntegerQ]];
        removeDownValues[junc1232[_?IntegerQ, _?IntegerQ]];
        removeDownValues[junc1213[_?IntegerQ, _?IntegerQ]];
    ]

clearIdentityConciseMemory[] :=
    Module[ {},
        removeDownValues[phi12[_?IntegerQ, _?IntegerQ]];
        removeDownValues[R12[_?IntegerQ]];
        removeDownValues[junc1122[_?IntegerQ, _?IntegerQ]];
    ]
    

pedIdentitySummary[pedigree_,founderFGL_,isAutosome_,pairlist_] :=
    Module[ {pairlist2,ibd12,Rab,j1112,j1121,j1122, j1211, j1213, j1222, j1232,rho,ls,res},
        If[ Count[pedigree[[2 ;;, 4]], {0, 0}]=!=Length[founderFGL],
            Print["Inconsistent number of founders between pedigree and founderFGL!"];
            Abort[]
        ];
        clearIdentityMemory[];        
        (*setstage calculate internal package variables: $isautosome,$mothers,$fathers,$genders*)
        pairlist2 = setstage[pedigree,founderFGL,isAutosome,pairlist];
        {ibd12,Rab, j1122,j1232,j1213} = {phi12 @@ # & /@ pairlist2, Map[R12, pairlist2, {2}],
            junc1122@@ # & /@ pairlist2,junc1232@@ # & /@ pairlist2,junc1213@@ # & /@ pairlist2};
        rho = Transpose[Join[Rab[[All, 1, 1]] + # & /@ Transpose[Rab[[All, 2]]],Rab[[All, 1, 2]] + # & /@ Transpose[Rab[[All, 2]]]]] - j1122;
        clearIdentityMemory[];
        ls = Transpose[{Rab[[All, 1, 1]], Rab[[All, 1, 1]], Rab[[All, 1, 2]], Rab[[All, 1, 2]]}];
        j1121 = j1222 = (ls - j1122 - j1232)/2;
        ls = Transpose[{Rab[[All, 2, 1]], Rab[[All, 2, 2]], Rab[[All, 2, 1]], Rab[[All, 2, 2]]}];
        j1112 = j1211 = (ls - j1122 - j1213)/2;
        res = Transpose[{pairlist, ibd12,Rab[[All,1]],Rab[[All,2]],rho,j1112,j1121,j1122, j1211, j1213, j1222, j1232}];
        res = Join[{{"{a,b}", "phi(12)_{ab}^{mm,mp,pm,pp}","R_a^{m,p}","R_b^{m,p}","rho_{ab}^{mm,mp,pm,pp}",
             "J1112_{ab}^{mm,mp,pm,pp}","J1121_{ab}^{mm,mp,pm,pp}","J1122_{ab}^{mm,mp,pm,pp}",
             "J1211_{ab}^{mm,mp,pm,pp}","J1213_{ab}^{mm,mp,pm,pp}","J1222_{ab}^{mm,mp,pm,pp}","J1232_{ab}^{mm,mp,pm,pp}"}}, res];
        res
    ]    

pedIdentityConciseSummary[pedigree_,founderFGL_,isAutosome_,pairlist_] :=
    Module[ {pairlist2,ibd12,Rab,j1122,rho,res},
        If[ Count[pedigree[[2 ;;, 4]], {0, 0}]=!=Length[founderFGL],
            Print["Inconsistent number of founders between pedigree and founderFGL!"];
            Abort[]
        ];
        clearIdentityConciseMemory[];        
        (*setstage calculate internal package variables: $isautosome,$mothers,$fathers,$genders*)
        pairlist2 = setstage[pedigree,founderFGL,isAutosome,pairlist];
        {ibd12, Rab, j1122} = {phi12 @@ # & /@ pairlist2, Map[R12, pairlist2, {2}], junc1122 @@ # & /@ pairlist2};
        rho = Transpose[Join[Rab[[All, 1, 1]] + # & /@ Transpose[Rab[[All, 2]]],Rab[[All, 1, 2]] + # & /@ Transpose[Rab[[All, 2]]]]] - j1122;
        clearIdentityMemory[];
        res = Transpose[{pairlist, ibd12,Rab[[All,1]],Rab[[All,2]],rho}];
        res = Join[{{"{a,b}", "phi(12)_{ab}^{mm,mp,pm,pp}","R_a^{m,p}","R_b^{m,p}",
             "rho_{ab}^{mm,mp,pm,pp}"}}, res];
        res
    ]  
            

pedIdentityPrior[pedigree_, founderFGL_, isAutosome_, indlist_] :=
    Module[ {ped,pairlist, res, lab, pos,pos2,str},
        pos = Flatten[Position[pedigree[[2 ;;, 2]], #, {1}, 1, Heads -> False] & /@indlist];
        If[ pedigree[[pos + 1, 2]] == indlist,
            ped = pedigree[[pos + 1, ;;4]],
            Print["pedIdentityPrior: Wrong indlist!"];
            Abort[];
        ];
        pairlist = Transpose[{indlist, indlist}];
        res = pedIdentitySummary[pedigree, founderFGL, isAutosome, pairlist];
        lab = {"phi(12)_{aa}^mp","rho_{aa}^{mp}","J1112_{aa}^mp","J1121_{aa}^mp","J1122_{aa}^mp", "J1211_{aa}^mp", "J1213_{aa}^mp", "J1222_{aa}^mp", "J1232_{aa}^mp"};
        pos = Flatten[Table[Position[res[[1]], _?(StringMatchQ[#, str ~~ ___] &), {1},Heads -> False], {str, StringTake[#, 5] & /@ lab}]];
        If[ Length[pos] != Length[lab],
            Print["pedIdentityPrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        pos2 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "R_a" ~~ ___] &), {1},Heads -> False]];
        If[ Length[pos2] != 1,
            Print["pedIdentityPrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        res = Join[{Join[pedigree[[1, ;; 4]], {"R_a^m", "R_a^p"}, lab]}, 
           Join[ped, res[[2 ;;, pos2[[1]]]], res[[2 ;;, pos, 2]], 2]];
        res[[All, 5 ;; 7]] = res[[All, {7, 5, 6}]];        
        (*change non-IBD probability into IBD probability*)
        pos = Dimensions[ped][[2]]+1; (*the column of the IBD probability*)
        res[[1, pos]] = StringReplace[res[[1, pos]], "12" -> "11"];
        res[[2 ;;, pos]] = 1 - res[[2 ;;, pos]];
        res
    ]

pedIdentityConcisePrior[pedigree_, founderFGL_, isAutosome_, indlist_] :=
    Module[ {ped,pairlist, res, lab, pos,pos2,str},
        pos = Flatten[Position[pedigree[[2 ;;, 2]], #, {1}, 1, Heads -> False] & /@indlist];
        If[ pedigree[[pos + 1, 2]] == indlist,
            ped = pedigree[[pos + 1, ;; 4]],
            Print["pedIdentityPrior: Wrong indlist!"];
            Abort[];
        ];
        pairlist = Transpose[{indlist, indlist}];
        res = pedIdentityConciseSummary[pedigree, founderFGL, isAutosome, pairlist];
        lab = {"phi(12)_{aa}^mp", "rho_{aa}^mp"};
        pos = Flatten[Table[Position[res[[1]], _?(StringMatchQ[#, str ~~ ___] &), {1}, 
             Heads -> False], {str, StringTake[#, 4] & /@ lab}]];
        If[ Length[pos] != Length[lab],
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        pos2 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "R_a" ~~ ___] &), {1},Heads -> False]];
        If[ Length[pos2] != 1,
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        res = Join[{Join[pedigree[[1, ;; 4]], {"R_a^m", "R_a^p"}, lab]}, 
           Join[ped, res[[2 ;;, pos2[[1]]]], res[[2 ;;, pos, 2]], 2]];
        res[[All, 5 ;; 7]] = res[[All, {7, 5, 6}]];
        (*change non-IBD probability into IBD probability*)
        pos = Dimensions[ped][[2]]+1; (*the column of the IBD probability*)
        res[[1, pos]] = StringReplace[res[[1, pos]], "12" -> "11"];
        res[[2 ;;, pos]] = 1 - res[[2 ;;, pos]];
        res
    ]        

fi[a_?(#<=$nfounder&)] :=
    fi[a] = Map[SparseArray, Replace[$fgllabset, {# -> 1, _ -> 0}, {1}] & /@ $founderfgl[[a]]]
fi[a_] :=
    fi[a] =
    Module[ {m, p},
        m = Mean[fi[$mothers[[a]]]];
        p = If[ $isautosome,
                Mean[fi[$fathers[[a]]]],
                fi[$fathers[[a]]][[$genders[[a]]]]
            ];
        {m, p}
    ];

(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*phiij=two-gene ancestry probability; each element ={mm,mp,pm,pp}*)
(*2 refers to two copies of genes at each loci of an diploid individual *)
phiij[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    phiij[a, b] = Map[SparseArray, Flatten[Table[Outer[Times, fi[a][[i]], fi[b][[j]]], {i, 2}, {j, 2}], 1]]
phiij[a_, b_] :=
    phiij[a, b] =
    Module[ {mm, mp, pm, pp, i},
        Which[
         a == b,
         {mm, pp} = Table[DiagonalMatrix[fi[a][[i]]], {i, 2}];
         pm = Mean[phiij[a, $mothers[[a]]][[{3, 4}]]];
         mp = Transpose[pm];
         {mm, mp, pm, pp},
         a > b,
         mm = Mean[phiij[$mothers[[a]], b][[{1, 3}]]];
         mp = Mean[phiij[$mothers[[a]], b][[{2, 4}]]];
         If[ $isautosome,
             pm = Mean[phiij[$fathers[[a]], b][[{1, 3}]]];
             pp = Mean[phiij[$fathers[[a]], b][[{2, 4}]]],
             If[ $genders[[a]] == 1,
                 pm = phiij[$fathers[[a]], b][[1]];
                 pp = phiij[$fathers[[a]], b][[2]],
                 pm = phiij[$fathers[[a]], b][[3]];
                 pp = phiij[$fathers[[a]], b][[4]];
             ];
         ];
         {mm, mp, pm, pp},
         a < b,
         Transpose[#]& /@ phiij[b, a][[{1, 3, 2, 4}]]
         ]
    ]
    
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"},{"m","p"}}]*)
 (*phiijk=three-gene non-ibd probability; phiijk[a,b,c]={mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}*)
phiijk[a_?(#<=$nfounder&), b_?(#<=$nfounder&),  c_?(#<=$nfounder&)] :=
    phiijk[a, b,c] = Map[SparseArray, Flatten[Table[Outer[Times, fi[a][[i]], fi[b][[j]], fi[c][[k]]], {i, 2}, {j, 2}, {k, 2}], 2]]    
phiijk[a_, b_, c_] :=
    phiijk[a, b, c] =
    Module[ {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp, i, j, k, ch,nfgl = Length[First[fi[1]]]},
        Which[
         a == b == c,
         {mmm, ppp} = Map[SparseArray, Table[Boole[i == j == k] fi[a][[ch, i]], {ch, 2}, {i,nfgl}, {j, nfgl}, {k, nfgl}]];
         {mmp, ppm} = Map[SparseArray, Table[Boole[i == j] phiij[a, a][[ch,i]], {ch, 2, 3}, {i, nfgl}, {j, nfgl}]];
         {mpm, pmp} = Transpose[#, {1, 3, 2}] & /@ {mmp, ppm};
         {mpp, pmm} = Transpose[#, {2, 3, 1}] & /@ {ppm, mmp};
         {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp},
         a == b > c,
         {mmm, mmp, ppm, ppp} = Map[SparseArray,Table[Boole[i == j] phiij[a, c][[ch, i]], {ch, 4}, {i,nfgl}, {j, nfgl}]];
         (*pmm=5,ppm=7*)
         pmm = Mean[phiijk[a, $mothers[[a]], c][[{5, 7}]]];
         mpm = Transpose[pmm, {2, 1, 3}];
         (*pmp=6,ppp=8*)
         pmp = Mean[phiijk[a, $mothers[[a]], c][[{6, 8}]]];
         mpp = Transpose[pmp, {2, 1, 3}];
         {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp},
         a > b == c,
         {mmm, mpp, pmm, ppp} = Map[SparseArray,Table[phiij[a, b][[ch, i]] IdentityMatrix[nfgl], {ch, 4}, {i,nfgl}]];
         {mmp, mpm} = (Mean[phiijk[$mothers[[a]], b, c][[{0, 4} + #]]] & /@ {2,3});
         If[ $isautosome,
             {pmp, ppm} = (Mean[phiijk[$fathers[[a]], b, c][[{0, 4} + #]]] & /@ {2,3}),
             If[ $genders[[a]] == 1,
                 {pmp, ppm} = phiijk[$fathers[[a]], b, c][[{2, 3}]],
                 {pmp, ppm} = phiijk[$fathers[[a]], b, c][[{6, 7}]]
             ];
         ];
         {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp},
         a > b > c,
         {mmm, mmp, mpm, mpp} = (Mean[phiijk[$mothers[[a]], b, c][[{0, 4} + #]]] & /@Range[4]);
         If[ $isautosome,
             {pmm, pmp, ppm,ppp} = (Mean[phiijk[$fathers[[a]], b, c][[{0, 4} + #]]] & /@Range[4]),
             If[ $genders[[a]] == 1,
                 {pmm, pmp, ppm, ppp} = phiijk[$fathers[[a]], b, c][[;; 4]],
                 {pmm, pmp, ppm, ppp} = phiijk[$fathers[[a]], b, c][[5 ;;]]
             ];
         ];
         {mmm, mmp, mpm, mpp, pmm, pmp, ppm, ppp},
         (*{{a,b,c},{a,c,b},{b,a,c},{b,c,a},{c,a,b},{c,b,a}}*)
         a > c > b || a == c > b || a > c == b,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp} 
         {mmm,mpm,mmp,mpp,pmm,ppm,pmp,ppp}*)
         Transpose[#, {1, 3, 2}] & /@ phiijk[a, c, b][[{1, 3, 2, 4, 5, 7, 6, 8}]],
         b > a > c || b == a > c || b > a == c,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp} 
         {mmm,mmp,pmm,pmp,mpm,mpp,ppm,ppp}*)
         Transpose[#, {2, 1, 3}] & /@ phiijk[b, a, c][[{1, 2, 5, 6, 3, 4, 7, 8}]],
         b > c > a || b == c > a || b > c == a,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp} 
         {mmm,mpm,pmm,ppm,mmp,mpp,pmp,ppp}*)
         Transpose[#, {2, 3, 1}] & /@ phiijk[b, c, a][[{1, 3, 5, 7, 2, 4, 6, 8}]],
         c > a > b || c == a > b || c > a == b,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp} 
         {mmm,pmm,mmp,pmp,mpm,ppm,mpp,ppp}*)
         Transpose[#, {3, 1, 2}] & /@ phiijk[c, a, b][[{1, 5, 2, 6, 3, 7, 4, 8}]],
         c > b > a || c == b > a || c > b == a,
         (*{mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}
          {mmm,pmm,mpm,ppm,mmp,pmp,mpp,ppp}*)
         Transpose[#, {3, 2, 1}] & /@ phiijk[c, b, a][[{1, 5, 3, 7, 2, 6, 4, 8}]],
         True,
         Print["Missing scenario: {a,b,c} = ", {a, b, c}];
         Abort[];
         ]
    ];
         
(*Rij= {Rm,Rp} map expansion for maternally and patenrally derived chromosomes*)
Rij[a_?(#<=$nfounder&)] :=
    Rij[a] = Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab}]], {2}]
Rij[a_] :=
    Rij[a] = Module[ {nondiag2,nfgl = Length[First[fi[1]]]},
                 nondiag2 = 1 - SparseArray[{x_, x_} -> 1, {nfgl,nfgl}];
                 If[ $isautosome,
                     {Mean[Rij[$mothers[[a]]]] +Mean[phiij[$mothers[[a]], $mothers[[a]]][[{2, 3}]]] nondiag2, 
                      Mean[Rij[$fathers[[a]]]] + Mean[phiij[$fathers[[a]], $fathers[[a]]][[{2, 3}]]] nondiag2},
                     {Mean[Rij[$mothers[[a]]]] + Mean[phiij[$mothers[[a]], $mothers[[a]]][[{2, 3}]]] nondiag2, 
                      Rij[$fathers[[a]]][[$genders[[a]]]]}
                 ]
             ];

(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*junciijj=two-locus junction density;junciijj[a,b]={mm,mp,pm,pp}*)
junciijj[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    junciijj[a, b] =  Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab}]], {4}]
junciijj[a_, b_] :=
    junciijj[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a == b,
         {mm, pp} = Rij[a];
         mp = pm = Mean[junciijj[a, $mothers[[a]]][[{3, 4}]]];
         {mm, mp, pm, pp},
         a > b,
         {mm, mp} = (Mean[junciijj[$mothers[[a]], b][[{0, 2} + #]]] & /@Range[2]);
         If[ $isautosome,
             {pm, pp} = (Mean[junciijj[$fathers[[a]], b][[{0, 2} + #]]] & /@Range[2]),
             If[ $genders[[a]] == 1,
                 {pm, pp} = junciijj[$fathers[[a]], b][[;; 2]],
                 {pm, pp} = junciijj[$fathers[[a]], b][[3 ;;]]
             ];
         ];
         {mm, mp, pm, pp},
         a < b,
         junciijj[b, a][[{1, 3, 2, 4}]]
         ]
    ];
        
(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*juncijkj=two-locus junction density;juncijkj[a,b]={mm,mp,pm,pp}*)
(*phiijk=three-gene non-ibd probability;phiijk[a,b,c]={mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}*)
juncijkj[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    juncijkj[a, b] = Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab,$nfgllab}]], {4}]
juncijkj[a_, b_] :=
    juncijkj[a, b] =
    Module[ {mm, mp, pm, pp,nfgl,nondiag3,temp},
        nfgl = Length[First[fi[1]]];
        nondiag3 = 1 - SparseArray[{x_, x_, y_} | {x_, y_, x_} | {y_, x_, x_} -> 1, {nfgl,nfgl,nfgl}];
        Which[
         a >= b,
         temp = Map[Transpose,phiijk[$mothers[[a]], $mothers[[a]], b],{2}];
         temp = nondiag3 SparseArray[Mean[temp[[{3,5}+#]]]]&/@{0,1};
         {mm, mp} = (Mean[juncijkj[$mothers[[a]], b][[{1,3} + #]]] & /@ {0,1})+temp;
         If[ $isautosome,
             temp = Map[Transpose,phiijk[$fathers[[a]], $fathers[[a]], b],{2}];
             temp = nondiag3 SparseArray[Mean[temp[[{3,5}+#]]]]&/@{0,1};
             {pm, pp} = (Mean[juncijkj[$fathers[[a]], b][[{1, 3} + #]]]& /@ {0,1})+temp,
             If[ $genders[[a]] == 1,
                 {pm, pp} = juncijkj[$fathers[[a]], b][[;; 2]],
                 {pm, pp} = juncijkj[$fathers[[a]], b][[3 ;;]]
             ];
         ];
         {mm, pp} = Boole[a=!=b] {mm, pp};
         {mm, mp, pm, pp},
         a < b,
         Transpose[#, {2, 1, 3}] & /@ juncijik[b, a][[{1, 3, 2, 4}]]
         ]
    ];

juncijik[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    juncijik[a, b] = Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab,$nfgllab}]], {4}]    
juncijik[a_, b_] :=
    juncijik[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[         
         a >= b,
         {mm, mp} = (Mean[juncijik[$mothers[[a]], b][[{0, 2} + #]]] & /@Range[2]);
         If[ $isautosome,
             {pm, pp} = (Mean[juncijik[$fathers[[a]], b][[{0, 2} + #]]] & /@Range[2]),
             If[ $genders[[a]] == 1,
                 {pm, pp} = juncijik[$fathers[[a]], b][[;; 2]],
                 {pm, pp} = juncijik[$fathers[[a]], b][[3 ;;]]
             ];
         ];
         {mm, pp} = Boole[a=!=b] {mm, pp};
         {mm, mp, pm, pp},
         a < b,
         Transpose[#, {2, 1, 3}] & /@ juncijkj[b, a][[{1, 3, 2, 4}]]
         ]
    ];


(*StringJoin[#]&/@Tuples[{{"m","p"},{"m","p"}}]*)
(*juncijjj=two-locus junction density;juncijjj[a,b]={mm,mp,pm,pp}*)
(*phiijk=three-gene non-ibd probability;phiijk[a,b,c]={mmm,mmp,mpm,mpp,pmm,pmp,ppm,ppp}*)
juncijjj[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    juncijjj[a, b] =  Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab}]], {4}]    
juncijjj[a_, b_] :=
    juncijjj[a, b] =
    Module[ {mm, mp, pm, pp,nfgl,nondiag2,temp},
        nfgl = Length[First[fi[1]]];
        nondiag2 = 1 - SparseArray[{x_, x_} -> 1, {nfgl,nfgl}];
        Which[
         a >= b,     
         temp = SparseArray[#]& /@ Map[Diagonal, phiijk[$mothers[[a]], $mothers[[a]], b], {2}];
         temp = {nondiag2 Mean[temp[[{3,5}]]],nondiag2 Mean[temp[[{4,6}]]]};
         {mm, mp} = (Mean[juncijjj[$mothers[[a]], b][[{1,3} + #]]] & /@ {0,1})+temp;
         If[ $isautosome,
             temp = SparseArray[#]& /@ Map[Diagonal, phiijk[$fathers[[a]], $fathers[[a]], b], {2}];
             temp = {nondiag2 Mean[temp[[{3,5}]]],nondiag2 Mean[temp[[{4,6}]]]};
             {pm, pp} = (Mean[juncijjj[$fathers[[a]], b][[{1,3} + #]]] & /@ {0,1})+temp,
             If[ $genders[[a]] == 1,
                 {pm, pp} = juncijjj[$fathers[[a]], b][[;; 2]],
                 {pm, pp} = juncijjj[$fathers[[a]], b][[3 ;;]]
             ];
         ];
         {mm, pp} = Boole[a=!=b] {mm, pp};
         {mm, mp, pm, pp},
         a < b,
         Transpose[#]& /@ juncijii[b, a][[{1, 3, 2, 4}]]
         ]
    ]

juncijii[a_?(#<=$nfounder&), b_?(#<=$nfounder&)] :=
    juncijii[a, b] =  Table[SparseArray[ConstantArray[N[0], {$nfgllab, $nfgllab}]], {4}]    
juncijii[a_, b_] :=
    juncijii[a, b] =
    Module[ {mm, mp, pm, pp},
        Which[
         a >= b,
         {mm, mp} = (Mean[juncijii[$mothers[[a]], b][[{0, 2} + #]]] & /@Range[2]);
         If[ $isautosome,
             {pm, pp} = (Mean[juncijii[$fathers[[a]], b][[{0, 2} + #]]] & /@Range[2]),
             If[ $genders[[a]] == 1,
                 {pm, pp} = juncijii[$fathers[[a]], b][[;; 2]],
                 {pm, pp} = juncijii[$fathers[[a]], b][[3 ;;]]
             ];
         ];
         {mm, pp} = Boole[a=!=b] {mm, pp};
         {mm, mp, pm, pp},
         a < b,
         Transpose[#]& /@ juncijjj[b, a][[{1, 3, 2, 4}]]
         ]
    ];
           
clearAncestryMemory[] :=
    Module[ {},
        removeDownValues[fi[_?IntegerQ]];
        removeDownValues[phiij[_?IntegerQ, _?IntegerQ]];
        removeDownValues[phiijk[_?IntegerQ, _?IntegerQ,_?IntegerQ]];
        removeDownValues[Rij[_?IntegerQ]];
        removeDownValues[junciijj[_?IntegerQ, _?IntegerQ]];
        removeDownValues[juncijjj[_?IntegerQ, _?IntegerQ]];
        removeDownValues[juncijii[_?IntegerQ, _?IntegerQ]];
        removeDownValues[juncijkj[_?IntegerQ, _?IntegerQ]];
        removeDownValues[juncijik[_?IntegerQ, _?IntegerQ]];
    ]
                             
clearAncestryConciseMemory[] :=
    Module[ {},
        removeDownValues[fi[_?IntegerQ]];
        removeDownValues[phiij[_?IntegerQ, _?IntegerQ]];
        removeDownValues[Rij[_?IntegerQ]];
        removeDownValues[junciijj[_?IntegerQ, _?IntegerQ]];
    ]
                                                             
pedAncestrySummary[pedigree_,founderFGL_,isAutosome_,pairlist_] :=
    Module[ {pairlist2,ibdij,fab,Rab,Jjijj,Jikjk,Jkikj, Jiiij,Jiiji,Jiijj, Jijii, Jijik, Jijjj, Jijkj,rho,ls,res},
        If[ Count[pedigree[[2 ;;, 4]], {0, 0}]=!=Length[founderFGL],
            Print["Inconsistent number of founders between pedigree and founderFGL!"];
            Abort[]
        ];
        clearAncestryMemory[];        
        (*setstage calculate internal package variables: $isautosome,$mothers,$fathers,$genders*)
        pairlist2 = setstage[pedigree,founderFGL,isAutosome,pairlist];
        {fab,ibdij,Rab,Jiijj,Jijjj,Jijii,Jijkj,Jijik} = {Map[fi, pairlist2, {2}], phiij @@ # & /@ pairlist2, Map[Rij, pairlist2, {2}], junciijj@@ # & /@ pairlist2,
            juncijjj@@ # & /@ pairlist2,juncijii@@ # & /@ pairlist2, juncijkj@@ # & /@ pairlist2,juncijik@@ # & /@ pairlist2};
        rho = Transpose[Join[Rab[[All, 1, 1]] + # & /@ Transpose[Rab[[All, 2]]],Rab[[All, 1, 2]] + # & /@ Transpose[Rab[[All, 2]]]]] -  Jiijj;
        clearAncestryMemory[];
        Jiiij = Jijii;
        Jiiji = Map[Transpose, Jijjj, {2}];
        ls = Transpose[{Rab[[All, 1, 1]], Rab[[All, 1, 1]], Rab[[All, 1, 2]], Rab[[All, 1, 2]]}];
        Jikjk = Map[Transpose, Jijkj, {3}];
        If[ ls!=Jiijj+Jiiji+Jijjj+Map[Total, Jikjk, {4}],
            Print["Junction densities are inconsistent with R_a!"]
        ];
        ls = Transpose[{Rab[[All, 2, 1]], Rab[[All, 2, 2]], Rab[[All, 2, 1]], Rab[[All, 2, 2]]}];
        Jkikj = Map[Transpose[#, {3, 1, 2}] &, Jijik, {2}];
        Jjijj = Map[Transpose, Jiiij, {2}];
        If[ ls !=Jiijj+Jiiij+Jjijj+Map[Total, Jkikj, {4}],
            Print["pedAncestrySummary: junction densities are inconsistent with R_b!"]
        ];
        res = Transpose[{pairlist,fab[[All, 1]],fab[[All, 2]],ibdij,Rab[[All,1]],Rab[[All,2]],rho,Jiiij,Jiiji,Jiijj, Jijii, Jijik, Jijjj, Jijkj}];
        res = Join[{{"{a,b}", "f(i)_a^{m,p}","f(i)_b^{m,p}","phi(ij)_{ab}^{mm,mp,pm,pp}","R(ij)_a^{m,p}","R(ij)_b^{m,p}","rho(ij)_{ab}^{mm,mp,pm,pp}",
             "Jiiij_{ab}^{mm,mp,pm,pp}","Jiiji_{ab}^{mm,mp,pm,pp}","Jiijj_{ab}^{mm,mp,pm,pp}","Jijii_{ab}^{mm,mp,pm,pp}",
             "Jijik_{ab}^{mm,mp,pm,pp}","Jijjj_{ab}^{mm,mp,pm,pp}","Jijkj_{ab}^{mm,mp,pm,pp}"}}, res];
        res
    ]

pedAncestryConciseSummary[pedigree_,founderFGL_,isAutosome_,pairlist_] :=
    Module[ {pairlist2,fab,ibdij, Rab,Jiijj,rho,res},
        If[ Count[pedigree[[2 ;;, 4]], {0, 0}]=!=Length[founderFGL],
            Print["Inconsistent number of founders between pedigree and founderFGL!"];
            Abort[]
        ];
        clearAncestryConciseMemory[];        
        (*setstage calculate internal package variables: $isautosome,$mothers,$fathers,$genders*)
        pairlist2 = setstage[pedigree,founderFGL,isAutosome,pairlist];
        {fab,ibdij,Rab, Jiijj} = {Map[fi, pairlist2, {2}], phiij @@ # & /@ pairlist2,Map[Rij, pairlist2, {2}], junciijj@@ # & /@ pairlist2};
        rho = Transpose[Join[Rab[[All, 1, 1]] + # & /@ Transpose[Rab[[All, 2]]],Rab[[All, 1, 2]] + # & /@ Transpose[Rab[[All, 2]]]]] -  Jiijj;
        clearAncestryConciseMemory[];
        res = Transpose[{pairlist,fab[[All,1]],fab[[All,2]],ibdij,Rab[[All,1]],Rab[[All,2]],rho}];
        res = Join[{{"{a,b}", "f(i)_a^{m,p}","f(i)_b^{m,p}","phi(ij)_{ab}^{mm,mp,pm,pp}","R(ij)_a^{m,p}","R(ij)_b^{m,p}","rho(ij)_{ab}^{mm,mp,pm,pp}"}}, res];
        res
    ]
    
pedAncestryPrior[pedigree_, founderFGL_, isAutosome_, indlist_List] :=
    Module[ {ped,pairlist, res, lab, pos,pos2,pos3,str},
        pos = Flatten[Position[pedigree[[2 ;;, 2]], #, {1}, 1, Heads -> False] & /@indlist];
        If[ pedigree[[pos + 1, 2]] == indlist,
            ped = pedigree[[pos + 1, ;;4]],
            Print["pedIdentityPrior: Wrong indlist!"];
            Abort[];
        ];
        pairlist = Transpose[{indlist, indlist}];
        res = pedAncestrySummary[pedigree, founderFGL, isAutosome, pairlist];
        lab = {"phi(ij)_{aa}^mp","rho(ij)_{aa}^mp", "Jiiij_{aa}^mp",  "Jiiji_{aa}^mp",  "Jiijj_{aa}^mp", "Jijii_{aa}^mp", "Jijik_{aa}^mp", "Jijjj_{aa}^mp", "Jijkj_{aa}^mp"};
        pos = Flatten[Table[Position[res[[1]], _?(StringMatchQ[#, str ~~ ___] &), {1},Heads -> False], {str, StringTake[#, 5] & /@ lab}]];
        If[ Length[pos] != Length[lab],
            Print["pedAncestryPrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        pos2 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "f(i)_a" ~~ ___] &), {1},Heads -> False]];
        If[ Length[pos2] != 1,
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        pos3 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "R(ij)_a"~~ ___] &), {1},Heads -> False]];
        If[ Length[pos3] != 1,
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        res = Join[{Join[pedigree[[1, ;; 4]], {"f(i)_a^m","f(i)_a^p","R(ij)_a^m", "R(ij)_a^p"}, lab]}, 
           Join[ped, res[[2 ;;, pos2[[1]]]],res[[2 ;;, pos3[[1]]]],res[[2 ;;, pos, 2]], 2]];
        res[[All, 7;;9]] = res[[All, {9,7,8}]];
        res
    ]

pedAncestryConcisePrior[pedigree_, founderFGL_, isAutosome_, indlist_List] :=
    Module[ {ped,pairlist, res, lab, pos,pos2,pos3,str},
        pos = Flatten[Position[pedigree[[2 ;;, 2]], #, {1}, 1, Heads -> False] & /@indlist];
        If[ pedigree[[pos + 1, 2]] == indlist,
            ped = pedigree[[pos + 1, ;;4]],
            Print["pedIdentityPrior: Wrong indlist!"];
            Abort[];
        ];
        pairlist = Transpose[{indlist, indlist}];
        res = pedAncestryConciseSummary[pedigree, founderFGL, isAutosome, pairlist];
        lab = {"phi(ij)_{aa}^mp","rho(ij)_{aa}^mp"};
        pos = Flatten[Table[Position[res[[1]], _?(StringMatchQ[#, str ~~ ___] &), {1},Heads -> False], {str, StringTake[#, 5] & /@ lab}]];
        If[ Length[pos] != Length[lab],
            Print["pedIdentityPrior: unexpected results obtained from pedAncestryConcisePrior!"];
            Abort[]
        ];
        pos2 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "f(i)_a" ~~ ___] &), {1},Heads -> False]];
        If[ Length[pos2] != 1,
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        pos3 = Flatten[Position[res[[1]], _?(StringMatchQ[#, "R(ij)_a"~~ ___] &), {1},Heads -> False]];
        If[ Length[pos3] != 1,
            Print["pedIdentityConcisePrior: unexpected results obtained from pedIdentitySummary!"];
            Abort[]
        ];
        res = Join[{Join[pedigree[[1, ;; 4]], {"f(i)_a^m","f(i)_a^p","R(ij)_a^m", "R(ij)_a^p"}, lab]}, 
           Join[ped, res[[2 ;;, pos2[[1]]]],res[[2 ;;, pos3[[1]]]],res[[2 ;;, pos, 2]], 2]];
        res[[All, 7;;9]] = res[[All, {9,7,8}]];
        res
    ]


ancestryJointDensity[junctions_] :=
    Module[ {Jiiij, Jiiji, Jiijj, Jijii, Jijik, Jijjj, Jijkj, nFgl, ls0, 
      density, ii, ix, xi, z, x},
        {Jiiij, Jiiji, Jiijj, Jijii, Jijik, Jijjj, Jijkj} = junctions;
        nFgl = Length[Jiiij];
        ls0 = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1];
        (*rate=Outer[List,ls0,ls0,1,1];*)
        density = SparseArray[{}, {Length[ls0], Length[ls0]}];
        ii = Flatten[Position[ls0, {x_, x_}]];
        ix = Flatten[Position[ls0, {#, x_}]] & /@ Range[nFgl];
        xi = Transpose[ix];
        Do[density[[ix[[z]], ix[[z]]]] = Jijik[[z]], {z, nFgl}];
        Do[density[[xi[[z]], xi[[z]]]] = Jijkj[[All, z]], {z, nFgl}];
        Do[density[[ii[[z]], ix[[z]]]] = Jiiij[[z]], {z, nFgl}];
        Do[density[[ii[[z]], xi[[z]]]] = Jiiji[[z]], {z, nFgl}];
        Do[density[[ix[[z]], ii[[z]]]] = Jijii[[z]], {z, nFgl}];
        Do[density[[xi[[z]], ii[[z]]]] = Jijjj[[z]], {z, nFgl}];
        density[[ii, ii]] = Jiijj;
        density
    ]

ancestryRateMatix[initial_List, densityMtx_] :=
    Module[ {rateQ, pos, pos2, z},
        rateQ = densityMtx;
        pos = Flatten[Position[initial, _?Positive]];
        rateQ[[pos]] /= initial[[pos]];
        pos2 = Complement[Range[Length[initial]], pos];
        If[ ! (Total[Flatten[rateQ[[pos2]]]] == 0),
            Print["ancestryMarkovProcess: Inconsistent results between ancestry coefficents and junction densities!"]
        ];
        If[ ! (Total[Diagonal[rateQ]] == 0),
            Print["ancestryMarkovProcess: Nonzero diasognal junction desensities!"]
        ];
        Do[rateQ[[z, z]] = -Total[rateQ[[z]]], {z, Length[rateQ]}];
        rateQ
    ]

stationaryMarkovProcessQ[initial_, rateMtx_] :=
    Abs[Total[initial.rateMtx]] < 10^(-10.)
 
pedAncestryPrior[pedigree_, founderFGL_, indlist_, isJointModel_?(MatchQ[#, True | False] &)] :=
    Module[ {isOogamy, isautosome,isautosomels, res, ls, initial, densityMtx,i,k},
        isOogamy = Total[Union[pedigree[[2 ;;, 3]]]] != 0;
        isautosomels = If[ isOogamy,
                           {True, False},
                           {True}
                       ];
        res = ConstantArray[0, {Length[indlist] + 1, 5 + Boole[isOogamy]}];
        If[ isJointModel,
            Do[
             (*ls: {"Generation","MemberID","Female=1/Male=2/hermaphrodite=0",{"MotherID", "FatherID"},
             "f(i)_a^m","f(i)_a^p","phi(ij)_{aa}^mp","R(ij)_a^m","R(ij)_a^p",
             "rho(ij)_{aa}^mp","Jiiij_{aa}^mp","Jiiji_{aa}^mp","Jiijj_{aa}^mp",
             "Jijii_{aa}^mp","Jijik_{aa}^mp","Jijjj_{aa}^mp","Jijkj_{aa}^mp"}*)
             ls = pedAncestryPrior[pedigree, founderFGL, isautosome, indlist];
             res[[2 ;;, If[ isautosome,
                            5,
                            6
                        ]]] = Table[
                If[ ls[[i,3]]==2&&(!isautosome),
                    initial = Normal[ls[[i, 5]]];
                    densityMtx = ls[[i, 8]],                        
                    (*inital=phiij of dimension 1xL^2, densityMtx of dimension L^2xL^2*)
                    initial = Normal[Flatten[ls[[i, 7]]]];
                    densityMtx = ancestryJointDensity[ls[[i, 11 ;;]]];
                ];
                {SparseArray[initial], ancestryRateMatix[initial, densityMtx]}, {i, 2, Length[ls]}];
             0, {isautosome, isautosomels}];
            res[[All, ;; 4]] = ls[[All, ;; 4]];
            res[[1, 5]] = {"initial_AA", "rateMtx_AA"};
            If[ isOogamy,
                res[[1, 6]] = {"initial_sex", "rateMtx_sex"};
            ],
            Do[
             (*ls: {"Generation","MemberID","Female=1/Male=2/hermaphrodite=0",
             "f(i)_a^m","f(i)_a^p","phi(ij)_{aa}^mp","R(ij)_a^m","R(ij)_a^p",
             "rho(ij)_{aa}^mp"}*)
             ls = pedAncestryConcisePrior[pedigree, founderFGL, isautosome, indlist];
             res[[2 ;;, If[ isautosome,
                            5,
                            6
                        ]]] = Table[
               (*inital=fi of dimension 1xL,densityMtx=Rij of dimension LxL*)
               initial = Normal[ls[[i, 5 + k]]];
               densityMtx = ls[[i, 8 + k]];
               {SparseArray[initial], ancestryRateMatix[initial, densityMtx]}, {i, 2, Length[ls]}, {k, {0, 1}}];
             0, {isautosome, isautosomels}];
            res[[All, ;; 4]] = ls[[All, ;; 4]];
            res[[1,5]] = {{"initial_A(m)", "rateMtx_A(m)"}, {"initial_A(p)","rateMtx_A(p)"}};
            If[ isOogamy,
                res[[1,6]] = {{"initial_sex(m)", "rateMtx_sex(m)"}, {"initial_sex(p)","rateMtx_sex(p)"}};
            ];
        ];
        If[ ! (And @@ Flatten[Map[stationaryMarkovProcessQ @@ # &,res[[2 ;;, 5 ;;]], {Depth[res[[2 ;;, 5 ;;]]] - 4}]]),
            Print["pedAncestryPrior: the initial is a non-stationary distribution!"]
        ];
        res
    ]
          
pedAncestryCompatibleQ[pedigree_] :=
    Module[ {isOogamy, nFounder, founderFGL, typels, type, inds, 
      isautosome, residentity, resancestry,j},
        isOogamy = (Union[pedigree[[2 ;;, 3]]] =!= {0});
        nFounder = Count[pedigree[[2 ;;, -1]], {0, 0}];
        founderFGL = Transpose[{Range[nFounder], Range[nFounder]}];
        typels = If[ isOogamy,
                     {"AA", "XX", "XY"},
                     {"AA"}
                 ];
        Table[
         inds = Switch[type,
                    "AA",
                    pedigree[[2 ;;, 2]],
                    "XX",
                    Select[pedigree[[2 ;;, ;; 3]], #[[3]] == 1 &][[All, 2]],
                    "XY",
                    Select[pedigree[[2 ;;, ;; 3]], #[[3]] == 2 &][[All, 2]]
                    ];
         isautosome = If[ type == "AA",
                          True,
                          False
                      ];
         residentity = pedIdentityPrior[pedigree, founderFGL, isautosome, inds];
         (*{"Generation","phi(11)_{aa}^mp","R_a^m","R_a^p","rho_{aa}^{mp}",
         "J1112_{aa}^mp","J1121_{aa}^mp","J1122_{aa}^mp","J1211_{aa}^mp",
         "J1213_{aa}^mp","J1222_{aa}^mp","J1232_{aa}^mp"}*)
         residentity = residentity[[All, Prepend[Range[5, Dimensions[residentity][[2]]], 1]]];
         resancestry = pedAncestryPrior[pedigree, founderFGL, isautosome, inds];
         (*{"Generation","phi(ij)_{aa}^mp","R(ij)_a^m","R(ij)_a^p",
         "rho(ij)_{aa}^mp","Jiiij_{aa}^mp","Jiiji_{aa}^mp","Jiijj_{aa}^mp",
         "Jijii_{aa}^mp","Jijik_{aa}^mp","Jijjj_{aa}^mp","Jijkj_{aa}^mp"}*)
         resancestry = resancestry[[All, Prepend[Range[7, Dimensions[resancestry][[2]]], 1]]];
         resancestry[[2 ;;, 2]] = Total[Diagonal[#]] & /@ resancestry[[2 ;;, 2]];
         Do[resancestry[[2 ;;, j]] = Total[Flatten[#]] & /@ resancestry[[2 ;;, j]], {j, 3, Dimensions[resancestry][[2]]}];
         resancestry[[2 ;;]] == residentity[[2 ;;]], {type, typels}]
    ]       
       
End[]

EndPackage[]