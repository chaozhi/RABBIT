(* Mathematica Package *)

BeginPackage["MagicReconstruct`MagicLikelihood`",{"MagicDataPreprocess`"}]
(* Exported symbols added here with SymbolName::usage *)  

(*For called genotype data from array chip or GBS*)

lineMagicLikelihood::usage = "lineMagicLikelihood  "

siteMagicLikelihood::usage = "siteMagicLikelihood  "

siteParentHaploPrior::usage = "siteParentHaploPrior  "

posteriorTrueDiplo::usage = "posteriorTrueDiplo  "

posteriorTrueHaplo::usage = "posteriorTrueHaplo  "

(*For allelic depth data of GBS*)

lineMagicLikelihoodGBS::usage = "lineMagicLikelihoodGBS  "

siteMagicLikelihoodGBS::usage = "siteMagicLikelihoodGBS  "

siteParentHaploPriorGBS::usage = "siteParentHaploPriorGBS  "

Begin["`Private`"] (* Begin Private Context *) 

(*diploPrior returns P(Z|D, O, epsF) at a given site, 
 where Z=true genotype, D=derived genotype,
 O= ancestral origins, epsF = founder allelic error probability*)
(*likelihoodDiplo returns l=Pr(Y|D(H,O), O, eps,epsF) the probability of observed Y given the derived D and IBD status*)
(*There are only two alleles 1 and 2, and missing alllels are denoted by N(Null).*)
(*likeli, dimension 6 x4, refers to the probability of observed genotype given latent genotpe*)
(*The 6 observed genotype is defined as {"NN" -> 1, "N1" -> 2, "1N" -> 2, "N2" -> 3, "2N" -> 3, "11" -> 4, "12" -> 5, "21" -> 5, "22" -> 6}     *)
(*likeliFibd0 or likeliFibd1, dimension 9 x4, refers to the probability of derived genotype given latent genotype*)
(*The 9 derived genotypes are defined as {"NN" -> 1, "N1" -> 2, "1N" -> 3, "N2" -> 4, "2N" -> 5, "11" -> 6, "12" -> 7, "21" -> 8, "22" -> 9}*)
(*assuming that prior for latent genotype (1,1), (1,2), (2,1), and (22) are equally probable*)    
diploPrior[epsF_] :=
    Module[ {likeliF,likeliFibd0, likeliFibd1},
        likeliFibd0 = {{1, 1, 1, 1}, 
                       {1 - epsF, epsF, 1 - epsF, epsF}, 
                       {1 - epsF, 1 - epsF, epsF, epsF}, 
                       {epsF, 1 - epsF, epsF, 1 - epsF}, 
                       {epsF, epsF, 1 - epsF, 1 - epsF}, 
                       {(1 - epsF)^2, (1 - epsF) epsF, (1 - epsF) epsF, epsF^2}, 
                       {(1 - epsF) epsF, (1 - epsF)^2, epsF^2, (1 - epsF) epsF}, 
                       {(1 - epsF) epsF, epsF^2, (1 -epsF)^2, (1 - epsF) epsF}, 
                       {epsF^2, (1 - epsF) epsF, (1 - epsF) epsF, (1 - epsF)^2}};
        likeliFibd1 = {{1, 0, 0, 1}, {0, 0, 0, 0}, {0, 0, 0, 0}, 
                       {0, 0, 0, 0}, {0, 0, 0, 0}, {1 - epsF, 0, 0, epsF}, 
                       {0, 0, 0, 0}, {0, 0, 0, 0}, {epsF, 0, 0, 1 - epsF}};
        Simplify[Table[Transpose[Normalize[#, Total] & /@ likeliF], {likeliF, {likeliFibd0, likeliFibd1}}]]
    ]

diploLikeli[eps_] :=
    {{1, 1, 1, 1}, {1 - eps, 1/2, 1/2, eps}, {eps, 1/2, 1/2, 1 - eps}, 
    {(1 - eps)^2, (1 - eps) eps, (1 - eps) eps, eps^2}, 
    {2 (1 - eps) eps, (1 - eps)^2 + eps^2, (1 - eps)^2 + eps^2, 2 (1 - eps) eps}, 
    {eps^2, (1 - eps) eps, (1 - eps) eps, (1 - eps)^2}}

marginLikeliDiplo[epsF_, eps_] :=
    Module[ {likeli, prior},
        likeli = diploLikeli[eps];
        prior = diploPrior[epsF];
        Transpose[likeli.#&/@prior]
    ]

(*diploPrior returns P(Z|D, O, epsF) at a given site, where Z=true genotype, D=derived genotype, 
O= ancestral origins, epsF = founder allelic error probability*)
(*likelihoodHaplo returns l=Pr(Y|D(H,O), O, eps,epsF) the probability of observed Y given the derived D and IBD status*)
(*There are only two alleles 1 and 2, and missing alllels are denoted by N(Null).*)
(*likeli, dimension 2x3, refers to the probability of observed allele given latent true allele*)
(*The 3 observed allele is defined as {"N" -> 1, "1" -> 2, "2" -> 3}*)
(*likeliF, dimension 2x3,, refers to the probability of derived allele given latent true allele*)
(*The 3 derived genotypes are defined as {"N" -> 1, "1" -> 2, "2" -> 3}*)
(*assuming that prior for the two latent true alleles  1 and 2 are equally probable*)
haploPrior[epsF_] :=
    Module[ {likeliF,prior},
        likeliF = {{1, 1},{1-epsF,epsF},{epsF,1-epsF}};
        prior = Transpose[Normalize[#, Total] & /@ likeliF];
        prior
    ]

haploLikeli[eps_] :=
    {{1, 1},{1-eps,eps},{eps,1-eps}}

marginLikeliHaplo[epsF_, eps_] :=
    Module[ {likeli, prior},
        likeli = haploLikeli[eps];
        prior = haploPrior[epsF];
        likeli.prior
    ]

    
(*y=observed genotype, z=true genotype, o=true origin state, D=derived genotype*)
(*likeli dimensions {6,4} ={ny,nz}*)
(*prior dimesnions {2, 4,9} ={nibdstate (non-ibd, ibd),nz,nD}*)
(*marginal dimensions {6,2,9} ={ny, nibdstate, nD}}, margin = P(y|o,D)*)
(*posterior = P(z|y,o,D)*)
(*posterior dimensions {4,6,2,9} = {nz,ny,nibdstate,nD}*)
diploJointPosterior[epsF_, eps_] :=
    Module[ {likeli, prior, marginal, posterior,ibd,z},
        likeli = diploLikeli[eps];
        prior = diploPrior[epsF];
        marginal = Transpose[likeli.# & /@ prior];
        marginal = marginal /. {0 | 0. -> Indeterminate};
        posterior = Transpose[#] & /@ Table[KroneckerProduct[likeli[[All, z]], prior[[ibd, z]]], {z, Dimensions[likeli][[2]]}, {ibd, Length[prior]}];
        posterior = #/marginal & /@ posterior;
        posterior /. {Indeterminate -> 0}
    ]    

(*returns data probability (likelohood) for one sampled indivdiual(magic line) at all SNP sites*)
lineMagicLikelihood[model_,founderHaplo_, obsgeno_, epsF_, eps_, posA_, posX_, gender_] :=
    Module[ {likelidiplo,likelihaplo, eps2 = eps,posDiplo, posHaplo, nFounder, states, ibd12, 
      flatstates, derivedRule, obsRule, obsSeq, dataProb, hap, der, i, j},
        If[ eps==0,
            eps2 = 10^(-10.)
        ];
        likelidiplo = marginLikeliDiplo[epsF, eps2];
        likelihaplo = marginLikeliHaplo[epsF, eps2];
        {posDiplo, posHaplo} = If[ gender == "Male",
                                   {posA, posX},
                                   {Join[posA, posX], {}}
                               ];
        (*There are only two alleles 1 and 2,and missing alllels are denoted by N(Null).*)
        (*The second part of obsRule is for the X chromosome of a male*)
        obsRule = Join[{"NN" -> 1, "N1" -> 2, "1N" -> 2, "N2" -> 3, "2N" -> 3, 
           "11" -> 4, "12" -> 5, "21" -> 5, "22" -> 6}, {"N" -> 1, "1" -> 2,"2" -> 3}];
        obsSeq = Replace[obsgeno, obsRule, {2}];
        dataProb = Table[0, {Length[posDiplo] + Length[posHaplo]}];
        nFounder = Length[founderHaplo];
        (*caluclate dataProb for diplo*)
        If[ posDiplo=!={},
            If[MatchQ[ToLowerCase[model], "depmodel"],
            	states = Transpose[Table[Range[nFounder], {2}]],
            	states = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1]
            ];
            ibd12 = Boole[Equal @@ # & /@ states] + 1;
            flatstates = Flatten[states];
            derivedRule = {{"N", "N"} -> 1, {"N", "1"} -> 2, {"1", "N"} -> 
               3, {"N", "2"} -> 4, {"2", "N"} -> 5, {"1", "1"} -> 
               6, {"1", "2"} -> 7, {"2", "1"} -> 8, {"2", "2"} -> 9};
            dataProb[[posDiplo]] = Table[
              hap = founderHaplo[[All, i, j]];
              der = Replace[Partition[hap[[flatstates]], 2], derivedRule, {1}];
              Extract[likelidiplo[[obsSeq[[i, j]]]], Transpose[{ibd12, der}]]
              , {i, posDiplo}, {j, Dimensions[founderHaplo[[All, i]]][[2]]}];
        ];
        (*caluclate dataProb for haplo*)
        If[ posHaplo=!={},
            dataProb[[posHaplo]] = Table[
                der = founderHaplo[[All, i, j]] /. {"N" -> 1, "1" -> 2, "2" -> 3};
                likelihaplo[[obsSeq[[i, j]], der]], {i, posHaplo}, {j, Dimensions[founderHaplo[[All, i]]][[2]]}];
        ];
        dataProb
    ]   
    

getDiploLikeliRule[epsF_] :=
    Module[ {pair, ww, eps, count, countrule, weightrule},
        (*pair of true diplotype and obsevered genotype*)
        pair = Outer[List, {"11", "12", "21", "22"}, {"NN", "1N", "2N", "11", "12", "22"}];
        ww = Transpose[diploLikeli[eps]];
        count = 1 - Sign[ww /. eps -> 0];
        countrule = Thread[Flatten[pair, 1] -> Flatten[count]];
        ww = ww /. eps -> epsF;
        weightrule = Thread[Flatten[pair, 1] -> Flatten[ww]];
        {countrule, weightrule}
    ]

(*return a prior for imputation HMM for parental haplotype at a given site
which is calcuated as a posterior after accouting for parental data at the site. *) 
siteParentHaploPrior[obsHt_, epsF_,isfoundermaleX_,isfounderinbred_,maxfoundererror_] :=
    Module[ {posmale,missing = "N", alleles = {"1", "2"},epsF2,maxfoundererror2,hh, ww,hh2,obsHt2,countmissing,countrule, weightrule,count,pos},
        posmale = Flatten[Position[isfoundermaleX, True]];
        epsF2 = Sign[maxfoundererror] epsF;
        maxfoundererror2 = Sign[epsF] maxfoundererror;
        If[ isfounderinbred,
            If[ epsF2==0,
                hh = Distribute[obsHt /. {missing -> alleles}, List];
                If[ (!isfounderinbred)&& (posmale=!={}),
                    hh[[All, 2 posmale - 1]] = hh[[All, 2 posmale]];
                    hh = DeleteDuplicates[hh];
                ];
                ww = Table[1./Length[hh], {Length[hh]}],
                hh = Distribute[Table[alleles, {Length[obsHt]}], List];
                If[ (!isfounderinbred)&& (posmale=!={}),
                    hh[[All, 2 posmale - 1]] = hh[[All, 2 posmale]];
                    hh = DeleteDuplicates[hh];
                ];
                countmissing = Count[obsHt, missing];
                count = Length[obsHt] - countmissing - (Count[# - obsHt, 0] & /@ hh);
                pos = Flatten[Position[count, _?(# <= maxfoundererror2 &), {1}, Heads -> False]];
                count = count[[pos]];
                ww = Normalize[(1 - epsF2)^(Length[obsHt]-countmissing-count) epsF2^count, Total];
                hh = hh[[pos]];
            ],
            hh = Distribute[Table[alleles, {Length[obsHt]}], List];
            If[ (!isfounderinbred)&& (posmale=!={}),
                hh[[All, 2 posmale - 1]] = hh[[All, 2 posmale]];
                hh = DeleteDuplicates[hh];
            ];
            hh2 = Map[StringJoin @@ # &, Partition[#, 2] & /@ hh, {2}];
            obsHt2 = StringJoin @@ # & /@ (Sort[#] & /@ Partition[obsHt, 2]);
            hh2 = Transpose[{#, obsHt2}] & /@ hh2;
            {countrule, weightrule} = getDiploLikeliRule[epsF2];
            count = Total[Replace[hh2, countrule, {2}], {2}];
            pos = Flatten[Position[count, _?(# <= maxfoundererror2 &), {1}, Heads -> False]];
            ww = Replace[hh2[[pos]], weightrule, {2}];
            ww = Normalize[Times @@ # & /@ ww,Total];
            hh = hh[[pos]];
        ];
        Transpose[{hh, ww}]
    ]
          
(*sitefhaplo dimenions: {ngeno,nfgl}, ngeno= number of possible founder haplotypes at a given site*) 
(*Set epsF=0 when sitefhaplo refers to the set of lattent haplotypes.*)
(*return dataprob with diemsnions {ngeno, noffspring,nstate}*)
(*returns data probability (likelohood) for all sampled invidual at a single SNP site, 
given each of the possible founder haplotype at that site*)

siteMagicLikelihood[model_,sitefhaplo_, sitegeno_, epsF_,eps_, ismaleX_] :=
    Module[ {eps2 = eps, obsRule, obsSeq, nfgl, noffspring,states, 
      ibd12, derivedRule, likelidiplo, likelihaplo,maleX, nonmaleX, dataprob, der, i},
        If[ eps == 0,
            eps2 = 10^(-10.)
        ];
        (*There are only two alleles 1 and 2,and missing alllels are denoted by N(Null).*)
        (*The second part of obsRule is for the X chromosome of a male*)
        obsRule = Join[{"NN" -> 1, "N1" -> 2, "1N" -> 2, "N2" -> 3, "2N" -> 3, 
            "11" -> 4, "12" -> 5, "21" -> 5, "22" -> 6}, {"N" -> 1, "1" -> 2, "2" -> 3}];
        obsSeq = Replace[sitegeno, obsRule, {1}];
        nfgl = Length[First[sitefhaplo]];
        noffspring = Length[sitegeno];
        If[ Length[ismaleX] != noffspring,
            Print["siteMagicLikelihood: The number of offspring = ", noffspring, " is unequal to length of ismaleX = ", Length[ismaleX],"!"];
            Abort[]
        ];
        likelihaplo = marginLikeliHaplo[epsF, eps2];
        likelidiplo = marginLikeliDiplo[epsF, eps2];
        maleX = Flatten[Position[ismaleX, True]];
        nonmaleX = Complement[Range[noffspring], maleX];
        dataprob = ConstantArray[0, {Length[sitefhaplo], noffspring}];
        dataprob[[All, maleX]] = Table[
           der = sitefhaplo[[i]] /. {"N" -> 1, "1" -> 2, "2" -> 3};
           likelihaplo[[obsSeq[[maleX]], der]], {i, Length[sitefhaplo]}];        
        If[MatchQ[ToLowerCase[model], "depmodel"],
        	states = Transpose[Table[Range[nfgl], {2}]],
        	states = Flatten[Outer[List, Range[nfgl], Range[nfgl]], 1]
        ];        
        ibd12 = Boole[Equal @@ # & /@ states] + 1;
        states = Flatten[states];
        derivedRule = {{"N", "N"} -> 1, {"N", "1"} -> 2, {"1", "N"} -> 
            3, {"N", "2"} -> 4, {"2", "N"} -> 5, {"1", "1"} -> 
            6, {"1", "2"} -> 7, {"2", "1"} -> 8, {"2", "2"} -> 9};
        dataprob[[All, nonmaleX]] = Table[
            der = Replace[Partition[sitefhaplo[[i, states]], 2], derivedRule, {1}];
            Extract[#, Transpose[{ibd12, der}]] & /@likelidiplo[[obsSeq[[nonmaleX]]]], {i, Length[sitefhaplo]}];
        dataprob
    ]
    
    
(*ancestralgenoprob dimensions: {nind,nchr,nsnp,nstate=nfgl (nfgl+1)/2}*)    
(*to calculate diploPosterior, dimensions {nind, nchr, nsnp, ntruegenotype=4}*)
posteriorTrueDiplo[ancestralgenoprob_, derivedgeno_, genotypes_, epsF_, eps_] :=
    Module[ {nind, nchr, posibd, posnonibd, prior,priornonibd,prioribd,
      nonibdAncestral, ibdAncestral, nonibdDerived, ibdDerived, res,i,ind,ch},
        {nind, nchr} = Take[Dimensions[ancestralgenoprob], 2];
        posibd = Flatten[Position[genotypes, {i_, i_}]];
        posnonibd = Complement[Range[Length[genotypes]], posibd];
        {nonibdAncestral, ibdAncestral} = ancestralgenoprob[[All, All, All, #]] & /@ {posnonibd, posibd};
        {nonibdDerived, ibdDerived} = derivedgeno[[All, All, #]] & /@ {posnonibd, posibd};
        (*prior=P(true genotype|derived genotype,ancestral origin,epsF); prior[[1]] for non-ibd,and prior[[2]] for ibd.*)
		(*Dimensions[prior]={2,4,9},where 4=number of true genotypes,9=number of derived genotypes*)
		(*9 derived genotypes ={NN,N1,1N,N2,2N,11,12,21,22}*)
		prior = diploPrior[epsF];
		(*Dimensions[priornonibd] ={nchr,nsnp,#ibdstate=6,4}*)
		(*Dimensions[prioribd] ={nchr,nsnp,#non-ibdstate=4,4}*)
		priornonibd = Transpose[prior[[1]]];
		priornonibd = Map[priornonibd[[#]] &, nonibdDerived, {2}];
		prioribd = Transpose[prior[[2]]];
		prioribd = Map[prioribd[[#]] &, ibdDerived, {2}];
		Table[
		   res = Total[nonibdAncestral[[ind, ch]] priornonibd[[ch]], {2}];
		   res += Total[ibdAncestral[[ind, ch]] prioribd[[ch]], {2}];
		   res, {ind, nind}, {ch, nchr}]
    ]       


(*ancestralhaploprob dimensions: {nind,nchr,nsnp,nstate=nFGL}*)
(*derivedhaplo dimensions: {nchr,nsnp,nstate}, dervied genotype correspoindg to each haplotype state*)    
posteriorTrueHaplo[ancestralhaploprob_, derivedhaplo_, epsF_, eps_] :=
    Module[ {nind, nchr, prior,ind,ch},
        {nind, nchr} = Take[Dimensions[ancestralhaploprob], 2];
        prior = Transpose[haploPrior[epsF]];
        prior = Map[prior[[#]] &, derivedhaplo, {2}];
        Table[Total[ancestralhaploprob[[ind, ch]] prior[[ch]], {2}], {ind,nind}, {ch, nchr}]
    ]   
    
(**********************Likelihood for aleleic depth data from genotyping-by-sequecing****************************)

(*scenario: for all chromosomes and all snps of a single offspring*)
(*observedgeno dimensions {nchr, nsnp, nallele=2}*)
(*return dimensions {nchr,nsnp, ntruegenotype=4}*)
diploLikeliGBS[observedgeno_?(Depth[#] == 4 &), eps_,minphredscore_] :=
    Module[ {baseerrorprob,n1, n2, ntot, p11, p12, p21, p22,pp11, pp12, pp21, pp22,epsls},
        baseerrorprob = 10^(-minphredscore/10.);
        n1 = observedgeno[[All, All, 1]];
        n2 = observedgeno[[All, All, 2]];
        ntot = n1 + n2;
        p12 = p21 = (1./2)^ntot;
        p11 = (1 - baseerrorprob)^n1 baseerrorprob^n2;
        p22 = (baseerrorprob)^n1 (1 - baseerrorprob)^n2;
        epsls = {(1 - eps)^2, (1 - eps) eps, (1 - eps) eps, eps^2};
        pp11 = Total[{p11, p12, p21, p22} epsls];
        pp12 = Total[{p12, p11, p22, p21} epsls];
        pp21 = Total[{p21, p11, p22, p12} epsls];
        pp22 = Total[{p22, p12, p21, p11} epsls];
        Transpose[#] & /@ Transpose[{pp11, pp12, pp21, pp22}]
    ]

(*scenario: for all individual at a single site*)
(*it applied to both offspring and outbred funders*)
(*observedgeno dimensions {nind, nallele=2}*)
(*return dimensions {nind, ntruegenotype=4}*)
diploLikeliGBS[observedgeno_?(Depth[#] == 3 &), eps_,minphredscore_] :=
    First[diploLikeliGBS[{observedgeno}, eps,minphredscore]]
 
(*likeli dimensions {nchr,nsnp, ntruegenotype=4}*)
(*prior dimesnions {nibdstate=2,ntruegenotype=4,nderivedgenotype=9}*)
(*derivedgeno dimensions {nchr, nsnp, ngenotypes=nFgl^2}*)
(*observedgeno dimensions {nchr, nsnp, nallele=2}*)
(*likelidiplo dimensions {nibdstate=2,nsnp,nderivedgenotype}*)
(*return dimensions {nchr, nsnp,ngenotypes}*)
lineMarginLikeliDiploGBS[model_,derivedgeno_, observedgeno_, epsF_, eps_,minphredscore_] :=
    Module[ {likeli,prior,nFounder,likelidiplo, genotypes, ibd12,i},
        likeli = diploLikeliGBS[observedgeno, eps,minphredscore];
        (*prior[[1]] for non-ibd, and prior[[2]] for ibd. *)
        prior = diploPrior[epsF];
        likelidiplo = Map[{#.prior[[1]], #.prior[[2]]} &, likeli, {2}];
        (*genotypes = origGenotype[nFounder][[1, 2]];*)
        If[MatchQ[ToLowerCase[model], "depmodel"],
        	nFounder = Length[derivedgeno[[1, 1]]];
        	genotypes = Transpose[Table[Range[nFounder], {2}]],
        	nFounder = Sqrt[Length[derivedgeno[[1, 1]]]];
        	genotypes = Flatten[Outer[List, Range[nFounder], Range[nFounder]], 1]
        ];
        ibd12 = Boole[Equal @@ # & /@ genotypes] + 1;
        Table[MapThread[Extract[#1, Transpose[{ibd12, #2}]] &, {likelidiplo[[i]], derivedgeno[[i]]}], {i, Length[likelidiplo]}]
    ]

(*scenario: for all chromosomes and all snps of a single offspring*)
(*Three possible derived haplotype: 1, 2, N*)
(*Two possible true haplotype: 1, 2*)
(*observedgeno dimensions {nchr, nsnp, nallele=2}*)
(*likeli dimensions {nchr,nsnp, ntruehaplotype=2}*)
(*prior dimesnions {ntruegenotype=2,nderivedhaplotype=3}*)
(*return dimensions {nchr, nsnp,nderivedhaplotye=3}*)
haploLikeliGBS[observedhaplo_?(Depth[#] == 4 &), eps_,minphredscore_] :=
    Module[ {baseerrorprob,n1, n2, ntot, p1, p2,pp1,pp2},
        baseerrorprob = 10^(-minphredscore/10.);
        n1 = observedhaplo[[All, All, 1]];
        n2 = observedhaplo[[All, All, 2]];
        ntot = n1 + n2;
        p1 = (1 - baseerrorprob)^n1 baseerrorprob^n2;
        p2 = (baseerrorprob)^n1 (1 - baseerrorprob)^n2;
        pp1 = (1-eps) p1+eps p2;
        pp2 =  (1-eps) p2+eps p1;
        Transpose[{pp1,pp2}, {3, 1, 2}]
    ]

(*scenario: for all offspring at a single site*)
(*observedgeno dimensions {noffspring, nallele=2}*)
(*return dimensions {noffspring, ntruehaplotype=2}*)
haploLikeliGBS[observedhaplo_?(Depth[#] == 3 &), eps_,minphredscore_] :=
    First[haploLikeliGBS[{observedhaplo}, eps,minphredscore]]

(*likelihaplo dimensions {nchr, nsnp,nderivedgenotype=3}*)
(*return dimensions {nchr,nsnp,nfgl}*)
lineMarginLikeliHaploGBS[fhaplo_, observedhaplo_, epsF_, eps_,minphredscore_] :=
    Module[ {likeli,prior,likelihaplo, der,i,j},
        likeli = haploLikeliGBS[observedhaplo, eps,minphredscore];
        prior = haploPrior[epsF];
        likelihaplo = Map[#.prior &, likeli, {2}];
        Table[
         der = fhaplo[[All, i, j]] /. {"N" -> 1, "1" -> 2, "2" -> 3};
         likelihaplo[[i, j, der]], {i, Length[fhaplo[[1]]]}, {j,Length[fhaplo[[1, i]]]}]
    ]

 
lineMagicLikelihoodGBS[model_,founderhaplo_,derivedgeno_, observedgeno_, epsF_, eps_, minphredscore_,posA_, posX_, gender_] :=
    Module[ {posDiplo, posHaplo, dataProb},
        {posDiplo, posHaplo} = 
         If[ gender == "Male",
             {posA, posX},
             {Join[posA, posX], {}}
         ];
        dataProb = Table[0, {Length[posDiplo] + Length[posHaplo]}];
        If[ posDiplo=!={},
            dataProb[[posDiplo]] = lineMarginLikeliDiploGBS[model, derivedgeno[[posDiplo]], observedgeno[[posDiplo]],epsF, eps,minphredscore];
        ];
        If[ posHaplo=!={},
            dataProb[[posHaplo]] = lineMarginLikeliHaploGBS[founderhaplo[[All, posHaplo]],observedgeno[[posHaplo]], epsF, eps,minphredscore];
        ];
        dataProb
    ]    

                
siteParentHaploPriorGBS[obsHt_, epsF_, minphredscore_,maxfoundererror_,genothreshold_,isfoundermaleX_,isfounderinbred_] :=
    Module[ {calledHt,hh, ww, likeli,rule},
        If[ isfounderinbred,
            calledHt = Replace[Sign[obsHt],{{0, 0} -> "N", {1, 0} -> "1", {0, 1} -> "2", {1, 1} -> "N"},{1}];
            likeli = haploLikeliGBS[obsHt, epsF,minphredscore];
            rule = Thread[Range[2] -> {"1", "2"}],
            calledHt = Flatten[Characters[rawGenotypeCall[obsHt, minphredscore,genothreshold]]];
            (*calledHt = Flatten[Replace[Sign[obsHt], {{0, 0} -> {"N", "N"}, {1, 0} -> {"1", "N"}, {0, 1} -> {"2", "N"}, {1, 1} -> {"1", "2"}}, {1}]];*)
            likeli = diploLikeliGBS[obsHt, epsF,minphredscore];
            rule = Thread[Range[4] -> {{"1", "1"}, {"1", "2"}, {"2", "1"}, {"2", "2"}}];
        ];
        hh =  siteParentHaploPrior[calledHt, epsF+10^(-minphredscore/10.), isfoundermaleX,isfounderinbred,maxfoundererror][[All, 1]];
        If[ isfounderinbred,
            hh = Transpose[Replace[hh, Reverse[#] & /@ rule, {2}]],
            hh = Transpose[Replace[Partition[#, 2] & /@ hh, Reverse[#] & /@ rule, {2}]];
        ];
        ww = Times @@ # & /@ Transpose[MapThread[#1[[#2]] &, {likeli, hh}]];
        ww = Normalize[ww, Total];
        hh = Transpose[hh] /.rule;
        hh = Flatten[#] & /@ hh;
        Transpose[{hh, ww}]
    ]

 
(*likeli dimensions {noffspring, ntruegenotype=4}*)
(*prior dimesnions {nibdstate=2,ntruegenotype=4,nderivedgenotype=9}*)
(*likelidiplo dimensions {noffspring,nibdstate=2,nderivedgenotype}*)
(*pos dimensions: {nparentalhaplo, ngenotypes=nfgl^2, 2}}*)
(*return dimensions {nparentalhaplo, noffspring,ngenotypes}*)
siteMarginLikeliDiploGBS[model_,sitefhaplo_, sitegeno_, epsF_, eps_,minphredscore_] :=
    Module[ {likeli, prior, likelidiplo, nFgl, genotypes, ibd12, derivedRule, pos},
        likeli = diploLikeliGBS[sitegeno, eps,minphredscore];
        (*prior[[1]] for non-ibd,and prior[[2]] for ibd.*)
        prior = diploPrior[epsF];
        likelidiplo = Transpose[likeli.# & /@ prior];
        nFgl = Dimensions[sitefhaplo][[2]];        
        If[MatchQ[ToLowerCase[model], "depmodel"],
        	genotypes = Transpose[Table[Range[nFgl], {2}]],
        	genotypes = Flatten[Outer[List, Range[nFgl], Range[nFgl]], 1]
        ]; 
        ibd12 = Boole[Equal @@ # & /@ genotypes] + 1;
        derivedRule = {{"N", "N"} -> 1, {"N", "1"} -> 2, {"1", "N"} -> 
           3, {"N", "2"} -> 4, {"2", "N"} -> 5, {"1", "1"} -> 
           6, {"1", "2"} -> 7, {"2", "1"} -> 8, {"2", "2"} -> 9};
        pos = Replace[Partition[#, 2] & /@ sitefhaplo[[All, Flatten[genotypes]]],derivedRule, {2}];
        pos = Transpose[{ibd12, #}] & /@ pos;
        Outer[Extract[#2, #1] &, pos, likelidiplo, 1, 1]
    ]
    

               

(*likeli dimensions {noffspring, ntruehaplo=2}*)
(*prior dimesnions {ntruehaplo=2,nderivedhaplo=3}*)
(*likelihaplo dimensions {noffspring,nderivedhaplo=3}*)    
(*return dimensions {nparentalhaplo, noffspring,nfgl}*)  
siteMarginLikeliHaploGBS[sitefhaplo_, sitegeno_, epsF_, eps_,minphredscore_] :=
    Module[ {likeli, prior, likelihaplo, derhaplo},
        likeli = haploLikeliGBS[sitegeno, eps,minphredscore];
        prior = haploPrior[epsF];
        likelihaplo = likeli.prior;
        derhaplo = sitefhaplo /. {"N" -> 1, "1" -> 2, "2" -> 3};
        Outer[#2[[#1]] &, derhaplo, likelihaplo, 1, 1]
    ]
        
(*Set epsF=0 when sitefhaplo refers to the set of lattent haplotypes.*)
(*return dimensions {nparentalhaplo, noffspring,nstate},nstate=nfgl^2 for diploid, nfgl for haploid*) 
siteMagicLikelihoodGBS[model_,sitefhaplo_, sitegeno_, epsF_,eps_, minphredscore_,ismaleX_] :=
    Module[ {noffspring, maleX, nonmaleX, dataprob},
        noffspring = Length[sitegeno];
        maleX = Flatten[Position[ismaleX, True]];
        nonmaleX = Complement[Range[noffspring], maleX];
        dataprob = ConstantArray[0, {Length[sitefhaplo], noffspring}];
        If[ maleX =!= {},
            dataprob[[All, maleX]] = siteMarginLikeliHaploGBS[sitefhaplo, sitegeno[[maleX]],epsF, eps,minphredscore];
        ];
        If[ nonmaleX =!= {},
            dataprob[[All, nonmaleX]] = siteMarginLikeliDiploGBS[model,sitefhaplo, sitegeno[[nonmaleX]],epsF,eps,minphredscore];
        ];
        dataprob
    ]
     
    
End[] (* End Private Context *)

EndPackage[]