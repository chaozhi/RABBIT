(* Mathematica Package *)

(* Created by the Wolfram Workbench Jan 22, 2017 *)

BeginPackage["SpectralOrdering`"]
(* Exported symbols added here with SymbolName::usage *) 

spectralClustering::usage = "spectralClustering  "

spectralOrdering::usage = "spectralOrdering  "

spetralEigenPlot::usage = "spetralEigenPlot  "

toKNNSimilarity::usage = "toKNNSimilarity  "

estimateKNN::usage = "estimateKNN  "

graphLaplacian::usage = "graphLaplacian is an option to specify one of three graph Laplacians: \"unNormalized\",\"rwNormalized\",or \"symNormalized\"."

toRWTransition::usage = "toRWTransition  "

toSymSimilarity::usage = "toSymSimilarity  "

eigenVectorSelection::usage = "eigenVectorSelection is an option to specify how to select eigen vectors for clustering. It can take \"fixed\",\"eigenval\", \"eigendiff\", and \"eigenratio\""

Begin["`Private`"]
(* Implementation of the package *)

Options[spectralClustering] = Options[spectralOrdering] = Options[spetralEigenPlot] = 
	{graphLaplacian -> "rwNormalized",
	 eigenVectorSelection ->"eigenratio"
	}
        
getEigenval[laplacian_, n_] :=
    Module[ {opt},
        opt = If[ Length[laplacian]>=5&&n<Length[laplacian],
                  Method -> "Arnoldi",
                  Method -> Automatic
              ];
        Eigenvalues[laplacian, -n, opt]
    ]        

spetralEigenPlot[similarity_, n_,isprint_:False] :=
    Module[ {laptype,laplacian,ev,figev},
        laptype = OptionValue[graphLaplacian];
        Switch[ToLowerCase[laptype],
            "unnormalized",
            laplacian = DiagonalMatrix[Total[similarity, {2}]] - similarity,
            "rwnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toRWTransition[similarity],
            "symnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toSymSimilarity[similarity],
            _,
            Print["spectralOrdering: wrong option value of graphLaplacian ",laptype,"!"];
            Abort[];
        ];
        ev = Reverse[getEigenval[laplacian, Min[2 n,Length[similarity]]]];
        figev = spetralEigenPlot[ev,n,isprint];
        {ev,figev}
    ]
    
spetralEigenPlot[eigenvalues_,ncluster_,evsel_,isprint_:False] :=
    Module[ {k,n,ev},
    	If[Depth[eigenvalues]==3,
    		k=Length[First[eigenvalues]];
    		ev=Flatten[eigenvalues],
	        ev = Union[Flatten[eigenvalues]];	        
	        Switch[evsel,
	        	"fixed",k = First[splitEigenval[ev,ncluster,isprint]],
	        	"eigenval",k = First[splitEigenval2[ev,ncluster,isprint]],
	        	"eigendiff",k = First[splitEigenval3[ev,ncluster,isprint]],
	        	"eigenratio",k = First[splitEigenval4[ev,ncluster,isprint]],
	        	_, 
	        	Print["Unknown eigenvector selection: ",evsel];Abort[]
	        ];
    	];
    	n=Length[ev];
    	ev = Transpose[{Range[n], ev}];
    	If[ k < n,
	            ev = {ev[[;; k]], ev[[k + 1 ;;]]}
	        ];
        ListPlot[ev, PlotMarkers -> {Automatic,12}, PlotRange -> All, 
          Filling -> Axis, Frame -> {{True, False}, {True, False}}, 
          FrameLabel -> {"i-th smallest", "Eigenvalue"}, 
          LabelStyle -> Directive[FontSize -> 13, FontFamily -> "Helvetica", Black],
          PlotStyle -> {Blue, Gray}, 
          PlotLabel -> ToString[k]<>" eigenvectors (blue) are used for clustering."]
    ]
                  
toRWTransition[aa_] :=
    Module[ {idd},
        idd = 1/Total[aa, {2}];
        aa idd
    ]

toSymSimilarity[aa_] :=
    Module[ {idd},
        idd = Total[aa, {2}]^(-1/2);
        Transpose[Transpose[aa] idd] idd
    ]       


spectralClustering[similarity_?MatrixQ, ncluster_Integer?(# > 1 &), OptionsPattern[]] :=
    Module[ {laptype,evsel,laplacian, eigenval0,eigenval, eigenvec,kls, k, res,nc,order, cluster, neigen,opt},
        (*laplacian = DiagonalMatrix[Total[similarity, {2}]] - similarity;*)
        laptype = OptionValue[graphLaplacian];
        evsel=OptionValue[eigenVectorSelection];
        Switch[ToLowerCase[laptype],
            "unnormalized",
            laplacian = DiagonalMatrix[Total[similarity, {2}]] - similarity,
            "rwnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toRWTransition[similarity],
            "symnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toSymSimilarity[similarity],
            _,
            Print["spectralClustering: wrong option of graphLaplacian ",laptype,"!"];
            Abort[];
        ];
        If[evsel=="fixed",
	        neigen = If[ Length[similarity]>ncluster+1,
	                     ncluster+1,
	                     Length[similarity]
	                 ],
	       neigen = If[ Length[similarity]>2 ncluster,
	                     2 ncluster,
	                     Length[similarity]
	                 ];
        ];
        opt = If[ Length[similarity]>=5&& 2 neigen<Length[similarity],
                  Method -> "Arnoldi",
                  Method -> Automatic
              ];
        {eigenval0, eigenvec} = Eigensystem[laplacian, -neigen, opt];
        eigenval0 = Reverse[eigenval0];
        Switch[evsel,
        	"fixed",kls = splitEigenval[eigenval0,ncluster],
        	"eigenval",kls = splitEigenval2[eigenval0,ncluster],
        	"eigendiff",kls = splitEigenval3[eigenval0,ncluster],
        	"eigenratio",kls = splitEigenval4[eigenval0,ncluster],
        	_, 
        	Print["Unknown eigenvector selection: ",evsel];Abort[]
        ];
        order = Ordering[eigenvec[[-2]]];
        k=kls[[1]];
        cluster = FindClusters[Transpose[eigenvec[[-k ;;, order]]] -> order,ncluster, 
            	DistanceFunction -> CosineDistance, Method -> {"Agglomerate",ClusterDissimilarityFunction -> "Average"}];
        If[Length[kls]>1&&Length[cluster]<ncluster,
	        res={{Length[cluster],k,cluster}};
	        Print["#selected eigenvec = ",k, "; #cluster = ", Length[cluster]];
	        Do[cluster = FindClusters[Transpose[eigenvec[[-k ;;, order]]] -> order,k, 
            		DistanceFunction -> CosineDistance, Method -> {"Agglomerate",ClusterDissimilarityFunction -> "Average"}];
               nc = Length[cluster];
	           Print["#selected eigenvec = ",k, "; #cluster = ", nc];
	           If[nc>ncluster,Break[]];
	           res=Join[res,{{nc,k,cluster}}],{k,kls}];
	        res=SortBy[res,{#[[1]],-#[[2]]}&];
	        {k,cluster}=res[[-1,2;;3]];
        ];    
        eigenval={eigenval0[[;; k]], eigenval0[[k + 1 ;;]]};
        {Sort[Sort[#] & /@ cluster, Length[#1] > Length[#2] &], eigenval, eigenvec}
    ]  

splitEigenval[eigenval_,ncluster_,isprint_:False] :={ncluster}
    
splitEigenval2[eigenval_,ncluster_,isprint_:False] :=
    Module[ {k,ev},
        (*The smallest eigval should be 0, but it is numerically close to 0*)
        (*Clustering with the smallest eigval is sensitive to shifting and scaling of similarity matrix*)
        ev = SortBy[FindClusters[Rest[eigenval]], First];
        ev[[1]] = Prepend[ev[[1]], First[eigenval]];
        If[isprint,Print["clustering of eigenvalues: ", ev]];
        k = Accumulate[Length[#] & /@ ev];
        Select[k, # >= ncluster &]
    ]
    
splitEigenval3[eigenval_,ncluster_,isprint_:False] :=
    Module[ {diff,peak},
    	diff=Differences[eigenval[[ncluster ;;]]];
        peak = Append[PeakDetect[diff], 1];
        peak=Flatten[Position[PeakDetect[peak], 1]];
        If[isprint, Print["differences between eigenvalues[[", ncluster,";;]]: ",diff,"; peak=",peak]];
        If[peak==={},{ncluster},ncluster + peak - 1]     
    ]    

splitEigenval4[eigenval_, ncluster_, isprint_: False] := 
 Module[{ratio, peak},
  ratio = Abs[eigenval[[ncluster ;;]]]+10^(-10.);
  ratio = ratio[[2 ;;]]/ratio[[;; -2]];  
  peak=Flatten[Position[PeakDetect[ratio], 1]];
  If[isprint, Print["Ratios between eigenvalues[[", ncluster,";;]]: ",ratio,"; peak=",peak]];
  If[peak==={},{ncluster},ncluster + peak - 1]
  ] 
  
(*splitEigenval4[eigenval_, ncluster_, isprint_: False] := 
 Module[{ratio, k, ev},
  ratio = Abs[eigenval[[ncluster ;;]]]+10^(-10.);
  ratio = ratio[[2 ;;]]/ratio[[;; -2]];
  If[isprint, Print["ratios between eigenvalues[[", ncluster,";;]]: ",ratio]];
  k = ncluster + First[Ordering[ratio, -1]] - 1;
  ev = {eigenval[[;; k]], eigenval[[k + 1 ;;]]};
  {k, ev}
  ]*)            
        
spectralOrdering[similarity_?MatrixQ, OptionsPattern[]] :=
    Module[ {laptype,laplacian,opt,n = 2,k,order,ev,fiedler},
        laptype = OptionValue[graphLaplacian];
        Switch[ToLowerCase[laptype],
            "unnormalized",
            laplacian = DiagonalMatrix[Total[similarity, {2}]] - similarity,
            "rwnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toRWTransition[similarity],
            "symnormalized",
            laplacian = IdentityMatrix[Length[similarity]] - toSymSimilarity[similarity],
            _,
            Print["spectralOrdering: wrong option value of graphLaplacian ",laptype,"!"];
            Abort[];
        ];
        opt = If[ Length[similarity]>=5,
                  Method -> "Arnoldi",
                  Method -> Automatic
              ];
        While[True,
          ev = Eigensystem[laplacian, -n, opt];
          k = Count[ev[[1]], _?(# <= 10^(-10.) &)];          
          If[ k <n,
              fiedler = If[ ToLowerCase[laptype]=="symnormalized",
                            (Total[similarity, {2}]^(-1/2))  ev[[2, -(k + 1)]],
                            ev[[2, -(k + 1)]]
                        ];
              order = Ordering[fiedler];
              Break[],
              If[ 2 n < Length[laplacian],
                  n *= 2,
                  Abort[]
              ]
          ];
        ];
        order
    ]     
      
spectralOrdering[similarity_?MatrixQ, cluster_List, opts:OptionsPattern[]] :=
    Module[ {ncluster,indices,order,i},
        If[ Complement[Union[Flatten[cluster]],Range[Length[similarity]]]=!={},
            Print["spectralOrdering: wrong input clusters!"];
            Abort[]
        ];
        ncluster = Length[cluster];
        If[ ncluster==1,
            {spectralOrdering[similarity,opts]},
            Table[
                indices = cluster[[i]];
                If[ Length[indices]>1,
                    PrintTemporary[Style["Starting ordering for cluster " <> ToString[i] <> 
                       " of size " <> ToString[Length[indices]] <> "", Blue]];
                    order = spectralOrdering[similarity[[indices, indices]],opts];
                    indices[[order]],
                    indices
                ],{i,Length[cluster]}]
        ]
    ]
                
adjConnectedComponents[similarity_] :=
    ConnectedComponents[AdjacencyGraph[Abs[Sign[similarity]]]]
 
estimateKNN[similarity_,knnmin_:2] :=
    Module[ {delt = 0,history, ncmax, knnmax, ncmin, knn, nc,aa,count = 0},
        knnmax = Length[similarity];
        ncmin = Length[adjConnectedComponents[similarity]];
        aa = toKNNSimilarity[similarity, knnmin,delt];
        ncmax = Length[adjConnectedComponents[aa]];
        history = {{knnmin, ncmax}, {knnmax, ncmin}};
        (*If[isprint,PrintTemporary["initial history = ", history]];*)
        If[ ncmax == ncmin,
            {knnmin, history},
            While[True,
                 knn = 2 history[[-2, 1]];
                 If[ knn >= history[[-1, 1]],
                     knn = Floor[Mean[history[[-2 ;;, 1]]]]
                 ];
                 aa = toKNNSimilarity[similarity, knn,delt];
                 nc = Length[adjConnectedComponents[aa]];
                 count++;
                 (*If[isprint,PrintTemporary["iteration = ",count,", knn = ", knn, ", nc = ", nc]];*)
                 If[ nc == history[[-1, 2]],
                     PrependTo[history, history[[-1]]];
                     history[[-1]] = {knn, nc},
                     history = Insert[history, {knn, nc}, -2];
                 ];
                 If[ history[[-1, 1]] == history[[-2, 1]] + 1,
                     Break[]
                 ];
             ];
            {history[[-1, 1]], SortBy[history, First]}
        ]
    ] 

toKNNSimilarity[similarity_, knn_String,delt_:0] :=
    Module[ {knn2, knnhis = {},figknn},
        Switch[knn,
            "Auto",
            {knn2, knnhis} = estimateKNN[similarity];
            figknn = ListLogLogPlot[knnhis, PlotRange->All,ImageSize -> 350,
              PlotMarkers -> {Automatic, 10}, 
              Filling -> Axis, Frame -> {{True, False}, {True, False}}, 
              FrameLabel -> {"# nearest neighbors (k)", "#connected components"}, 
              PlotLabel -> "Keep k = " <> ToString[knn2] <> " nearest neighbors"];
            Print[figknn],
            "Sqrt",
            knn2 = Round[Sqrt[Length[similarity]]],
            "All",
            knn2 = Length[similarity],
            _, 
            Print["toKNNSimilarity: wrong "<>knn];
            Abort[]
        ];
        toKNNSimilarity[similarity, knn2,delt]
    ]
    
toKNNSimilarity[similarity_, knn_Integer,delt_:0] :=
    Module[ {aa, max, knn2, res, pos, i},
        If[ !(2<=knn<=Length[similarity]),
            Print["kNNsimilarity: wrong number of nearest neighbor (knn)!"];
            Abort[]
        ];
        If[ knn==Length[similarity],
            SparseArray[similarity],
            aa = N[Normal[similarity]];
            max = N[Max[aa]];
            knn2 = (Count[Normal[#], 0.] & /@ Quiet[aa - max]) /. {_?(# <= knn &) -> knn};
            (*knn2 = Table[knn,{Length[aa]}];*)
            res = delt Sign[aa];
            Do[
             pos = Ordering[aa[[i]], -knn2[[i]]];
             res[[i, pos]] = aa[[i, pos]], {i, Length[res]}];
            res = MapThread[Max, {res, Transpose[res]}, 2];
            SparseArray[res, Dimensions[similarity], delt]
        ]
    ] 
          

End[]

EndPackage[]

