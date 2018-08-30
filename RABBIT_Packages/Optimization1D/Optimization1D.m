(* Mathematica Package *)

(* Created by the Wolfram Workbench Jul 6, 2016 *)

BeginPackage["Optimization1D`"]

brentLocalMin::usage = "brentLocalMin[f[#]&,lowbound,upbound] returs {xmin,f[xmin],updatehistory}. "

brentLocalMax::usage = "brentLocalMax  "

(* Exported symbols added here with SymbolName::usage *) 

Begin["`Private`"]
(* Implementation of the package *)

Options[brentLocalMin] = {
    AccuracyGoal -> 10, 
    PrecisionGoal -> 10,
    MaxIterations -> 100
    }

Options[brentLocalMax] = Options[brentLocalMin]

brentLocalMax[f_Function, start_?(NumberQ[#] || # == Automatic &),lowbound_?NumericQ, upbound_?NumericQ, 
  opts: OptionsPattern[]] :=
    Module[ {x, fx, his},
        {x, fx, his} = brentLocalMin[-f[#] &, start,lowbound, upbound, opts];
        fx *= -1;
        his[[All, -1]] *= -1;
        {x, fx, his}
    ]
    
brentLocalMin[f_Function, start_?(NumberQ[#] || # == Automatic &),lowbound_?NumericQ, upbound_?NumericQ, 
  OptionsPattern[]] :=
    Module[ {accuracygoal, precisiongoal, itmax, cgold,a, b, v, w, x=start, e, fv, fw,
       fx, iter, xm, tol1, tol2, p, q, r, d, u, fu,his},
        {accuracygoal, precisiongoal, itmax} = 
         OptionValue@{AccuracyGoal, PrecisionGoal, MaxIterations};
        cgold = 0.381966;(*(3-Sqrt[5])/2*)
        a = Min[lowbound, upbound];
        b = Max[lowbound, upbound];
        If[x===Automatic||(!(a<x<b)),
        	x=a + cgold (b - a)
        ];
        v = w =  a + cgold (b - a);
        e = 0;
        fv = fw = fx = f[x];
        his = ConstantArray[0,itmax];
        iter = 1;
        While[iter <= itmax,
         xm = 0.5 (a + b);
         tol1 = 10^(-precisiongoal) Abs[x] + 10^(-accuracygoal);
         tol2 = 2 tol1;
         (*Abs[x - xm] <= tol2 - 0.5 (b - a) <==> Max[x-a,b-x]<=tol2*)
         If[ Abs[x - xm] <= tol2 - 0.5 (b - a),
             his = Take[his,iter-1];
             Break[]
         ];
         p = q = r = 0.;
         If[ Abs[e] > tol1,
          (*Parabolic fit from x, v, and w*)
             r = (x - w) (fx - fv);
             q = (x - v) (fx - fw);
             p = (x - v) q - (x - w) r;
             q = 2 (q - r);
             (*result in proposal step size -p/q with always q\[GreaterEqual]0*)
             If[ q > 0,
                 p = -p,
                 q = -q
             ];
             (*r denotes the movement of the step before last*)
             r = e; 
             (*e denotes the movement of the last step*)
             e = d;
             ];
         If[ Abs[p] < Abs[0.5 q r] && p > q (a - x) && p < q (b - x),
             (*Parabolic interpolation step*)
             d = p/q;
             u = x + d;
             (*f must not be evalued too close to bracketing ends a or b*)
             If[ u - a < tol2 || b - u < tol2,
                 d = If[ x < xm,
                         tol1,
                         -tol1
                     ]
             ],
             (*Golden section step*)
             e = If[ x < xm,
                     b,
                     a
                 ] - x;
             d = cgold e;
         ];
         (*f must not be evalued too close to x*)
         u = x + If[ Abs[d] >= tol1,
                     d,
                     If[ d > 0,
                         tol1,
                         -tol1
                     ]
                 ];
         fu = f[u];
         (*housekeeping*)
         If[ fu <= fx,
             If[ u < x,
                 b = x,
                 a = x
             ];
             v = w;
             fv = fw;
             w = x;
             fw = fx;
             x = u;
             fx = fu,
             (**)
             If[ u < x,
                 a = u,
                 b = u
             ];
             If[ fu <= fw || w == x,
                 v = w;
                 fv = fw;
                 w = u;
                 fw = fu,
                 If[ fu <= fv || v == x || v == w,
                     v = u;
                     fv = fu
                 ];
             ];
         ];
         his[[iter]] = {iter,x,fx};
         iter++;
         ];
        If[ iter == itmax + 1,
            his = Take[his,iter-1];
            Print[Style[
              "brentLocalMin: exceed maximum iterations of " <> 
               ToString[itmax] <> "!", Red]]
        ];
        {x, fx, his}
    ]
  

End[]

EndPackage[]

