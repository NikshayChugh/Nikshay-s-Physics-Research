(* ::Package:: *)

BeginPackage["QSeries`"];

(* Global LaTeX output option *)
$QSeriesLaTeXOutput = False;

(* Global q variable *)
$q = q; (* Default value, can be changed *)

(* Declare symbols that should be in Global context *)
q; \[Eta]; \[Tau]; x; (* Add any other symbols you want in Global context *)

(* Ensure these symbols are in Global context *)
Global`q; Global`\[Eta]; Global`\[Tau]; Global`x;

(*Exported functions - now support both modes*)
aqprod::usage=
	"aqprod[a, n] or aqprod[a, q, n] computes the product (1 - a)(1 - aq)...(1 - aq^(n - 1)). Returns 1 if n = 0.";
etaq::usage = 
	"etaq[a, T] or etaq[q, a, T] returns the q-series expansion to O(q^T) of the eta-product.";
qbin::usage = 
	"qbin[m, n] or qbin[q, m, n] computes the q-binomial coefficient (n choose m) in terms of q.";
theta::usage = 
  "theta[z, T] or theta[z, q, T] computes the theta function, defined as the sum of z^i * q^(i^2) from i = -T to T.";
theta2::usage = 
  "theta2[T] or theta2[q, T] computes the theta function, defined as the sum of q^((i+1)^2/2) from i = -T to T.";
theta3::usage =
  "theta3[T] or theta3[q, T] computes the theta function, defined as the sum of q^(i^2) from i = -T to T.";
theta4::usage =
  "theta4[T] or theta4[q, T] computes the theta function, defined as the sum of (-1)^i * q^(i^2) from i = -T to T.";
jacprod::usage = 
  "jacprod[a, b, T] or jacprod[a, b, q, T] returns the q-series expansion to O(q^T) of the Jacobi-type infinite product.";
prodmake::usage = 
  "prodmake[f, T] or prodmake[f, q, T] converts the q-series f into a product expansion that agrees with f to O(q^T).";
qfactor::usage = 
  "qfactor[f, t] or qfactor[f, q, t] first factorises f and then applies prodmake to O(q^t).";
T::usage = 
 "T[r, j] computes Andrew's T function.";
etamake::usage = 
	"etamake[f, T] or etamake[f, q, T] converts the q-series f into an eta-product expansion that agrees with f to O(q^T).";
findPeriod::usage = 
  "findPeriod[A_list, T] finds the period of a list, if there exists some period."
jacprodmake::usage =
  "jacprobmake[f, T] or jacprodmake[f, q, T] converts the q-series f into a Jacobi product if one exists to order O(q^T)."
findcong::usage = 
	"findcong[QS_, T_, LM_:Null, XSET_:{}] computes a set of congruence relations for the q-series QS up to order T.";
findhom::usage = 
	"findhom[L_, n_, topshift_] or findhom[L_, q_, n_, topshift_] finds homogeneous relations of order n between the q-series in list L.";
checkprod::usage = "checkprod[f, M, Q] or checkprod[f, q, M, Q] checks if the q-series f is a 'nice' product up to O(q^Q), where 'nice' means the exponents in the product form have a max absolute value less than M.";
findprod::usage = "findprod[FL, T, M, Q] or findprod[FL, q, T, M, Q] searches for Z-linear combinations of the q-series in list FL that are probable products. T is the max coefficient in the linear combination, M is the max exponent size for a 'nice' product, and Q is the series order.";

(* LaTeX output control functions *)
SetLaTeXOutput::usage = "SetLaTeXOutput[True|False] enables or disables LaTeX output for all QSeries functions.";
QSLaTeX::usage = "QSLaTeX[expr] returns the LaTeX form of an expression.";

(* q variable control functions *)
SetQ::usage = "SetQ[value] sets the global q variable to the specified value.";
GetQ::usage = "GetQ[] returns the current value of the global q variable.";

Begin["`Private`"];

(* Helper function to format output *)
formatOutput[expr_] := If[$QSeriesLaTeXOutput, TeXForm[expr], expr];

(* LaTeX control functions *)
SetLaTeXOutput[bool_] := ($QSeriesLaTeXOutput = bool);
QSLaTeX[expr_] := TeXForm[expr];

(* q variable control functions *)
SetQ[value_] := ($q = value);
GetQ[] := $q;

(*aqprod - supports both 2 and 3 argument forms*)
aqprod[a_, n_] := aqprod[a, q, n];
aqprod[a_, q_, n_] := Module[{x = 1, i},
   If[IntegerQ[n] && n >= 0,
    Do[x = x*(1 - a*q^(i - 1)), {i, 1, n}],
    x = Subscript[{a, q}, n]];
   formatOutput[x]
];

(*etaq - supports both 2 and 3 argument forms*)
etaq[i, T] := etaq[q, i, T];
(* Optimized version that returns a SeriesData object directly *)
etaq[i_, T_] := etaq[q, i, T];
etaq[q_, i_, trunk_] := Module[{k, z1, z, w, terms},
  z1 = (1/6)*(i + Sqrt[i*i + 24*trunk*i])/i;
  z = 1 + Floor[N[z1]];
  
  (* Generate only the terms needed for the series *)
  terms = Table[
    w = (1/2)*i*k*(3*k - 1);
    If[w <= trunk, (-1)^k q^w, Nothing],
    {k, -z, z}
  ];
  
  (* Sum the terms and add the O(q) term to create the series *)
  Sum[t, {t, terms}] + O[q]^(trunk + 1)
];
(*qbin - supports both 2 and 3 argument forms*)
qbin[m_, n_] := qbin[q, m, n];
qbin[q_, m_, n_] := 
  If[IntegerQ[m] && IntegerQ[n],
    If[0 <= m <= n,
      output = Normal[aqprod[x, x, n]/(aqprod[x, x, m] * aqprod[x, x, n - m])] //FullSimplify //Expand;
      output = output /. x -> q // Expand;
      formatOutput[output],
      formatOutput[0]
    ],
    Message[qbin::argerr, m, n]
  ];

qbin::argerr = "m and n must be integers.";

(*theta functions - support both 1/2 and 2/3 argument forms*)
theta[z_, t_] := theta[z, q, t];
theta[z_, q_, t_] := formatOutput[Sum[z^i * q^(i^2), {i, -t, t}]];

theta2[t_] := theta2[q, t];
theta2[q_, t_] := formatOutput[theta[q,q,Floor[Sqrt[t]]+2]*q^(1/4)];

theta3[t_] := theta3[q, t];
theta3[q_, t_] := formatOutput[theta[1,q,Floor[Sqrt[t]]+1]];

theta4[t_] := theta4[q, t];
theta4[q_, t_] := formatOutput[theta[-1,q,Floor[Sqrt[t]]+1]];

(*jacprod - supports both 3 and 4 argument forms*)
jacprod[a_, b_, t_] := jacprod[a, b, q, t];
jacprod[a_, b_, q_, t_] := Module[{p1, p2, p3, maxIter},
  
  (* This corresponds to JAC(0,b,inf) *)
  If[a == 0, Return[etaq[q, b, t]]];

  (* Heuristic for iteration limit *)
  maxIter = Floor[t/b] + 2;

  (* Calculate each infinite product as a series *)
  p1 = Product[1 - q^a (q^b)^(i-1), {i, 1, maxIter}] + O[q]^(t+1);
  p2 = Product[1 - q^(b-a) (q^b)^(i-1), {i, 1, maxIter}] + O[q]^(t+1);
  p3 = etaq[q, b, t];
  
  (* Multiplying SeriesData objects is much faster *)
  formatOutput[p1 * p2 * p3]
];

(*prodmake - supports both 2/3 and 3/4 argument forms - combined best version*)
prodmake[f_, T_] := prodmake[f, q, T, False];
prodmake[f_, T_, returnList_] := prodmake[f, q, T, returnList];
prodmake[f_, q_, T_] := prodmake[f, q, T, False];
prodmake[f_, q_, T_, returnList_] := Module[
    {ft, f0, b, B, A, sum1, sum2, divj, divjb, prd, coeffs, result, n, m},
    
    (* Handle trivial cases early *)
    If[f === 1, Return[formatOutput[1]]];
    If[PolynomialQ[f, q] && Exponent[f, q] == 0, Return[formatOutput[f]]];
    
    (* Create series expansion *)
    ft = Series[f, {q, 0, T + 5}];
    
    If[Head[ft] =!= SeriesData,
        Message[prodmake::series, "f must be a series"];
        Return[$Failed]
    ];
    
    (* Check constant term *)
    f0 = SeriesCoefficient[ft, 0];
    
    If[f0 != 1,
        Message[prodmake::coeff, "Coefficient of q^0 must be 1"];
        Return[f]
    ];
    
    (* Extract coefficients with bounds checking *)
    coeffs = Table[SeriesCoefficient[ft, n], {n, 0, T}];
    
    If[Length[coeffs] < 2,
        Return[formatOutput[1]]
    ];
    
    (* Safely extract B coefficients *)
    B = coeffs[[2 ;; Min[T, Length[coeffs]]]];
    
    If[Length[B] == 0,
        Return[formatOutput[1]]
    ];
    
    (* Check if first coefficient is an integer *)
    If[!NumericQ[B[[1]]],
        Return[f]
    ];
    
    (* Initialize Association *)
    A = Association[1 -> B[[1]]];
    
    (* Main algorithm with robust bounds checking *)
    Do[
        If[n > Length[B], Break[]];
        
        sum2 = 0;
        Do[
            If[j < n && j <= Length[B],
                divj = Divisors[j];
                sum1 = Sum[If[KeyExistsQ[A, d], d * A[d], 0], {d, divj}];
                sum2 += B[[n - j]] * sum1
            ],
            {j, 1, n - 1}
        ];
        
        (* Use proper divisor calculation *)
        divjb = Select[Divisors[n], # < n &];
        sum1 = Sum[If[KeyExistsQ[A, d], d * A[d], 0], {d, divjb}];
        
        (* Calculate new coefficient *)
        sum2 = Expand[n * B[[n]] - sum2 - sum1];
        
        If[NumericQ[sum2/n],
            A[n] = sum2/n,
            Return[f] (* Return original if not integer *)
        ],
        {n, 2, Min[T - 1, Length[B]]}
    ];
    
    (* Build result *)
    If[!returnList,
        (* Build the product with proper bounds *)
        If[Length[Keys[A]] == 0,
            result = 1,
            prd = 1;
            Do[
                If[KeyExistsQ[A, m] && A[m] != 0,
                    prd *= (1 - q^m)^(-A[m])
                ],
                {m, 1, Min[T - 1, Max[Keys[A]]]}
            ];
            result = prd
        ],
        (* Return list of exponents *)
        result = Table[If[KeyExistsQ[A, m], -A[m], 0], {m, 1, T - 1}]
    ];
    
    formatOutput[result]
];

prodmake::args = "`1`";
prodmake::series = "`1`";
prodmake::coeff = "`1`";

(*qfactor - corrected implementation to match Maple version*)
qfactor[f_] := qfactor[f, Global`q];
qfactor[f_, q_] := Module[{d, T},
    (* Calculate default T value: 4d + 3 where d is max degree *)
    d = Max[Exponent[Numerator[f], q], Exponent[Denominator[f], q]];
    T = 4*d + 3;
    qfactor[f, q, T]
];

qfactor[f_, T_] := qfactor[f, Global`q, T];

qfactor[f_, q_, T_] := Module[{factored, result},
    (* First factor the function as usual *)
    factored = Factor[f];
    
    (* Handle different cases based on the type of result *)
    Which[
        (* If it's a symbol or rational number, return as is *)
        AtomQ[factored] || NumberQ[factored],
        formatOutput[factored],
        
        (* If it's a polynomial *)
        PolynomialQ[factored, q],
        bigqpfac[factored, q, T],
        
        (* If it's a rational function (fraction) *)
        !PolynomialQ[factored, q] && PolynomialQ[Numerator[factored], q] && PolynomialQ[Denominator[factored], q],
        result = bigqpfac[Numerator[factored], q, T] / bigqpfac[Denominator[factored], q, T];
        formatOutput[result],
        
        (* Otherwise, return an error or the original *)
        True,
        Message[qfactor::ratpoly, "f must be a rational polynomial"];
        formatOutput[f]
    ]
];

(* Helper function that does the actual q-factorization work *)
bigqpfac[f_, q_, T_] := Module[{fullFactors, prod, i, factor, power, constTerm, normalizedFactor, qProdResult},
    (* Factor the polynomial and apply prodmake to each factor *)
    fullFactors = FactorList[f];
    prod = 1;
    
    (* Process each factor *)
    Do[
        factor = fullFactors[[i, 1]];
        power = fullFactors[[i, 2]];
        
        Which[
            (* Handle constant factors *)
            NumberQ[factor],
            prod *= factor^power,
            
            (* Handle pure powers of q *)
            factor === q,
            prod *= q^power,
            
            (* Handle polynomial factors *)
            PolynomialQ[factor, q],
            (* Check if the factor has constant term 1 *)
            constTerm = factor /. q -> 0;
            
            If[constTerm == 1,
                (* Can apply prodmake directly *)
                qProdResult = prodmake[factor, q, T];
                If[qProdResult =!= factor, (* prodmake succeeded *)
                    prod *= qProdResult^power,
                    prod *= factor^power (* prodmake failed, keep original *)
                ],
                (* Factor has non-unit constant term, need to normalize *)
                If[constTerm != 0,
                    (* Normalize to make constant term 1 *)
                    normalizedFactor = Expand[factor/constTerm];
                    qProdResult = prodmake[normalizedFactor, q, T];
                    If[qProdResult =!= normalizedFactor, (* prodmake succeeded *)
                        prod *= constTerm^power * qProdResult^power,
                        prod *= factor^power (* prodmake failed, keep original *)
                    ],
                    (* Factor has no constant term, keep as is *)
                    prod *= factor^power
                ]
            ],
            
            (* For other cases, keep the factor as is *)
            True,
            prod *= factor^power
        ],
        {i, 1, Length[fullFactors]}
    ];
    
    prod
];

qfactor::ratpoly = "`1`";
(*Andrew's T function*)
ClearAll[T]
T[r_, j_] := T[r, j] = Module[{x},
    x = 0;
    If[j == 0 || j == 1,
        (j - 1)^2,
        Do[
            x = x - qbin[k, r + 2 k] * T[r + 2 k, j - 2 k],
            {k, 1, Floor[j/2]}
        ];
        Normal[Series[x, {q, 0, 1000}]]
    ]
];

(*etamake - supports both 2 and 3 argument forms*)
(*etamake - supports both 2 and 3 argument forms*)
etamake[f_, T_] := etamake[f, q, T];
etamake[f_, q_, T_] := Module[
  {
    qvar,                    (* q-variable for expansion *)
    seriesOrder,             (* Order of truncation *)
    seriesF,                 (* Truncated input series *)
    leadingCoeff,            (* Leading coefficient *)
    lowestDeg,               (* Lowest exponent of q *)
    g,                       (* Normalized series *)
    qExponent,               (* Accumulated q-exponent *)
    etaProduct,              (* Accumulated eta product *)
    correctionSeries,        (* g - 1, series to cancel *)
    correctionDegree,        (* Degree of first nonzero term *)
    correctionCoeff,         (* Coefficient at that degree *)
    maxOrder
  },

  (* Decide which q variable to use *)
  qvar = If[q === Automatic, q, q];

  (* Expand f to a high enough degree *)
  seriesOrder = T + 10;
  seriesF = Normal[Series[f, {qvar, 0, seriesOrder}]];

  (* Extract leading term info *)
  lowestDeg = Exponent[seriesF, qvar, Min];
  leadingCoeff = Coefficient[seriesF, qvar, lowestDeg];

  (* Normalize to make leading term q^0 with coeff 1 *)
  g = Expand[seriesF/(leadingCoeff * qvar^lowestDeg)];
  qExponent = lowestDeg;
  etaProduct = leadingCoeff;

  (* Final expansion length after removing q^lowestDeg *)
  maxOrder = T - lowestDeg;

  (* Begin the correction loop *)
  While[True,
    correctionSeries = Normal[Series[g - 1, {qvar, 0, maxOrder}]];

    (* Stop if fully corrected *)
    If[correctionSeries === 0, Break[]];

    (* Extract first nonzero term *)
    correctionDegree = Exponent[correctionSeries, qvar, Min];
    correctionCoeff = Coefficient[correctionSeries, qvar, correctionDegree];

    (* Update q exponent and eta product *)
    qExponent += (correctionDegree * correctionCoeff)/24;
    (* Use symbols from Global context *)
    etaProduct *= \[Eta][correctionDegree*\[Tau]]^(-correctionCoeff);

    (* Multiply normalized g by etaq^correctionCoeff *)
    g *= Expand[
      etaq[qvar, correctionDegree, maxOrder]^correctionCoeff
    ];
  ];

  (* Return eta product form *)
  qvar^qExponent * etaProduct
];

(* ::Section::, Closed:: *)
(*findPeriod*)
(* A more efficient, Mathematica-style findPeriod function *)
findPeriod[list_List, T_Integer] := Module[{p, sublist, partitions, check},

  sublist = list[[;; Min[T, Length[list]]]];

  For[p = 2, p <= Floor[Length[sublist]/2], p++,
    (* Partition the list into chunks of size p *)
    partitions = Partition[sublist[[;; Floor[Length[sublist]/p]*p]], p];

    (* Check if all chunks are identical *)
    check = Equal @@ partitions;

    If[check, Return[p]]; (* If they are, we found the period *)
  ];

  Return[$Failed]; (* If the loop finishes, no period was found *)
];

(* ::Section::, Closed:: *)
(*jacprodmake*)
jacprodmake[f_, T_] := jacprodmake[f, q, T];
jacprodmake[f_, q_, T_] := Module[
  {ft, ft0, newT, LT, seriesCoeffs, A, pp, jacCheckResult, y, mp, i},
  
  Print["Step 1: Starting jacprodmake..."];
  ft0 = Series[f, {q, 0, T + 5}];
  LT = Series[SeriesCoefficient[ft0, SeriesOrder[ft0]]*q^SeriesOrder[ft0], {q, 0, T + 5}];
  
  Print["Step 2: Calling prodmake to get exponents..."];
  seriesCoeffs = prodmake[f/Normal[LT], q, T, True];
  If[seriesCoeffs === $Failed || !ListQ[seriesCoeffs],
    Message[jacprodmake::noproduct];
    Return[f];
  ];
  Print["Step 3: Got exponents. Now finding the period..."];
  A = AssociationThread[Range[Length[seriesCoeffs]] -> -seriesCoeffs];
  
  (* This is the slow part *)
  pp = findPeriod[Values[A], T - 1];
  Print["Step 4: Period found: ", pp];
  
  If[pp === $Failed,
    Message[jacprodmake::noperiod];
    Return[f];
  ];
  
  Print["Step 5: Checking symmetry..."];
  jacCheckResult = And @@ Table[A[i] == A[pp - i], {i, 1, Floor[(pp - 1)/2]}];
  
  If[!jacCheckResult,
    Message[jacprodmake::nosymmetry];
    Return[f];
  ];

  Print["Step 6: Reconstructing product..."];
  y = 1;
  mp = Mod[pp, 2];
  If[mp == 0,
    Do[y *= (jacprod[i, pp, q, T]/jacprod[0, pp, q, T])^A[i], {i, 1, pp/2 - 1}];
    y *= (jacprod[pp/2, pp, q, T]/jacprod[0, pp, q, T])^(A[pp/2]/2);
  ,
    Do[y *= (jacprod[i, pp, q, T]/jacprod[0, pp, q, T])^A[i], {i, 1, (pp - 1)/2}];
  ];
  y *= jacprod[0, pp, q, T]^A[pp];
  
  Print["Step 7: Done."];
  formatOutput[Normal[LT]*y^(-1)]
];

jacprodmake::noproduct = "The series could not be converted to a product form.";
jacprodmake::noperiod = "No period found in the exponent sequence.";
jacprodmake::nosymmetry = "The exponent sequence does not have theta function symmetry.";

(* ::Section::, Closed:: *)
(*findcong*)
findcong[QS_, T_] := findcong[QS, T, Floor[Sqrt[T]], {}];
findcong[QS_, T_, LM_] := findcong[QS, T, LM, {}];
findcong[QS_, T_, LM_, XSET_] := Module[
  {S = {}, M, r, xxT, CFS, IGCD, primePowers, p, R, newCong, implied},
  
  For[M = 2, M <= LM, M++,
    For[r = 0, r < M, r++,
      xxT = Floor[(T - r)/M];
      If[xxT < 0, Continue[]];
      
      CFS = Table[SeriesCoefficient[QS, {q, 0, M*n + r}], {n, 0, xxT}];
      If[CFS === {0}, Continue[]];
      
      IGCD = GCD @@ CFS;
      
      If[Abs[IGCD] > 1,
        primePowers = FactorInteger[IGCD];
        For[i = 1, i <= Length[primePowers], i++,
          p = primePowers[[i, 1]];
          R = p^primePowers[[i, 2]];
          If[MemberQ[XSET, R], Continue[]];
          
          (* Check if this congruence is implied by an existing one *)
          implied = False;
          Do[
            If[p^primePowers[[i,2]] == S[[j,3]] && Mod[M, S[[j,2]]] == 0 && Mod[r, S[[j,2]]] == S[[j,1]],
              implied = True; Break[]
            ],
            {j, Length[S]}
          ];
          
          If[!implied,
            newCong = {r, M, R};
            Print[newCong];
            AppendTo[S, newCong];
          ]
        ]
      ]
    ]
  ];
  S
];


(* ::Section::, Closed:: *)
(*findhom*)
findhom[L_List, n_Integer, topshift_Integer] := findhom[L, q, n, topshift];
findhom[L_List, q_, n_Integer, topshift_Integer] := Module[
  {D = Length[L], X, monomials, numMons, order, coeffMatrix,
   nullVecs, relations, exponents, var, vec},
   
  monomials = MonomialList[Sum[X[i], {i, 1, D}]^n, Table[X[i], {i, 1, D}]];
  numMons = Length[monomials];
  order = numMons + 20 + topshift;
  
  coeffMatrix = Table[
    exponents = CoefficientRules[mon, Table[X[i], {i, 1, D}]][[1, 1]];
    
    (* Series expansion of the monomial *)
    series = Series[Product[L[[i]]^exponents[[i]], {i, 1, D}], {q, 0, order}];
    
    (* Row of coefficients *)
    Table[SeriesCoefficient[series, j], {j, 0, order - 1}]
    , {mon, monomials}
  ];
  
  nullVecs = NullSpace[coeffMatrix];
  
  relations = Table[
    vec . monomials,
    {vec, nullVecs}
  ];
  
  formatOutput[relations]
];

(* ::Section::, Closed:: *)
(*checkprod and findprod*)

checkprod[f_, M_Integer, Q_Integer] := checkprod[f, q, M, Q];
checkprod[f_, q_, M_Integer, Q_Integer] := Module[
  {ss, ssp, LL, nss, c0, nss1, xx, mx},
  
  If[f == 0, Return[{$Failed, Infinity}]];
  ss = Series[f, {q, 0, Q + 1}];
  LL = SeriesOrder[f];
  
  nss1 = Series[f / (SeriesCoefficient[f, LL] q^LL), {q, 0, Q - LL}];
  
  xx = prodmake[nss1, q, Q - LL, True];
  If[!ListQ[xx], Return[{$Failed, Infinity}]];
  
  mx = Max[Abs[xx]];
  If[mx < M, Return[{LL, 1}], Return[{LL, mx}]]
];

findprod[FL_List, T_Integer, M_Integer, Q_Integer] := findprod[FL, q, T, M, Q];
findprod[FL_List, q_, T_Integer, M_Integer, Q_Integer] := Module[
  {COMBS = {}, nfl = Length[FL], coeffVectors, vec, mcomb, cc},
  
  coeffVectors = Tuples[Range[-T, T], nfl];
  
  Do[
    vec = item;
    If[GCD @@ vec != 1, Continue[]];
    
    mcomb = vec . FL;
    cc = checkprod[mcomb, q, M, Q];
    
    If[cc[[2]] == 1,
      AppendTo[COMBS, {cc[[1]], vec}]
    ],
    {item, coeffVectors}
  ];
  
  COMBS
];

End[];
EndPackage[];