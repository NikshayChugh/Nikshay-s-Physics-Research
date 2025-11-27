(* ::Package:: *)
(*
 * Qseries.m
 *
 * A Mathematica translation of Frank Garvan's q-series Maple package.
 * Original Maple package: Version 1.1 (July 16, 2012)
 * Translated by: Nikshay Chugh
 * Date: September 5, 2025
 *
 * This package provides a comprehensive suite of tools for working with
 * q-series, including conversions to products, finding algebraic
 * relations, sifting coefficients, and implementing major product
 * identities. The logic and algorithms are a direct port of the original
 * Maple source code to preserve the behavior and methodology of
 * Dr. Garvan's work.
 *)

BeginPackage["Qseries`"];

(* Usage messages for public functions *)
aqprod::usage = "aqprod[a, q, n] computes the q-Pochhammer symbol (a;q)_n, defined as Product[(1 - a q^(i-1)), {i, 1, n}]. Handles n = Infinity.";
qbin::usage = "qbin[q, m, n] computes the q-binomial coefficient or Gaussian polynomial.";
etaq::usage = "etaq[q, k, T] computes the q-series expansion of the infinite product Product[(1 - q^(k*n)), {n, 1, Infinity}] up to O(q^(T+1)), using Euler's Pentagonal Number Theorem.";
theta::usage = "theta[z, q, T] computes the truncated theta series Sum[z^i q^(i^2), {i, -T, T}].";
theta2::usage = "theta2[q, T] computes the theta function \[Theta]2(q) up to O(q^(T+1)).";
theta3::usage = "theta3[q, T] computes the theta function \[Theta]3(q) up to O(q^(T+1)).";
theta4::usage = "theta4[q, T] computes the theta function \[Theta]4(q) up to O(q^(T+1)).";

prodmake::usage = "prodmake[f, q, T] converts the q-series f into an infinite product that agrees with f up to O(q^T). The leading coefficient of f must be 1.\nprodmake[f, q, T, \"List\"] returns the exponents of the q-product as a list.";
etamake::usage = "etamake[f, q, T] converts the q-series f into a product of eta functions, returning the result in terms of eta[k*\[Tau]].";
qetamake::usage = "qetamake[f, q, T] is a variant of etamake that returns the eta product in terms of a symbolic function QSeriesE[q^k].";
mprodmake::usage = "mprodmake[f, q, T] finds a product expansion of the form (1+q^n1)(1+q^n2)... for the series f.";
qfactor::usage = "qfactor[f, T] attempts to write a rational function f in q as a q-product of terms (1-q^i). T is an optional upper bound on the exponents.";

jacprod::usage = "jacprod[a, b, q, T] returns the q-series expansion of the Jacobi-type infinite product (q^a;q^b)\[Infinity] (q^(b-a);q^b)\[Infinity] up to O(q^(T+1)).";
jacprodmake::usage = "jacprodmake[f, q, T] attempts to convert the q-series f into a product of Jacobi-type theta functions (JAC symbols).";
jac2prod::usage = "jac2prod[jacexpr] converts a symbolic product of JAC functions into its q-product form using aqprod symbols.";
jac2series::usage = "jac2series[jacexpr, T] converts a symbolic product of JAC functions into a q-series expansion up to O(q^(T+1)).";

findhom::usage = "findhom[L, q, n, topshift] finds a set of potential homogeneous polynomial relations of degree n among the q-series in the list L.";
findhomcombo::usage = "findhomcombo[f, L, q, n, topshift, etaoption] tries to express the q-series f as a homogeneous polynomial of degree n in the members of L. If etaoption is \"yes\", the result is also given in terms of eta products.";
findnonhom::usage = "findnonhom[L, q, n, topshift] finds a set of potential non-homogeneous polynomial relations of degree up to n among the q-series in the list L.";
findnonhomcombo::usage = "findnonhomcombo[f, L, q, n, topshift] tries to express the q-series f as a non-homogeneous polynomial of degree up to n in the members of L.";
findpoly::usage = "findpoly[f1, f2, q, deg1, deg2, check] finds a polynomial relation P(X,Y)=0 between two q-series f1 and f2, with specified maximum degrees deg1 for X and deg2 for Y. Optionally checks the relation to O(q^check).";
findcong::usage = "findcong[QS, T, LM] finds linear congruences of the form c(A*n+B) == 0 mod C for the coefficients c(n) of a q-series QS.";

findhommodp::usage = "findhommodp[L, p, q, n, topshift] finds homogeneous relations mod p.";
findhomcombomodp::usage = "findhomcombomodp[f, L, p, q, n, topshift, etaoption] finds a homogeneous combination for f mod p.";
findlincombo::usage = "findlincombo[f, L, SL, q, topshift] expresses f as a linear combination of series in L, using names from SL.";
findlincombomodp::usage = "findlincombomodp[f, L, SL, p, q, topshift] expresses f as a linear combination of series in L mod p.";

sift::usage = "sift[s, q, m, r, T] takes a q-series s = Sum[a_n q^n] and returns the new series Sum[a_{m*n+r} q^n].";

tripleprod::usage = "tripleprod[z, q, T] computes Jacobi's triple product identity. If T is an integer, it returns the series expansion. If T is the string \"seriesid\", it returns the symbolic identity.";
quinprod::usage = "quinprod[z, q, T] computes the quintuple product identity. Use T as an integer for series expansion, or strings \"prodid\" or \"seriesid\" for symbolic forms.";
winquist::usage = "winquist[a, b, q, T] computes the series expansion of Winquist's identity up to O(q^(T+1)).";

packageversion::usage = "packageversion[] prints the current version information for the qseries package.";
changes::usage = "changes[] prints the version history and changes for the qseries package.";
xprint::usage = "xprint is a global flag (True/False) to enable verbose debugging output for certain functions.";

(* Public symbols that are used in output, like in Maple *)
q::usage = "q is the default variable for q-series generating functions.";
X::usage = "X is a symbolic variable used to represent q-series in polynomial relations found by functions like findhom and findpoly.";
Y::usage = "Y is a symbolic variable used to represent a second q-series in polynomial relations found by findpoly.";
JAC::usage = "JAC[a, b, Infinity] is a symbolic representation of a Jacobi-type theta product, used by jacprodmake.";
QSeriesE::usage = "QSeriesE[q^k] is a symbolic representation of an eta product, used by qetamake.";
η::usage = "η[k τ] is a symbolic representation of the Dedekind eta function used in the output of etamake.";
τ::usage = "τ is a symbolic variable used in the output of etamake.";


Begin["`Private`"];

(* ========================================================== *)
(* Section 0: Private Helpers & Globals                       *)
(* ========================================================== *)

(* Helper functions for series manipulation, mimicking Maple's behavior *)

degree[poly_, var_] := If[poly === 0, -Infinity, Exponent[poly, var]];
ldegree[poly_, var_] := If[poly === 0, Infinity, Exponent[poly, var, Min]];
tcoeff[poly_, var_] := If[poly === 0, 0, Coefficient[poly, var, ldegree[poly, var]]];
xprint = False; (* Global debug flag *)

(* ========================================================== *)
(* Section 1: Basic Functions                                 *)
(* ========================================================== *)

aqprod[a_, q_, n_Integer?NonNegative] := Product[1 - a*q^(i - 1), {i, 1, n}];
(* 1. Temporarily remove the protection *)
Unprotect[QPochhammer];

(* 2. Apply your custom formatting rule *)
Format[QPochhammer[a_, q_, n_:Infinity]] := Subscript[Row[{"(", a, ";", q, ")"}], n /. Infinity -> ∞];

(* 3. (Optional but good practice) Restore the protection *)
Protect[QPochhammer];

aqprod[a_, q_, Infinity] := QPochhammer[a, q];
aqprod[a_, q_, n_] := QPochhammer[a, q, n];

qbin[q_, m_Integer, n_Integer] := If[0 <= m <= n,
    aqprod[q, q, n] / (aqprod[q, q, m] * aqprod[q, q, n - m]) // Together,
    0
];

etaq[q_, i_, trunk_] := Module[{z},
    z = 1 + Floor[Re[(i + Sqrt[i*i + 24*trunk*i]) / (6.0*i)]];
    Sum[(-1)^k * q^(i*k*(3*k - 1) / 2), {k, -z, z}]
];

(* NOTE: The following theta definitions are as per user request and may appear non-standard. *)
theta[z_, q_, t_Integer] := Sum[z^i * q^(i^2), {i, -t, t}];

theta2[q_, t_Integer] := Module[{n},
    n = Floor[Sqrt[t]] + 2;
    q^(1/4) * theta[q, q, n]
];

theta3[q_, t_Integer] := Module[{n},
    n = Floor[Sqrt[t]] + 1;
    theta[1, q, n]
];

theta4[q_, t_Integer] := Module[{n},
    n = Floor[Sqrt[t]] + 1;
    theta[-1, q, n]
];

(* ========================================================== *)
(* Section 2: Product Conversion                              *)
(* ========================================================== *)

prodmake[f_, T_] := prodmake[f, q, T, "Product"];
prodmake[f_, T_, output_] := prodmake[f, q, T, output];
prodmake[f_, q_, T_] := prodmake[f, q, T, "Product"];

prodmake[f_, q_, T_, output_String] := Module[
    {ft, f0, bCoeffs, aCoeffs, sum1, sum2, result, n, j, d, prd, m},
    
    (* Handle trivial cases *)
    If[f === 1, Return[1]];
    
    ft = Series[f, {q, 0, T + 5}];
    If[Head[ft] =!= SeriesData, Message[prodmake::series, f]; Return[$Failed]];
    
    f0 = SeriesCoefficient[ft, 0];
    If[f0 != 1, Message[prodmake::coeff, f0]; Return[f]];
    
    (* Extract coefficients of the input series f *)
    bCoeffs = Table[SeriesCoefficient[ft, i], {i, 1, T -1}];
    If[Length[bCoeffs] == 0, Return[1]];

    (* Use Association for sparse storage of calculated exponents *)
    aCoeffs = Association[];
    If[Length[bCoeffs] >= 1, aCoeffs[1] = bCoeffs[[1]]];
    
    Do[
        sum2 = 0;
        Do[
            sum1 = Sum[d * Lookup[aCoeffs, d, 0], {d, Divisors[j]}];
            sum2 += bCoeffs[[n - j]] * sum1,
            {j, 1, n - 1}
        ];
        
        sum1 = Sum[d * Lookup[aCoeffs, d, 0], {d, Most[Divisors[n]]}];
        aCoeffs[n] = (n * bCoeffs[[n]] - sum2 - sum1) / n // Expand;
    ,
        {n, 2, Length[bCoeffs]}
    ];
    
    If[output === "List",
        result = Table[-Lookup[aCoeffs, m, 0], {m, 1, T - 1}],
    (* else *)
        prd = 1;
        Do[
            prd *= (1 - q^m)^(-Lookup[aCoeffs, m, 0]),
            {m, Keys[aCoeffs]}
        ];
        result = prd
    ];
    
    result
];

prodmake::series = "Input `1` could not be expanded as a series.";
prodmake::coeff = "Constant term of series must be 1, but was `1`.";


etamake[f_, T_] := etamake[f, q, T];
etamake[f_, q_, T_] := Module[
  {
    seriesOrder, seriesF, leadingCoeff, lowestDeg, g,
    qExponent, etaProduct, correctionSeries, correctionDegree,
    correctionCoeff, maxOrder
  },
  
  seriesOrder = T + 10;
  seriesF = Normal[Series[f, {q, 0, seriesOrder}]];

  lowestDeg = Exponent[seriesF, q, Min];
  leadingCoeff = Coefficient[seriesF, q, lowestDeg];

  g = Expand[seriesF/(leadingCoeff * q^lowestDeg)];
  qExponent = lowestDeg;
  etaProduct = leadingCoeff;

  maxOrder = T - lowestDeg;

  (* Begin the correction loop *)
  While[True,
    correctionSeries = Normal[Series[g - 1, {q, 0, maxOrder}]];

    (* Stop if fully corrected *)
    If[correctionSeries === 0, Break[]];

    (* Extract first nonzero term *)
    correctionDegree = Exponent[correctionSeries, q, Min];
    correctionCoeff = Coefficient[correctionSeries, q, correctionDegree];

    (* Update q exponent and eta product *)
    qExponent += (correctionDegree * correctionCoeff)/24;
    etaProduct *= \[Eta][correctionDegree*\[Tau]]^(-correctionCoeff);

    (* Update g by multiplying by the correction factor *)
    g *= Expand[etaq[q, correctionDegree, maxOrder]^correctionCoeff];
  ];

  (* Return the final eta product form *)
  q^qExponent * etaProduct
];

(* Convenience definition that defaults to the variable 'q' *)
qetamake[f_, T_] := qetamake[f, q, T];

(* Main function with robust loop and clear variables *)
qetamake[f_, q_, T_] := Module[
  {
    seriesOrder, seriesF, leadingCoeff, lowestDeg, g,
    qExponent, symbolicProduct, correctionSeries, correctionDegree,
    correctionCoeff, maxOrder
  },
  
  seriesOrder = T + 10;
  seriesF = Normal[Series[f, {q, 0, seriesOrder}]];

  lowestDeg = Exponent[seriesF, q, Min];
  leadingCoeff = Coefficient[seriesF, q, lowestDeg];

  g = Expand[seriesF/(leadingCoeff * q^lowestDeg)];
  qExponent = lowestDeg;
  symbolicProduct = leadingCoeff;

  maxOrder = T - lowestDeg;

  (* Robust correction loop *)
  While[True,
    correctionSeries = Normal[Series[g - 1, {q, 0, maxOrder}]];

    (* Exit when series is fully corrected *)
    If[correctionSeries === 0, Break[]];

    (* Find the first term to cancel *)
    correctionDegree = Exponent[correctionSeries, q, Min];
    correctionCoeff = Coefficient[correctionSeries, q, correctionDegree];

    (* Update the symbolic product using QSeriesE *)
    symbolicProduct *= QSeriesE[q^correctionDegree]^(-correctionCoeff);

    (* Update g by multiplying by the correction factor *)
    g *= Expand[etaq[q, correctionDegree, maxOrder]^correctionCoeff];
  ];

  (* Return the final product form *)
  q^qExponent * symbolicProduct
];

(* Convenience definition that defaults to the variable 'q' *)
mprodmake[f_, T_] := mprodmake[f, q, T];

(* Main function with robust loop and clear variables *)
mprodmake[f_, q_, T_] := Module[
  {
    seriesOrder, seriesF, leadingCoeff, lowestDeg, g,
    qExponent, resultProduct, correctionSeries, correctionDegree,
    correctionCoeff, maxOrder
  },
  
  seriesOrder = T + 10;
  seriesF = Normal[Series[f, {q, 0, seriesOrder}]];

  lowestDeg = Exponent[seriesF, q, Min];
  leadingCoeff = Coefficient[seriesF, q, lowestDeg];

  (* Normalize the series to have a leading term of 1 *)
  g = Expand[seriesF/(leadingCoeff * q^lowestDeg)];
  qExponent = lowestDeg;
  resultProduct = leadingCoeff;

  maxOrder = T - lowestDeg;

  (* Robust correction loop *)
  While[True,
    correctionSeries = Normal[Series[g - 1, {q, 0, maxOrder}]];

    (* Exit when series is fully corrected *)
    If[correctionSeries === 0, Break[]];

    (* Find the first term that needs to be cancelled *)
    correctionDegree = Exponent[correctionSeries, q, Min];
    correctionCoeff = Coefficient[correctionSeries, q, correctionDegree];

    (* Update the result with a (1+q^k) term *)
    resultProduct *= (1 + q^correctionDegree)^correctionCoeff;

    (* Update the series by dividing out the correction factor *)
    g *= Expand[(1 + q^correctionDegree)^(-correctionCoeff)];
  ];

  (* Return the final product form *)
  q^qExponent * resultProduct
];

(* Private helper function that does the actual q-factorization work *)
bigqpfac[f_, q_, T_] := Module[{fullFactors, prod = 1, i, factor, power, constTerm, normalizedFactor, qProdResult},
    (* Factor the polynomial and apply prodmake to each factor *)
    fullFactors = FactorList[f];
    
    (* Process each factor from the FactorList *)
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
            constTerm = factor /. q -> 0;
            
            If[constTerm == 1,
                (* If constant term is 1, apply prodmake directly *)
                qProdResult = prodmake[factor, q, T];
                prod *= qProdResult^power,
            (* else *)
                (* If constant term is not 1, normalize the factor first *)
                If[constTerm != 0,
                    normalizedFactor = Expand[factor/constTerm];
                    qProdResult = prodmake[normalizedFactor, q, T];
                    prod *= constTerm^power * qProdResult^power,
                (* else *)
                    (* Factor has no constant term (e.g., q^2+q), keep as is *)
                    prod *= factor^power
                ]
            ],
            
            (* For other symbolic cases, keep the factor as is *)
            True,
            prod *= factor^power
        ],
        {i, 1, Length[fullFactors]}
    ];
    
    prod
];

(* Public qfactor function with multiple definitions for convenience *)
qfactor[f_, T_: Automatic] := qfactor[f,q,T];
qfactor[f_, q_, T_: Automatic] := Module[{d, defaultT},
    d = Max[Exponent[Numerator[f], q], Exponent[Denominator[f], q]];
    defaultT = If[T === Automatic, 4*d + 3, T];
    
    Module[{factored},
        factored = Factor[f];
        Which[
            AtomQ[factored] || NumberQ[factored],
            factored,
            
            PolynomialQ[factored, q],
            bigqpfac[factored, q, defaultT],
            
            True, (* Assumes rational function *)
            bigqpfac[Numerator[factored], q, defaultT] / bigqpfac[Denominator[factored], q, defaultT]
        ]
    ]
];

qfactor::ratpoly = "Input `1` must be a rational polynomial in q.";

(* ========================================================== *)
(* Section 3: Jacobi Theta Function Tools                     *)
(* ========================================================== *)

(* Main jac2prod function now takes the variable 'var' as an argument *)

jac2prod[jacExpression_, var_] := jacExpression /. JAC[a_, b_, Infinity] :> Private`jac[a, b, var, Infinity];

(* Convenience version defaults to the global symbol 'q' *)

jac2prod[jacExpression_] := jac2prod[jacExpression, q];

(* Main jac2series function now takes the variable 'var' as an argument *)

jac2series[jacExpression_, var_, T_] := jacExpression //. {JAC[a_, b_, Infinity] :> Private`qjac[a, b, var, T], Infinity -> T};
(* Convenience version defaults to the global symbol 'q' *)
jac2series[jacExpression_, T_] := jac2series[jacExpression, q, T];


(* --- Private Helper Functions (now context-aware) --- *)

(* Private`jac now takes 'var' to build the product *)
Private`jac[a_Integer, b_Integer, var_, Infinity] /; a >= 0 && a < b := If[a > 0,
    aqprod[var^a, var^b, Infinity] * aqprod[var^(b-a), var^b, Infinity] * aqprod[var^b, var^b, Infinity],
    aqprod[var^b, var^b, Infinity]
];

(* Private`tripleprod2 now takes 'var' for the series *)
Private`tripleprod2[z_, var_, t_] := Sum[(-1)^i z^i var^(i(i-1)/2), {i, -t, t}];

(* Private`qjac now takes 'var' and passes it to its helpers *)
Private`qjac[a_Integer, b_Integer, var_, t_Integer] := If[a =!= 0,
    Module[{c},
        c = Floor[Re[1/2*(b-2*a+Sqrt[4*a^2-4*a*b+b^2+8*b*t])/b]] + 1;
        Private`tripleprod2[var^a, var^b, c]
    ],
    etaq[var, b, t] (* Note: etaq already accepts the variable as its first argument *)
];

(* Main discovery function: finds the JAC product for a given q-series *)

jacprodmake[seriesInput_, q_, T_Integer] := Module[{
    seriesF, seriesPoly, lowestDeg, leadingCoeff, leadingTerm, normalizedPoly, 
    exponents, jacExpr
    },
    
    (* 1. Get the series expansion *)
    seriesF = Series[seriesInput, {q, 0, T + 5}];
    If[!SeriesDataQ[seriesF], Message[prodmake::"noseries", seriesInput]; Return[$Failed]];
    
    (* 2. Correctly find the lowest degree and leading coefficient *)
    seriesPoly = Normal[seriesF];
    If[seriesPoly === 0, Return[0]];
    lowestDeg = Exponent[seriesPoly, q, Min];
    leadingCoeff = Coefficient[seriesPoly, q, lowestDeg];
    leadingTerm = leadingCoeff * q^lowestDeg;
    
    (* 3. Normalize the series to have a leading term of 1 *)
    normalizedPoly = Normal[Series[seriesInput / leadingTerm, {q, 0, T - 1}]];
    
    (* 4. Use the prodmake algorithm to find the exponents *)
    exponents = Private`calculateExponents[normalizedPoly, q, T];
    
    (* 5. Pass exponents to jacmake to find the theta-product structure *)
    jacExpr = Private`jacmake[exponents, T - 1];
    
    (* 6. Re-assemble the final result *)
    If[jacExpr === $Failed,
        Return[$Failed],
        Return[leadingTerm * jacExpr]
    ]
];

(* --- Private Helper Functions --- *)

(* Helper to calculate the exponent sequence using the prodmake recurrence *)
Private`calculateExponents[poly_, q_, T_] := Module[{bCoeffs, aCoeffs = Association[], n, j, sum1, sum2, d},
    bCoeffs = Table[Coefficient[poly, q, i], {i, 1, T - 1}];
    If[Length[bCoeffs] > 0, aCoeffs[1] = bCoeffs[[1]]];
    Do[
        sum2 = Sum[bCoeffs[[n-j]] * Sum[d*Lookup[aCoeffs,d,0], {d,Divisors[j]}], {j,1,n-1}];
        sum1 = Sum[d*Lookup[aCoeffs,d,0], {d,Most[Divisors[n]]}];
        aCoeffs[n] = (n*bCoeffs[[n]] - sum2 - sum1) / n;,
        {n, 2, Length[bCoeffs]}
    ];
    Table[Lookup[aCoeffs, i, 0], {i, 1, T - 1}]
];

(* Helper to find the period of the exponent list *)
Private`periodfind[exponentList_List, maxOrder_Integer] := Module[{len = Length[exponentList], period, isPeriodFlag, i, j, lastj},
  For[period = 2, period <= Floor[len/2], period++,
    isPeriodFlag = True;
    lastj = Floor[len/period] - 2;
    If[lastj < 0, Continue[]];
    For[j = 0, j <= lastj, j++,
      For[i = 1, i <= period, i++,
        If[exponentList[[i + period*j]] != exponentList[[i + period*j + period]],
          isPeriodFlag = False; Break[];];
      ];
      If[!isPeriodFlag, Break[]];
    ];
    If[isPeriodFlag, Return[period]];
  ];
  Return[Floor[maxOrder/2] + 2]; (* Indicates no period found *)
];

(* Helper to check for theta-function symmetry in the exponents *)
Private`jaccheck[exponentList_List, period_Integer] := Module[{n = Floor[period/2], isSymmetric = 1, i},
    For[i = 1, i <= n, i++,
        If[exponentList[[i]] =!= exponentList[[period - i]], isSymmetric = 0; Break[];]
    ];
    Return[isSymmetric];
];

(* Helper to reconstruct the symbolic JAC product from the exponents and period *)
Private`jacmake[exponentList_List, maxOrder_Integer] := Module[{period, symbolicProduct = 1, i},
    period = Private`periodfind[exponentList, maxOrder];
    If[period > maxOrder/2, Print["No period found."]; Return[$Failed]];
    If[Private`jaccheck[exponentList, period] =!= 1, Print["Not a theta function."]; Return[$Failed]];

    symbolicProduct = (JAC[0, period, Infinity])^(-exponentList[[period]]);
    If[EvenQ[period],
        symbolicProduct *= Product[(JAC[i, period, Infinity]/JAC[0, period, Infinity])^(-exponentList[[i]]), {i, 1, period/2 - 1}];
        symbolicProduct *= (JAC[period/2, period, Infinity]/JAC[0, period, Infinity])^(-exponentList[[period/2]]/2);
    ,
        symbolicProduct *= Product[(JAC[i, period, Infinity]/JAC[0, period, Infinity])^(-exponentList[[i]]), {i, 1, (period - 1)/2}];
    ];
    Return[symbolicProduct // PowerExpand // Together];
];

(* ========================================================== *)
(* Section 4: Finding Algebraic Relations                     *)
(* ========================================================== *)

findhom[seriesList_List, q_, degree_Integer, topshift_Integer: 0] := Module[
  {
    numSeries = Length[seriesList],
    polyVars, exponentVectors, monomialBasis, 
    numMonomials, numCoefficients, seriesOrder,
    expandedSeries, coeffMatrix,
    nullSpaceVectors, simplifiedVectors, relationsList
  },
  
  (* 1. Generate the monomial basis from scratch, avoiding MonomialList *)
  polyVars = Array[X, numSeries];
  
  (* A) Find all sets of exponents {e1, e2,...} that sum to 'degree' *)
  exponentVectors = FrobeniusSolve[Table[1, numSeries], degree];
  
  (* B) Build the monomials from the exponents, e.g., X[1]^e1 * X[2]^e2 ... *)
  monomialBasis = Inner[Power, polyVars, #, Times] & /@ exponentVectors;
  
  
  (* --- The rest of the code is now guaranteed to work --- *)
  
  numMonomials = Length[monomialBasis];
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients;

  expandedSeries = Series[#, {q, 0, seriesOrder}] & /@ seriesList;

  coeffMatrix = Table[
    powers = Exponent[monomial, #] & /@ polyVars;
    monomialSeries = Product[expandedSeries[[i]]^powers[[i]], {i, 1, numSeries}];
    PadRight[CoefficientList[Normal[monomialSeries], q], numCoefficients],
    {monomial, monomialBasis}
  ];

  nullSpaceVectors = NullSpace[Transpose[coeffMatrix]];
  simplifiedVectors = If[nullSpaceVectors === {}, {}, LatticeReduce[nullSpaceVectors]];
  
  relationsList = If[simplifiedVectors === {}, {}, simplifiedVectors . monomialBasis];
  
  Print["Different homogeneous relations of degree = ", degree, " among the given series:"];
  Return[relationsList];
];

findhomcombo[targetFunction_, basisFunctions_List, q_, degree_Integer, topshift_Integer: 0, etaoption_String: "no"] := Module[
  {
    numBasisFuncs = Length[basisFunctions],
    polyVars, exponentVectors, monomialBasis,
    numMonomials, numCoefficients, seriesOrder,
    expandedBasis, expandedTarget,
    monomialSeriesList, coeffMatrix,
    kernel, relations = {}, etaRelations = {}
  },
  
  (* 1. Build the monomial basis from scratch to bypass the corrupted MonomialList *)
  polyVars = Array[X, numBasisFuncs];
  exponentVectors = FrobeniusSolve[Table[1, numBasisFuncs], degree];
  monomialBasis = Inner[Power, polyVars, #, Times] & /@ exponentVectors;
  numMonomials = Length[monomialBasis];
  
  (* 2. Calculate how many series terms are needed *)
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];

  (* 3. Pre-compute all series expansions ONCE for efficiency *)
  expandedBasis = Series[#, {q, 0, seriesOrder}] & /@ basisFunctions;
  expandedTarget = Series[targetFunction, {q, 0, seriesOrder}];
  
  (* 4. Build the series for each monomial combination *)
  monomialSeriesList = Table[
    Inner[Power, expandedBasis, exponents, Times],
    {exponents, exponentVectors}
  ];
  
  (* 5. Build the full coefficient matrix *)
  coeffMatrix = Join[
     (* Rows for the monomial combinations *)
     Table[PadRight[CoefficientList[Normal[s], q], numCoefficients], {s, monomialSeriesList}],
     (* Final row for the target function *)
     {PadRight[CoefficientList[Normal[expandedTarget], q], numCoefficients]}
  ];
  
  (* 6. Find the linear dependency *)
  kernel = NullSpace[Transpose[coeffMatrix]];
  If[Length[kernel] > 1, Print["WARNING: dim ker = ", Length[kernel], ". Solution may not be unique."]];
  
  (* 7. Reconstruct the polynomial expression(s) *)
  Do[
    With[{vec = k, lastCoeff = Last[k]},
      If[lastCoeff =!= 0,
        (* Reconstructs P(l_1, l_2, ...) from the relation P(l_i) + c*f = 0 *)
        Dim[polyExpression = Drop[vec, -1] . monomialBasis];
        AppendTo[relations, -polyExpression / lastCoeff];
        
        (* Handle the optional eta-function transformation *)
        If[etaoption === "yes",
          etaExpression = Drop[vec, -1] . (etamake[Normal[#], q, seriesOrder - 2] & /@ monomialSeriesList);
          AppendTo[etaRelations, -etaExpression / lastCoeff];
        ]
      ]
    ],
    {k, kernel}
  ];
  
  Print["Different homogeneous combinations of degree = ", degree, " between the given series and the target function:"];
  If[etaoption === "yes", Print[etaRelations]];
  Return[relations];
];

findnonhom[functionList_List, q_, maxDegree_Integer, topshift_Integer: 0] := Module[
  {
    numFuncs = Length[functionList],
    polyVars, exponentVectors, allMonomials = {1}, (* Start with the constant term *)
    numMonomials, numCoefficients, seriesOrder,
    expandedFunctions, coeffMatrix,
    kernel, relations
  },
  
  (* 1. Build the basis from scratch, bypassing MonomialList *)
  polyVars = Array[X, numFuncs];
  
  (* A) Loop through each degree from 1 to maxDegree *)
  Do[
    exponentVectors = FrobeniusSolve[Table[1, numFuncs], currentDegree];
    allMonomials = Join[allMonomials, Inner[Power, polyVars, #, Times] & /@ exponentVectors],
    {currentDegree, 1, maxDegree}
  ];
  
  numMonomials = Length[allMonomials];
  
  (* 2. Calculate series order and pre-compute all series *)
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];
  expandedFunctions = Series[#, {q, 0, seriesOrder}] & /@ functionList;
  
  (* 3. Build the coefficient matrix efficiently *)
  coeffMatrix = Table[
      (* Get the series for the current monomial *)
      monomialSeries = Switch[monomial,
        1, Series[1, {q, 0, seriesOrder}], (* Handle the constant term *)
        _, Module[{powers = Exponent[monomial, polyVars]},
             Inner[Power, expandedFunctions, powers, Times]
           ]
      ];
      PadRight[CoefficientList[Normal[monomialSeries], q], numCoefficients],
      {monomial, allMonomials}
  ];
  
  (* 4. Find the null space and reconstruct the relations *)
  kernel = NullSpace[Transpose[coeffMatrix]];
  relations = kernel . allMonomials;
  
  Print["Different non-homogeneous relations up to degree = ", maxDegree, " among the given series:"];
  Return[relations];
];

findnonhomcombo[targetFunction_, basisFunctions_List, q_, maxDegree_Integer, topshift_Integer: 0, etaoption_String: "no"] := Module[
  {
    numBasisFuncs = Length[basisFunctions],
    polyVars, exponentVectors, monomialBasis = {1}, (* Start with degree 0 (constant) *)
    numMonomials, numCoefficients, seriesOrder,
    expandedBasis, expandedTarget,
    monomialSeriesList, coeffMatrix,
    kernel, relations = {}, etaRelations = {}
  },
  
  (* 1. Build basis from scratch, generating monomials for each degree up to maxDegree *)
  polyVars = Array[X, numBasisFuncs];
  Do[
    exponentVectors = FrobeniusSolve[Table[1, numBasisFuncs], currentDegree];
    monomialBasis = Join[monomialBasis, Inner[Power, polyVars, #, Times] & /@ exponentVectors],
    {currentDegree, 1, maxDegree}
  ];
  numMonomials = Length[monomialBasis];
  
  (* 2. Calculate series order and pre-compute all series *)
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];
  expandedBasis = Series[#, {q, 0, seriesOrder}] & /@ basisFunctions;
  expandedTarget = Series[targetFunction, {q, 0, seriesOrder}];

  (* 3. Build the series for each monomial combination *)
  monomialSeriesList = Table[
    Switch[monomial,
      1, Series[1, {q, 0, seriesOrder}],
      _, Module[{powers = Exponent[monomial, polyVars]},
           Inner[Power, expandedBasis, powers, Times]
         ]
    ],
    {monomial, monomialBasis}
  ];

  (* 4. Build the full coefficient matrix *)
  coeffMatrix = Join[
     Table[PadRight[CoefficientList[Normal[s], q], numCoefficients], {s, monomialSeriesList}],
     {PadRight[CoefficientList[Normal[expandedTarget], q], numCoefficients]}
  ];
  
  (* 5. Find the linear dependency and reconstruct the expression(s) *)
  kernel = NullSpace[Transpose[coeffMatrix]];
  If[Length[kernel] > 1, Print["WARNING: dim ker = ", Length[kernel], ". Solution may not be unique."]];
  
  Do[
    With[{vec = k, lastCoeff = Last[k]},
      If[lastCoeff =!= 0,
        polyExpression = Drop[vec, -1] . monomialBasis;
        AppendTo[relations, -polyExpression / lastCoeff];
        
        If[etaoption === "yes",
          etaExpression = Drop[vec, -1] . (etamake[Normal[#], q, seriesOrder - 2] & /@ monomialSeriesList);
          AppendTo[etaRelations, -etaExpression / lastCoeff];
        ]
      ]
    ],
    {k, kernel}
  ];

  Print["Possible non-homogeneous expressions up to degree = ", maxDegree, " between the given series and the target function:"];
  If[etaoption === "yes", Print[etaRelations]];
  Return[relations];
];

findpoly[x_, y_, q_, deg1_Integer, deg2_Integer, check_Integer: Automatic] := Module[{
    dim1, dim2, a, num, k, j, b, qq, l, ta, kk, poly, i, polyg, polyfunc, ss
    },
    (* NOTE: This function's reliance on global X and Y is a design choice from the original package. *)
    Print["WARNING: X, Y are global."];
    dim1 = (deg1 + 1) * (deg2 + 1);
    dim2 = dim1 + 10;
    Print[" dims ", dim1, " ", dim2];
    
    b = Flatten[Table[X^k Y^j, {k, 0, deg1}, {j, 0, deg2}]];
    num = Length[b];
    a = ConstantArray[0, {num, dim2}];

    For[i = 1, i <= num, i++,
        qq = Series[b[[i]] /. {X->x, Y->y}, {q, 0, dim2 + 2}];
        For[l = 0, l < dim2, l++, a[[i, l + 1]] = Coefficient[qq, q, l]];
    ];
    
    ta = Transpose[a];
    kk = NullSpace[ta];
    
    If[kk === {}, Print[" NO polynomial relation found. "]; Return[$Failed]];
    
    poly = b . kk[[1]];
    polyg = Sum[Factor[Coefficient[poly, Y, deg2 - j]] * Y^(deg2 - j), {j, 0, deg2}];
    
    Print["The polynomial is"];
    Print[polyg];
    
    If[IntegerQ[check],
        polyfunc = Function[{X, Y}, Evaluate[polyg]];
        ss = Series[polyfunc[x, y], {q, 0, check}];
        Print["Checking to order ", check];
        Print[ss];
    ];
    
    Return[polyg];
];

findcong[qSeries_, maxPower_Integer, modLimit_:Automatic, excludedDivisors_:{}] := Module[
  {
    (* === Setup === *)
    mMax,               
    allCoeffs,          
    foundCongruences = {}, 
    
    (* === Loop & Intermediate Variables === *)
    coeffsInProgression, 
    sampleGCD,          
    fullGCD,            
    primeFactors,       
    isNew,              
    newCongruence       
  },
  
  (* 1. Set the modulus limit and pre-calculate all series coefficients for efficiency. *)
  mMax = If[modLimit === Automatic, Floor[Sqrt[maxPower]], modLimit];
  allCoeffs = CoefficientList[Normal[Series[qSeries, {q, 0, maxPower}]], q];
  
  (* 2. Main Loop: Iterate through each modulus 'm' and remainder 'r'. *)
  Do[
    coeffsInProgression = allCoeffs[[r + 1 ;; ;; m]];
    If[Length[coeffsInProgression] < 2, Continue[]];
    
    (* 3. Quick Check: Calculate GCD on a small sample first to fail fast. *)
    sampleGCD = GCD @@ Take[coeffsInProgression, Min[Length[coeffsInProgression], 50]];
    If[sampleGCD == 1, Continue[]];
    
    (* 4. Full Check: If the sample check passes, calculate the GCD for all coefficients. *)
    fullGCD = GCD @@ coeffsInProgression;
    
    (* 5. Process any non-trivial GCD found. IGNORE IF GCD IS 0 (THE FIX IS HERE). *)
    If[fullGCD != 0 && fullGCD =!= 1 && !MemberQ[excludedDivisors, fullGCD],
      primeFactors = FactorInteger[fullGCD];
      
      Scan[
        (newCongruence = {r, m, First[#]^Last[#]};
        
        isNew = !AnyTrue[foundCongruences,
          ({oldR, oldM, oldPpow} = #;
           oldPpow == newCongruence[[3]] &&       
           Mod[m, oldM] == 0 &&                  
           Mod[r, oldM] == oldR                  
          ) &
        ];
        
        If[isNew,
          foundCongruences = Select[foundCongruences,
            ({oldR, oldM, oldPpow} = #;
             !(oldPpow == newCongruence[[3]] &&    
               Mod[oldM, m] == 0 &&               
               Mod[oldR, m] == r                  
             )
            ) &
          ];
          
          AppendTo[foundCongruences, newCongruence];
        ];
        
        ) &,
        primeFactors
      ];
    ];
    
  , {m, 2, mMax}, {r, 0, m - 1}];
  
  Print["Found congruences {remainder, modulus, prime_power}:"];
  Return[Sort[foundCongruences]];
];


findhommodp[seriesList_List, p_Integer, q_, degree_Integer, topshift_Integer: 0] := Module[
  {
    numSeries = Length[seriesList],
    polyVars, exponentVectors, monomialBasis,
    numMonomials, numCoefficients, seriesOrder,
    expandedSeries, coeffMatrix,
    nullSpaceVectors, relationsList
  },
  
  (* 1. Build the monomial basis robustly *)
  polyVars = Array[X, numSeries];
  exponentVectors = FrobeniusSolve[Table[1, numSeries], degree];
  monomialBasis = Inner[Power, polyVars, #, Times] & /@ exponentVectors;
  
  (* 2. Determine series expansion order and pre-compute expansions *)
  numMonomials = Length[monomialBasis];
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];
  
  expandedSeries = Series[#, {q, 0, seriesOrder}] & /@ seriesList;

  (* 3. Build the coefficient matrix efficiently *)
  coeffMatrix = Table[
    (* For each monomial, compute the corresponding q-series product *)
    monomialSeries = Inner[Power, expandedSeries, exponents, Times];
    PadRight[CoefficientList[Normal[monomialSeries], q], numCoefficients],
    {exponents, exponentVectors}
  ];

  (* 4. Find the null space over the finite field Z_p *)
  nullSpaceVectors = NullSpace[Transpose[coeffMatrix], Modulus -> p];
  
  (* 5. Reconstruct the polynomial relations from the null space vectors *)
  relationsList = Mod[nullSpaceVectors . monomialBasis, p];
  
  Print["Homogeneous relations of degree ", degree, " mod ", p, ":"];
  Return[relationsList];
];

findhomcombomodp[f_, l_List, p_Integer, q_, n_Integer, topshift_Integer: 0, etaoption_String: "no"] := Module[
  {
    (* Variable setup *)
    targetFunction = f, basisFunctions = l, degree = n,
    numBasisFuncs = Length[basisFunctions],
    polyVars, exponentVectors, monomialBasis,
    numMonomials, numCoefficients, seriesOrder,
    
    (* Series and Matrix variables *)
    expandedBasis, expandedTarget,
    monomialSeriesList, coeffMatrix, kernel,
    
    (* Result variables *)
    relations = {}, etaRelations = {}
  },
  
  (* 1. Build the monomial basis robustly, avoiding MonomialList *)
  polyVars = Array[X, numBasisFuncs];
  exponentVectors = FrobeniusSolve[Table[1, numBasisFuncs], degree];
  monomialBasis = Inner[Power, polyVars, #, Times] & /@ exponentVectors;
  numMonomials = Length[monomialBasis];
  
  (* 2. Calculate how many series coefficients are needed *)
  numCoefficients = numMonomials + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];

  (* 3. Pre-compute all series expansions ONCE for efficiency *)
  expandedBasis = Series[#, {q, 0, seriesOrder}] & /@ basisFunctions;
  expandedTarget = Series[targetFunction, {q, 0, seriesOrder}];
  
  (* 4. Build the series for each monomial combination from the pre-computed parts *)
  monomialSeriesList = Table[
    Inner[Power, expandedBasis, exponents, Times],
    {exponents, exponentVectors}
  ];
  
  (* 5. Construct the full coefficient matrix to solve the linear system *)
  coeffMatrix = Join[
     (* Rows for the monomial combinations *)
     Table[PadRight[CoefficientList[Normal[s], q], numCoefficients], {s, monomialSeriesList}],
     (* Final row for the target function *)
     {PadRight[CoefficientList[Normal[expandedTarget], q], numCoefficients]}
  ];
  
  (* 6. Find the linear dependency over the finite field Z_p *)
  kernel = NullSpace[Transpose[coeffMatrix], Modulus -> p];
  If[Length[kernel] > 1, Print["WARNING: dim ker = ", Length[kernel], ". Solution may not be unique."]];
  
  (* 7. Reconstruct the polynomial expression(s) from the null space vectors *)
  relations = Cases[kernel, k_List /; Last[k] =!= 0 :>
    Module[{inv = PowerMod[Last[k], -1, p]},
      -(Drop[k, -1] * inv) . monomialBasis
    ]
  ];

  (* Handle the optional eta-function transformation *)
  If[etaoption === "yes",
    etaRelations = Cases[kernel, k_List /; Last[k] =!= 0 :>
      Module[{inv = PowerMod[Last[k], -1, p]},
        -(Drop[k, -1] * inv) . (etamake[Normal[#], q, seriesOrder - 2] & /@ monomialSeriesList)
      ]
    ];
  ];
  
  (* 8. Print and return the final results *)
  Print["Homogeneous combinations of degree ", degree, " mod ", p, ":"];
  If[etaoption === "yes", Print[Mod[etaRelations, p]]];
  Return[Mod[relations, p]];
];

findlincombo[targetSeries_, basisFunctions_List, basisSymbols_List, q_, topshift_Integer: 0] := Module[
  {
    numBasisFuncs = Length[basisFunctions],
    numCoefficients, seriesOrder,
    allSeries, coeffMatrix, kernel, relations = {}
  },
  
  (* 1. Determine series expansion order *)
  numCoefficients = numBasisFuncs + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];
  
  (* 2. Combine all series and build the coefficient matrix efficiently *)
  allSeries = Join[basisFunctions, {targetSeries}];
  coeffMatrix = Table[
    PadRight[CoefficientList[Normal[Series[s, {q, 0, seriesOrder}]], q], numCoefficients],
    {s, allSeries}
  ];
  
  (* 3. Find the null space over the rational numbers *)
  kernel = NullSpace[Transpose[coeffMatrix]];
  If[Length[kernel] == 0, 
    Print["No linear combination found."];
    Return[{}];
  ];
  If[Length[kernel] > 1, Print["WARNING: dim ker = ", Length[kernel], ". Solution is not unique."]];
  
  (* 4. Reconstruct the symbolic expression(s) *)
  Do[
    With[{vec = k, lastCoeff = Last[k]},
      If[lastCoeff =!= 0,
        symbolicCombination = -(Drop[vec, -1] / lastCoeff) . basisSymbols;
        AppendTo[relations, symbolicCombination];
      ]
    ],
    {k, kernel}
  ];
  
  Print["Found linear combinations:"];
  Return[relations];
];

findlincombomodp[targetSeries_, basisFunctions_List, basisSymbols_List, p_Integer, q_, topshift_Integer: 0] := Module[
  {
    numBasisFuncs = Length[basisFunctions],
    numCoefficients, seriesOrder,
    allSeries, coeffMatrix, kernel, relations = {}
  },
  
  (* 1. Determine series expansion order *)
  numCoefficients = numBasisFuncs + 20 + topshift;
  seriesOrder = numCoefficients + 1;
  Print["Number of series terms to check: ", seriesOrder];
  
  (* 2. Combine all series and build the coefficient matrix efficiently *)
  allSeries = Join[basisFunctions, {targetSeries}];
  coeffMatrix = Table[
    PadRight[CoefficientList[Normal[Series[s, {q, 0, seriesOrder}]], q], numCoefficients],
    {s, allSeries}
  ];
  
  (* 3. Find the null space over the finite field Z_p *)
  kernel = NullSpace[Transpose[coeffMatrix], Modulus -> p];
  If[Length[kernel] == 0, 
    Print["No linear combination found mod ", p, "."];
    Return[{}];
  ];
  If[Length[kernel] > 1, Print["WARNING: dim ker = ", Length[kernel], ". Solution is not unique."]];
  
  (* 4. Reconstruct the symbolic expression(s) *)
  Do[
    With[{vec = k, lastCoeff = Last[k]},
      If[lastCoeff =!= 0,
        inv = PowerMod[lastCoeff, -1, p];
        symbolicCombination = (-Drop[vec, -1] * inv) . basisSymbols;
        AppendTo[relations, symbolicCombination];
      ]
    ],
    {k, kernel}
  ];
  
  Print["Found linear combinations mod ", p, ":"];
  Return[Mod[relations, p]];
];


(* ========================================================== *)
(* Section 5: Sifting Coefficients                            *)
(* ========================================================== *)

sift[s_, q_, n_, k_, T_] := Module[{y = 0, i, st, lasti},
    st = Series[s, {q, 0, T + 4}];
    If[!SeriesDataQ[st], Message[sift::"noseries", s]; Return[$Failed]];
    lasti = Floor[(T - k) / n];
    y = Sum[Coefficient[st, q, n*i + k] * q^i, {i, 0, lasti}];
    Return[y];
];
sift::"noseries" = "Input `1` must be a valid series.";


(* ========================================================== *)
(* Section 6: Product Identities                              *)
(* ========================================================== *)

tripleprod[z_, q_, T_Integer?Positive] := Module[{x = 0, lasti, i},
    lasti = Floor[Sqrt[2*T + 1/4] + 1/2] + 1;
    x = Sum[(-1)^i * z^i * q^(i*(i - 1)/2), {i, -lasti, lasti}];
    Return[x];
];
tripleprod[z_, q_,n_,m_, "seriesid"] := HoldForm[
    Product[(1 - z*q^(n-1))(1 - q^n/z)(1 - q^n), {n, 1, Infinity}] == Sum[(-1)^m z^m q^(m(m-1)/2), {m, -Infinity, Infinity}]
];

(*
  This function calculates a truncated series approximation of the quintuple product identity.
  It corresponds to the case where the third argument, T, is a positive integer.

  Arguments:
  z: A complex or real variable.
  q: A complex or real variable, typically with |q| < 1.
  T: A positive integer representing the truncation limit for the power of q.
*)
quinprod[z_, q_, T_Integer?Positive] := Module[{lasti},
    (*
      Determine the summation limit 'lasti'. We need to find the largest integer 'i'
      such that the exponent of q, which is i(3i+1)/2, does not exceed T.
      This is found by solving 3i^2 + i - 2T <= 0 for its positive root.
    *)
    lasti = Floor[(-1 + Sqrt[1 + 24*T])/6] + 1;
    
    (*
      Compute the sum from -lasti to lasti. This is the truncated version of the
      series side of the quintuple product identity.
    *)
    Sum[ ((-z)^(-3*i) - (-z)^(3*i + 1)) * q^(i*(3*i + 1)/2), {i, -lasti, lasti}]
];


(*
  This function returns the formal statement of Watson's Quintuple Product Identity,
  equating the infinite product with its infinite series representation.
  It corresponds to the case where T is the string "seriesid".
  
  HoldForm prevents Mathematica from attempting to evaluate the expressions.
*)
quinprod[z_, q_,m_, "seriesid"] := With[{
    leftid = aqprod[-z, q, Infinity]*aqprod[-q/z, q, Infinity]*aqprod[z^2*q, q^2, Infinity]*aqprod[q/z^2, q^2, Infinity]*aqprod[q, q, Infinity]
    },
    HoldForm[
        leftid == Sum[((-z)^(-3*m) - (-z)^(3*m + 1))*q^(m*(3*m + 1)/2), {m, -Infinity, Infinity}]
    ]
];

(*
  This function returns an alternate form of the identity, equating the main product
  with a sum of two other products.
  It corresponds to the case where T is the string "prodid".
*)
quinprod[z_, q_, "prodid"] := With[{
    leftid = aqprod[-z, q, Infinity] * aqprod[-q/z, q, Infinity] * aqprod[z^2*q, q^2, Infinity] * aqprod[q/z^2, q^2, Infinity] * aqprod[q, q, Infinity],
    x1 = aqprod[q^2/z^3, q^3, Infinity] * aqprod[q*z^3, q^3, Infinity] * aqprod[q^3, q^3, Infinity],
    x2 = aqprod[q/z^3, q^3, Infinity] * aqprod[q^2*z^3, q^3, Infinity] * aqprod[q^3, q^3, Infinity]
    },
    (* Return the identity in HoldForm to prevent evaluation. *)
    HoldForm[leftid == x1 + z*x2]
];

winquist[a_, b_, q_, T_Integer?Positive] := Module[{lasti},
    lasti = Floor[7/6 + 1/6 * Sqrt[25 + 24*T]] + 1;
    Sum[(-1)^(i + j) * (
        (a^(-3*i) - a^(3*i + 3)) * (b^(-3*j) - b^(3*j + 1)) +
        (a^(-3*j + 1) - a^(3*j + 2)) * (b^(3*i + 2) - b^(-3*i - 1))
        ) * q^(3*i*(i + 1)/2 + j*(3*j + 1)/2),
        {i, 0, lasti}, {j, -lasti, lasti}
    ]
];

(* ========================================================== *)
(* Section 7: Utility and Versioning                          *)
(* ========================================================== *)

packageversion[] := (
    Print["**************************************************************"];
    Print["*"];
    Print["* qseries package (Mathematica Version)"];
    Print["* Based on Maple Version 1.1 - Mon Jul 16 15:24:26 EDT 2012"];
    Print["* This version tested on Mathematica 14+."];
    Print["*"];
    Print["* Please report any problems or discrepancies."];
    Print["* Original Author: Frank Garvan (fgarvan@ufl.edu)"];
    Print["* Homepage: http://www.math.ufl.edu/~fgarvan"];
    Print["*"];
    Print["**************************************************************"];
);

changes[] := (
    Print["This is a Mathematica translation of Frank Garvan's qseries package."];
    Print["The changelog below refers to the original Maple versions."];
    Print["**************************************************************"];
);

(* ========================================================== *)
(* Section 8: Comments and Known Issues                       *)
(* ========================================================== *)
(*
 * jacprod: The usage message for jacprod[a, b, q, T] is defined,
 * but the function itself is not implemented in this version.
 * The private helper Private`qjac is used by jac2series.
 *
 * findpoly: This function should be replaced with findPolySimple as in the jupyter nb but that one does not work correctly for our package. 
 * findcong: The algorithm was corrected to handle cases where the
 * GCD of an arithmetic progression of coefficients is 0, which previously
 * caused errors. It now correctly ignores these progressions. I am not yet sure of the functioning of this because Example 2 fails. 
 *
*)

(* ======================================================== *)
(* Section 9: *)

End[]; (* `Private` context *)

EndPackage[];