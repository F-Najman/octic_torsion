/* Derived verbatim from nt-lib/CurveArith commit0263b3435d34c5efe6a0fade37059cac317b7295, src/gonality.m. GPL-3.0; see experimental/CurveArith/LICENSE. Derickx--Terao arXiv2602.17417, Algorithm3.3. */
function UniformizingDifferential(omegas, x)
/*
Given a sequence omegas of differentials on a function field and a place x on the same function field,
returns a differential of omegas with valuation 0 at x.
*/
    for omega in omegas do
        if Valuation(Divisor(omega), x) eq 0 then
            return omega;
        end if;
    end for;
    error "No differential with valuation zero";
end function;

function UniformizingParameter(fs, x)
/*
Given a sequence fs of functions on a function field and a place x on the same function field,
returns an element of the function field with valuation 1 at x.

TODO: can remove fallback if X is hyperelliptic by computing a Weierstrass model uing IsHyperelliptic and working with functions on it
*/
    CF := ConstantField(FunctionField(x));
    
    for f in fs do
        ux := Evaluate(MinimalPolynomial(Evaluate(f, x), CF), f);
        if Valuation(ux, x) eq 1 then
            return ux;
        end if;
    end for;
    return UniformizingElement(x);
end function;

function DifferentialExpansionMatrices(omegas, x, d)
/*
Given a sequence omegas of differentials on a function field, a place x on the same function field,
and a precision d, returns a sequence of d matrices, the ith of which encodes the ith coefficient of
the Laurent series expansion of each differential of omegas at x.
*/
    F := FunctionField(x);
    CFx := ResidueClassField(x);
    CF := ConstantField(F);
    du := UniformizingDifferential(omegas, x);
    
    fs := [F!RationalFunction(omega/du) : omega in omegas]; // Why does RationalFunction make Evaluate much faster in certain cases?
    as := [Eltseq(Evaluate(f, x), CF) : f in fs];
    expansions := [Matrix(as)];
    
    if d gt 1 then
        u := UniformizingParameter(fs, x);
        uinv := 1/u;
        v := Lift(Generator(CFx, CF), x);
        
        for _ in [2..d] do
            for i -> f in fs do
                g := &+[c * v^(e-1) : e -> c in as[i]];
                fs[i] := F!RationalFunction((f - g)*uinv);
                as[i] := Eltseq(Evaluate(fs[i], x), CF);
            end for;
            Append(~expansions, Matrix(as));
        end for;
    end if;
    
    return expansions;
end function;

