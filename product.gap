StandardConjugate := function(G2, basis, isom)
    local target, twist, i;

    twist := IsomorphismToTwistedDiagonalSubgroup(G2, isom);
    for i in [1..Size(basis)] do
        target := IsomorphismToTwistedDiagonalSubgroup(G2, basis[i]);
        if IsConjugate(G2, target, twist) then
            return i;
        fi;
    od;

    Print("Failed to find basis for ", isom, " in", basis, "\n");
    return fail;
end;


# Pick a conjugate isomorphism of isom in G from the list isoms
# G: Group
# isoms: List of isomorphisms
# isom: Isomorphism
StandardConjugateOld := function(G, isoms, isom)
    local orbit, f, pos;
    orbit := IsomorphismConjugateOrbit(G, G, G, isom);

    for f in orbit do
        pos := Position(isoms, f);
        if pos <> fail then
            return pos;
        fi;
    od;

    Print("Failed to find basis for ", isom, " in", isoms, "\n");
    return fail;
end;

StarProduct := function(isom1, isom2)
    local composed, intersection, gens, images,
          aRestricted, bRestricted;
    intersection := Intersection(Source(isom1), Range(isom2));
    gens := GeneratorsOfGroup(intersection);

    images := Images(isom1, intersection);
    aRestricted := GroupHomomorphismByImagesNC(intersection,
                       images, gens,
                       List(gens, x -> Image(isom1, x)));

    images := PreImages(isom2, intersection);
    bRestricted := GroupHomomorphismByImagesNC(images, intersection,
                       List(gens, x -> PreImage(isom2, x)), gens);

    composed := CompositionMapping(aRestricted, bRestricted);
    return composed;
end;

# Tensor product of f*h in group G
# where f,h are isomorphisms
# basis: representatives of each conjugacy class
#           of twisted diagonal subgroups of GxG
MackeyProduct := function(G, G2, f, h, basis)
    local representatives, rep, result, results;

    results := ZeroList(basis);
    representatives := DoubleCosetRepsAndSizes(G, Range(f), Source(h));

    for rep in representatives do
        result := StarProduct(f, ConjugateIsomorphism(h, Identity(Source(h)), rep[1]));
        result := StandardConjugate(G2, basis, result);

        results[result] := results[result] + 1;
    od;

    return results;
end;

# Compute representatives of the outer automorphism group of G
OuterAutmorphismGroup := function(G)
    local automorphismGroup, innerAutomorphism;

    automorphismGroup := AutomorphismGroup(G);
    innerAutomorphism := InnerAutomorphismsAutomorphismGroup(automorphismGroup);

    return RightTransversal(automorphismGroup, innerAutomorphism);
end;

# Take the product of two lists of the basis in G
ProductByBasis := function(G, G2, l1, l2, basis)
    local i, j, result, results;

    results := ZeroList(basis);

    for i in [1..Size(l1)] do
        for j in [1..Size(l2)] do
            if l1[i] <> 0 and l2[j] <> 0 then
                result := l1[i] * l2[j] * MackeyProduct(G, G2, basis[i], basis[j], basis);
                results := results + result;
            fi;
        od;
    od;

    return results;
end;

# Take a list by basis and compute the orthogonal element
OrthogonalElement := function(G2, T, basis)
    local i, inverse, result;

    result := ZeroList(basis);

    for i in [1..Size(T)] do
        inverse := StandardConjugate(G2, basis, InverseGeneralMapping(basis[i]));
        result[inverse] := T[i];
    od;

    return result;
end;

# Compute the semi-direct product of a normal subgroup of orthogonal bifree double burnside ring of G, T,
# with the positive embedding of Out(G).
# basis: choice for the basis for the bifree subgroup
GroupProduct := function(G, G2, T, basis)
    local out, aut, t, result, results, i, element, product;

    out := OuterAutmorphismGroup(G);
    results := [];

    for t in T do
        for aut in out do
            result := ZeroList(basis);

            for i in [1..Size(t)] do
                if t[i] <> 0 then
                    product := t[i] * MackeyProduct(G, G2, basis[i], aut, basis);
                    result := result + product;
                fi;
            od;
            Add(results, result);
        od;
    od;

    return results;
end;
