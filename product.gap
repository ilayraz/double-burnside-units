# Pick a conjugate isomorphism of isom in G from the list isoms
# G: Group
# isoms: List of isomorphisms
# isom: Isomorphism
StandardConjugate := function(G, isoms, isom)
    local orbit, f, pos;
    orbit := IsomorphismConjugateOrbit(G, G, G, isom);

    for f in orbit do
        pos := Position(isoms, f);
        if pos <> fail then
            return pos;
        fi;
    od;

    Print("Failed to find standard for ", isom, " in", isoms, "\n");
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
# standard: representatives of each conjugacy class
#           of twisted diagonal subgroups of GxG
MackeyProduct := function(G, f, h, standard)
    local representatives, rep, result, results;

    results := ZeroList(standard);
    representatives := DoubleCosetRepsAndSizes(G, Range(f), Source(h));

    for rep in representatives do
        result := StarProduct(f, ConjugateIsomorphism(h, Identity(Source(h)), rep[1]));
        result := StandardConjugate(G, standard, result);

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

# Take the product of two lists of the basis of standard in G
ProductByBasis := function(G, l1, l2, standard)
    local i, j, result, results;

    results := ZeroList(standard);

    for i in [1..Size(l1)] do
        for j in [1..Size(l2)] do
            if l1[i] <> 0 and l2[j] <> 0 then
                result := l1[i] * l2[j] * MackeyProduct(G, standard[i], standard[j], standard);
                results := results + result;
            fi;
        od;
    od;

    return results;
end;

# Take a list by basis and compute the orthogonal element
OrthogonalElement := function(G, T, standard)
    local i, inverse, result;

    result := ZeroList(standard);

    for i in [1..Size(T)] do
        inverse := StandardConjugate(G, standard, InverseGeneralMapping(standard[i]));
        result[inverse] := T[i];
    od;

    return result;
end;

# Compute the semi-direct product of a normal subgroup of orthogonal bifree double burnside ring of G, T,
# with the positive embeddings of Out(G).
# standard: choice for the basis for the bifree subgroup
GroupProduct := function(G, T, standard)
    local out, aut, t, result, results, i, element, product;

    out := OuterAutmorphismGroup(G);
    results := [];

    for t in T do
        for aut in out do
            result := ZeroList(standard);

            for i in [1..Size(t)] do
                if t[i] <> 0 then
                    product := t[i] * MackeyProduct(G, standard[i], aut, standard);
                    result := result + product;
                fi;
            od;
            Add(results, result);
        od;
    od;

    return results;
end;
