# Import all submodules
Read("imports.gap");


# Compute all valid subgroups of G, up to normalizer/centralizer step
GroupPartition := function(G)
    local subgroups;

    subgroups := SubgroupPartitionByIsomorphism(G);
    return GroupsPartitionByNormalizerCentralizer(G, subgroups);
end;


# Compute sublist whose basis is composed only of identity morphisms
# Note that this function assumes that if an element in the basis
# is conjugate to an identity morphism, then the element choosen is the identity morphism itself
EmbeddedSubgroup := function(basis, G)
    local identityCoordinates;

    identityCoordinates := List(basis, f -> IdentityMapping(Source(f)) = f);
    return Filtered(G,
                    g -> ForAll([1..Size(basis)],
                                i -> g[i] = 0 or identityCoordinates[i]));
end;


Main := function(G)
    local G2, subgroups, isoms, twists, isomsTwists, tom, results, basis, el;

    G2 := DirectProduct(G,G);
    subgroups := GroupPartition(G);
    isoms := SubgroupIsoms(G, G2, subgroups);

    isomsTwists := FilterTwists(G, G2, isoms);
    isoms := isomsTwists[1];
    twists := isomsTwists[2];

    tom := ComputeTableOfMarks(G2, Flat(twists));

    Print("Aligning matrix...\n");
    tom := TransposedMat(tom);

    Print("Inverting matrix...\n");
    tom := Inverse(tom);

    results := StartWalk(G, isoms, tom);
    results := List(results, res -> tom * res);

    Print("Got results of size: ",Size(results), "\n");

    basis := Flat(isoms);

    return [basis, results];
end;
