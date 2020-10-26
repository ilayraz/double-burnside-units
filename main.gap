# Import all submodules
Read("imports.gap");


# Compute all valid subgroups of G, up to normalizer/centralizer step
GroupPartition := function(G)
    local subgroups;

    subgroups := SubgroupPartitionByIsomorphism(G);
    return GroupsPartitionByNormalizerCentralizer(G, subgroups);
end;


Main := function(G)
    local G2, subgroups, isoms, twists, isomsTwists, tom, subresults, results, standard, el;

    G2 := DirectProduct(G,G);
    subgroups := GroupPartition(G);
    standard := SubgroupIsoms(G, G2, subgroups);

    isomsTwists := FilterTwists(G, G2, standard);
    isoms := isomsTwists[1];
    twists := isomsTwists[2];

    tom := ComputeTableOfMarks(G2, Flat(twists));

    Print("Aligning matrix...\n");
    # tom := Reversed(TransposedMat(Reversed(TransposedMat(tom))));
    tom := TransposedMat(tom);

    Print("Inverting matrix...\n");
    tom := Inverse(tom);

    subresults := StartWalk(G, isoms, tom);
    subresults := List(subresults, res -> tom * res);

    Print("Got results of size: ",Size(subresults), "\n");

    standard := Flat(standard);

    Print("Computing semidirect product with outer automorphism subgroup...\n");
    results := GroupProduct(G, G2, subresults, standard);

    results := subresults;

    return [standard, results, tom];
end;
