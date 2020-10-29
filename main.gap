# Import all submodules
Read("imports.gap");


# Compute all valid subgroups of G, up to normalizer/centralizer step
GroupPartition := function(G)
    local subgroups;

    subgroups := SubgroupPartitionByIsomorphism(G);
    return GroupsPartitionByNormalizerCentralizer(G, subgroups);
end;


Main := function(G)
    local G2, subgroups, isoms, twists, isomsTwists, tom, results, standard, el;

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

    results := StartWalk(G, isoms, tom);
    results := List(results, res -> tom * res);

    Print("Got results of size: ",Size(results), "\n");

    standard := Flat(standard);

    return [standard, results];
end;
