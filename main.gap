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
    standard := SubgroupIsoms(G, subgroups);

    isomsTwists := FilterTwists(G, G2, standard);
    isoms := isomsTwists[1];
    twists := isomsTwists[2];

    tom := ComputeTableOfMarks(G2, Flat(twists));
    TransposedMatDestructive(tom);
    tom := Inverse(tom);

    subresults := StartWalk(G, isoms, tom);
    subresults := List(subresults, res -> tom * res);

    standard := Flat(standard);


    results := GroupProduct(G, subresults, standard);

    return [standard, results];
end;
