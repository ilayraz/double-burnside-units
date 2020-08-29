# Paritition all non-conjugate subgroups of G
# by isomorphism classes
SubgroupPartitionByIsomorphism := function(G)
    local subgroups, groupPartition, group, rep, id, idval, list;
    subgroups := ConjugacyClassesSubgroups(G);
    groupPartition := rec();

    for group in subgroups do
        rep := Representative(group);
        id := IdGroup(rep);
        idval := CantorPairing(id[1], id[2]);

        RecordListAdd(groupPartition, idval, rep);
    od;

    list := List(RecNames(groupPartition), name -> groupPartition.(name));
    Sort(list, function(u,v) return Size(u[1]) > Size(v[1]); end);
    return list;
end;

# Take a group G and a list of lists of subgroups of G
# and partition each sublist further into sublists by size of centralizer
# and normalizer in G.
# Two subgroups will be in the same partition if they have the same
# order of both normalizer and centralizer in G
GroupsPartitionByNormalizerCentralizer := function(G, groupsList)
    local partition, subPartition, groups, group, centralizer, normalizer, cantor;

    partition := [];

    for groups in groupsList do
        subPartition := rec();
        for group in groups do
            centralizer := Size(Centralizer(G, group));
            normalizer := Size(Normalizer(G, group));
            cantor := CantorPairing(centralizer, normalizer);

            RecordListAdd(subPartition, cantor, group);
        od;

        Add(partition, List(RecNames(subPartition), name -> subPartition.(name)));
    od;

    return partition;
end;


# Compute p_1(N_{GxG}(H))
# Where f: R->S is an isomorphism
# G is a direct product group
# TODO: Opposite Group shenenigans??
ProductNormalizerTest := function(G, H)
    local normalizer, p1;

    p1 := Projection(G, 1);
    normalizer := Normalizer(G, H);
    return Set(normalizer, g -> Image(p1, g));
end;


# Filter isoms by normalizer condition.
# Also return twisted diagonal subgroup list matching isoms.
FilterTwists := function(G, G2, isomsList)
    local isoms, isom, doubleNorm, norm, twist, twists, tempIsoms, isomsRet, twistsRet;

    isomsRet := [];
    twistsRet := [];

    for isoms in isomsList do
        twists := [];
        tempIsoms := [];

        for isom in isoms do
            twist := IsomorphismToTwistedDiagonalSubgroup(G2, isom);
            doubleNorm := ProductNormalizerTest(G2, twist);
            norm := Normalizer(G, Source(isom));

            if Size(doubleNorm) = Size(norm) then
                Add(twists, twist);
                Add(tempIsoms, isom);
            fi;
        od;

        if Size(twists) <> 0 then
            Add(isomsRet, tempIsoms);
            Add(twistsRet, twists);
        fi;
    od;

    return [isomsRet, twistsRet];
end;
