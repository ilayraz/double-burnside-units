# Compute (Source, f, Dest)
# Goes from isomorphism to twisted diagonal subgroup
# G is direct product group
# TODO: Opposite Group shenenigans??
IsomorphismToTwistedDiagonalSubgroup := function(G, f)
    local e1, e2, generators;

    # If trivial subgroup, just return it
    if Size(Source(f)) = 1 then
        return TrivialSubgroup(G);
    fi;


    e1 := Embedding(G, 1);
    e2 := Embedding(G, 2);
    generators := List(GeneratorsOfGroup(Source(f)), el -> Image(e2, el) * Image(e1, Image(f, el)));

    return Group(generators);
end;

# Compute all isomorphisms group1 -> group2
# if isomorphic, else return fail
GroupIsomorphisms := function(group1, group2)
    local iso, autGroup, isoList;

    iso := IsomorphismGroups(group1, group2);
    if iso = fail then
        return fail;
    fi;

    autGroup := AutomorphismGroup(group1);
    isoList := List(autGroup, aut -> CompositionMapping(iso, aut));

    return isoList;
end;

# Conjugate source(f) by g, range(f) by h.
# f: Isomorphism
# g: Element in some parent group of Source(f)
# h: Element in some parent group of Range(f)
ConjugateIsomorphism := function(f, g, h)
    local source, gens;
    source := Source(f);
    gens := GeneratorsOfGroup(source);
    return GroupHomomorphismByImagesNC(source^g, Range(f)^h,
                                       List(gens, x -> x^g),
                                       List(gens,
                                            x -> Image(f, x)^h));
end;

# Compute conjugacy orbit of (H, isom, K) in G
# G,H,K: Groups
# isom: Isomorphism
IsomorphismConjugateOrbit := function(G, H, K, isom)
    local conjugates, h, k;
    conjugates := [];

    for h in H do
        for k in K do
            Add(conjugates, ConjugateIsomorphism(isom, h, k));
        od;
    od;
    return conjugates;
end;


# Compute list of non-conjugate group isomorphisms H->K
# in G
NonconjugateIsomorphisms := function(G,H,K)
    local HCoset, KCoset, isoms, isom, conjugates, result, orbit;

    conjugates := [];
    result := [];
    isoms := GroupIsomorphisms(H, K);
    if isoms = fail then
        return [];
    fi;

    HCoset := RightTransversal(
                                Normalizer(G, H),
                                Centralizer(G, H));
    KCoset := RightTransversal(
                                Normalizer(G, K),
                                Centralizer(G, K));

    # If possible do the identity isomorphism instead of its conjugate
    if H = K then
        isom := IdentityMapping(H);
        orbit := IsomorphismConjugateOrbit(G, HCoset, KCoset, isom);
        Add(result, isom);
        UniteSet(conjugates, orbit);
    fi;

    for isom in isoms do
        if not isom in conjugates then
            orbit := IsomorphismConjugateOrbit(G, HCoset, KCoset, isom);
            Add(result, isom);
            UniteSet(conjugates, orbit);
        fi;
    od;

    return result;
end;

# Compute all the valid isomorphisms from list returned by GroupPartition
SubgroupIsoms := function(G, groupsSubs)
    local isoClass, subs, H, K, i, j, isoms, result, subresult;

    result := [];
    for isoClass in groupsSubs do
        for subs in isoClass do
            subresult := [];

            for i in [1..Size(subs)] do
                H := subs[i];
                for j in [i..Size(subs)] do
                    K := subs[j];

                    isoms := NonconjugateIsomorphisms(G,H,K);

                    if i = 1 then
                        Add(subresult, isoms);

                        if i <> j then
                            Add(subresult, List(isoms, f -> InverseGeneralMapping(f)));
                        fi;
                    else
                        Append(subresult[i], isoms);

                        if i <> j then
                            Append(subresult[j], List(isoms, f -> InverseGeneralMapping(f)));
                        fi;
                    fi;
                od;
            od;
            Append(result, subresult);
        od;
    od;

    return result;
end;
