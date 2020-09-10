#Compute M_(G/V)(U) = #{gV | g in G s.t. ugV = gV for all v in V}
# U,V subsets of G
ComputeMark := function(G, V, U)
    local cosets, coset, mark;

    # Print("Computing mark of ", U, "with ", V, "\n");

    mark := 0;
    cosets := DoubleCosetRepsAndSizes(G,V,U);

    for coset in cosets do
        if coset[2] = Size(V) then
            mark := mark + 1;
        fi;
    od;

    return mark;
end;

# Compute table of marks of subs in G
# Assumes the first entry in subs is G itself
ComputeTableOfMarks := function(G, subs)
    local i, j, U, V, mark, matrix, row, normalizers, count;

    Print("Computing table of marks... for ", Size(subs), " elements\n");

    # We need all the normalizers so compute them now
    normalizers := List(subs, sub -> Normalizer(G, sub));

    matrix := [];
    count := 0;

    for i in [1..Size(subs)] do
        row := List([1..i-1], _ -> 0);
        V := subs[i];

        for j in [i..Size(subs)] do
            U := subs[j];

            if Size(U) = 1 then
                mark := Index(G, V);
            elif Size(V) = Size(G) then
                mark := 1;
            elif Size(V) mod Size(U) <> 0 then
                mark := 0;
            elif IsNormal(G, U) or IsNormal(G,V) then
                if IsSubset(V,U) then
                    mark := Index(G,V);
                else
                    mark := 0;
                fi;
            elif i = j then
                mark := Index(normalizers[i], V);
            else
                mark := ComputeMark(G, V, U);
                count := count + 1;
            fi;

            Add(row, mark);
        od;

        Add(matrix, row);
    od;

    Print(count, " elements computed the hard way\n");

    return matrix;
end;
