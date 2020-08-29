#Compute M_(G/V)(U) = #{gV | g in G s.t. ugV = gV for all v in V}
# U,V subsets of G
ComputeMark := function(G, U, V)
    local cosets, coset, mark;

    mark := 0;
    cosets := DoubleCosetRepsAndSizes(G,U,V);

    for coset in cosets do
        if coset[2] = Size(U) then
            mark := mark + 1;
        fi;
    od;

    return mark;
end;

# Compute table of marks of subs in G
# Assumes the first entry in subs is G itself
ComputeTableOfMarks := function(G, subs)
    local i, j, U, V, mark, matrix, row;

    matrix := [];

    for i in [1..Size(subs)] do
        row := List([1..i-1], _ -> 0);
        U := subs[i];

        for j in [i..Size(subs)] do
            V := subs[j];

            if i = j then
                mark := Index(Normalizer(G, V), V);
            else
                mark := ComputeMark(G, U, V);
            fi;
            Add(row, mark);
        od;

        Add(matrix, row);
    od;

    return matrix;
end;
