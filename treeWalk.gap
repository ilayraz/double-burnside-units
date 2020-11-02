# Compute the result of multiplying a matrix with a vector at the last element of vector only
# mat: lower-triangular matrix
ElementMultiply := function(mat, vector)
    local result, length, row, i, element;

    result := 0;
    length := Size(vector);
    row := mat[length];

    for i in [1..length] do
        element := row[i] * vector[i];
        result := result + element;
    od;

    return result;
end;

# Copy vector, check if the product at element is an integer.
# If it is, return vector, else return fail.
# At the end, pad vector with numPad zeroes
CheckComponent := function(mat, vector, element, numPad)
    local newVector, product;

    newVector := ShallowCopy(vector);
    Add(newVector, element);
    product := ElementMultiply(mat, newVector);

    if IsInt(product) then
        Append(newVector, Zero([1..numPad]));
        return newVector;
    else
        # Print("Got nonintegral result (", product, ") for ", newVector, "\n");
        return fail;
    fi;
end;


# Walk the next level of the tree
# G: group in which marks are computed
# twistsList: List of lists of twisted diagonal isomorphisms to check.
#             Sorted by size. Every group should be a subgroup of the first entry.
# tom: table of marks used to verify
# usedGroups: Set of groups already used as targets
# vector: vector of the current chain in the tree
WalkTree := function(G, twistsList, tom, usedGroups, vector)
    local candidateList, twistsSize, i, twist, centralSize, newVector, newSet;

    # Base case
    if IsEmpty(twistsList) then
        return [vector];
    fi;


    candidateList := [];
    twistsSize := Size(twistsList[1]);

    for i in [1..twistsSize] do
        twist := twistsList[1][i];

        if not Range(twist) in usedGroups then
            centralSize := Size(Centralizer(G, Source(twist)));

            # Check positive
            newVector := CheckComponent(tom, vector, centralSize, twistsSize - i);
            if newVector <> fail then

                newSet := ShallowCopy(usedGroups);
                AddSet(newSet, Range(twist));

                Append(candidateList, WalkTree(G, twistsList{[2..Size(twistsList)]}, tom, newSet, newVector));
            fi;

            # Check negative
            newVector := CheckComponent(tom, vector, -1 * centralSize, twistsSize - i);
            if newVector <> fail then
                newSet := ShallowCopy(usedGroups);
                AddSet(newSet, Range(twist));

                Append(candidateList, WalkTree(G, twistsList{[2..Size(twistsList)]}, tom, newSet, newVector));
            fi;
        fi;

        # Add another 0 for each twist we passed
        Add(vector, 0);
    od;

    return candidateList;
end;

# Start descending tree
# G: The group in which marks are computed
# twistsList: List of lists of twisted diagonal isomorphisms to check.
#             Sorted by size.
# tom: table of marks used to verify

StartWalk := function(G, twistsList, tom)
    local center, usedGroups, vector;

    Print("Starting to walk tree...\n");

    center := Size(Center(G));
    usedGroups := [G];

    vector := ZeroList(twistsList[1]);
    vector[1] := center;

    return WalkTree(G, twistsList{[2..Size(twistsList)]}, tom, usedGroups, vector);
end;
