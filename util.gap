# Compute th Cantor Pairing function
CantorPairing := function(k1, k2)
    return (k1+k2)*(k1+k2+1)/2 + k2;
end;

# Check if record has key (val)
# If not, insert an a list at that value
# Then add elm to the list at that value
RecordListAdd := function(record, val, elm)
    if not IsBound(record.(val)) then
        record.(val) := [];
    fi;
    Add(record.(val), elm);
end;

# Returns a list of the same size with all zeroes
ZeroList := function(lst)
    return List(lst, x -> 0);
end;
