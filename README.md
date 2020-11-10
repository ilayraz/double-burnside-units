Source code for algorithm to find all orthogonal bifree units of the double Burnside ring of a group.
This algorithm makes heavy usage of [this](https://boltje.math.ucsc.edu/publications/OrthogonalUnits.pdf) text by
Robert Boltje, and Philipp Perepelitsky. Particularly on lemma 3.1.

The algorithm works as such:
Given a group `G`:
1. Find all non-conjugate subgroups `H1, H2, ..., Hn`  of `G`.
2. Partition subgroups in sets where all elements in the set are isomorphic, and have the same order of both their normalizer and centralizer in `G`.
3. Find all non-conjugate twisted diagonal subgroups of `G^2`. That is, subgroups that can be described as `Δ(Hi, φ, Hj)` where φ is an isomorphism `Hj->Hi`, and the subgroup itself is `{(φ(h), h) | h in Hj}`.
4. Remove all twisted diagonal subgroups `H=Δ(Hi, φ, Hj)` that don't satisfy the condition `|p_1(N_(G^2)(H)| = |N_G(Hj)|`
5. Compute the sub-table of the table of marks of `G^2` with only the remaining twisted diagonal subgroups for entries (ordered in the same way as our subgroups will be in step 6), then compute its inverse
6. Order the partitions by the size of the subgroups, consider all options as follows:
- Starting with the first partition, Pick an unpicked subgroup `H=Δ(Hi, φ, Hj)` in it.
- For the given subgroup, consider the element under the subset mark homomorphism in the same sense as the table of marks described in step 5. The value of the vector at the column for `H` must be `|N_G(Hj)|` or `-|N_G(Hj)|` (Consider both paths independently), and have `0` at the columns of every other subgroup of the form `Δ(Hk, φ, Hj)` for any `k`.
- If the vector constructed above constitutes a valid unit, then it would be an element in the double burnside ring, so considering it under the inverse table of marks would yield an integer number in each entry in the resulting vector. Furthermore, since the table of marks is a triangular matrix, we can check each step in the path only with the previous choices we made, so if at any points a path fails, we can terminate it early.
- Recursively repeat this process until all branches either reach the end or terminate early. All branches that reached the end are valid orthogonal bifree units, and in fact this gives us all such units.

We can make a slight optimization to step 6 by only considering the positive branch of the identity automorphism on `G`.
That is, only considering paths where the first entry in the mark vector is `|C(G)|` the center of `G`.
Using this, we won't find all units, but a normal subgroup `T` of the orthogonal bifree unit subgroup.
If we let `B(G)` be the subgroup of units from the single burnside ring embedded into `B(G,G)`, then the
the orthogonal bifree unit group is the semi-direct product of `T` with `B(G)`.
