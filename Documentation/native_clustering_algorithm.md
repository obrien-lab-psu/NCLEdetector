# Native NCLE Clustering Algorithm

This document describes, step by step, the native non-covalent lasso entanglement (NCLE) clustering algorithm as implemented in
`ClusterNativeEntanglements.Cluster_NativeEntanglements()` in [`NCLEdetector/clustering.py`](../NCLEdetector/clustering.py).

It is based on the clustering algorithm [ALG1] in the supplemental information of the original NCLE paper
([Supplemental_algorithms_and_information.pdf](https://github.com/obrien-lab-psu/Non-covalent-Lasso-Entanglements-in-Folded-Proteins-Prevalence-Functional-Implications-and-Evolut/blob/main/Supplemental_algorithms_and_information.pdf)).
Places where the current implementation differs from that description are noted as **Deviation from ALG1**.

## Inputs

- The HQ NCLE file (pipe-delimited) with columns
  `ID|chain|i|j|crossingsN|crossingsC|gn|gc|GLNn|GLNc|TLNn|TLNc|CCbond|ENT|Quality`.
- The spatial cutoff for Step 4, set by organism:

  | Organism | Cutoff |
  |----------|--------|
  | Human    | 52     |
  | Ecoli    | 57     |
  | Yeast    | 49     |
  | Any other organism | 53 (average of the three above) |

  A user-supplied `cut_off` overrides the organism default.

If the output file `<outdir>/<ID>.csv` already exists, the function loads and returns it without re-clustering.

## Internal representation

Each entanglement is carried through the algorithm as a tuple:

```
(n_raw, i, j, r1, ..., rn, "i-j;i-j;...")
```

- `n_raw`: the number of raw NCLEs this entanglement represents.
- `i, j`: the native contact that closes the representative loop.
- `r1 ... rn`: signed crossing residues, sorted by residue number (for example `-9`, `+1107`). The sign is the chirality.
- The last field lists every raw native contact (loop) this entanglement represents.

After Steps 2A, 2B, 3 and 4, a QC check confirms that the total `n_raw` still equals the number of raw NCLEs that were loaded. If it doesn't, a `ValueError` is raised.

## Step 1: Load and group by crossing set

1. Keep only rows where `ENT == True`. If a `Quality` column exists, also require `Quality == 'High'`.
2. Combine `crossingsN` and `crossingsC` into one list and sort it by residue number, ignoring the sign.
3. Group the raw NCLEs by `(ID, signed sorted crossings)`. Every loop `(i, j)` with exactly the same signed crossing vector goes into the same group.
4. For each raw NCLE, keep `gn, gc, GLNn, GLNc, TLNn, TLNc` (and the original `crossingsN`/`crossingsC`) for the output.
5. Record which native contacts are disulfide bonds (`CCbond == True`).

Note: grouping uses `ID` and crossings but not `chain`. A multi-chain input that shares one `ID` across chains could mix chains. Single-chain inputs such as AlphaFold monomers are not affected.

## Step 2A: Pick the minimal loop per crossing set

For each crossing-set group, the representative is the loop with the smallest `j − i`. The other loops in the group are stored in the contacts string, and `n_raw` is the group size.

## Step 2B: Merge entanglements with different numbers of crossings

For each pair `(A, B)` of Step 2A representatives (all combinations, no repeats), merge them if **all** of the following are true:

1. **Close crossing:** at least one pair of same-chirality crossings (one from A, one from B) is ≤ 3 residues apart.
2. **Different counts:** A and B have different numbers of crossings.
3. **Loops overlap:**
   `iA ∈ [iB, jB]`, or `jA ∈ [iB, jB]`, or `iB ∈ [iA, jA]`, or `jB ∈ [iA, jA]`.
4. **No crossing inside the combined loop span:** no crossing from A or B lies in
   `[min(iA, jA, iB, jB), max(iA, jA, iB, jB)]`.
5. **Crossing distance ≤ 20:**
   - Strip the signs from the crossings.
   - Use the Hungarian algorithm (`scipy.optimize.linear_sum_assignment`) to match each crossing of the entanglement with fewer crossings to a different crossing of the other one, minimizing the total squared difference.
   - The distance is √(sum of squared differences) for that matching.
   - When one entanglement has a single crossing, this reduces to the smallest |a − b| between that crossing and any crossing of the other entanglement.

When all five conditions hold, the entanglement with **fewer** crossings is removed and its loops and `n_raw` are added to the entanglement with **more** crossings.

Merges can chain: if A is absorbed into B and B into C, A ends up in C. This works because the pairs are fixed from the original list, so an entanglement that has already been absorbed can still absorb others.

**Deviation from ALG1:**
- ALG1 step 3e says the distance must be **< 20**. The code uses **≤ 20**.
- ALG1 describes a "double cartesian product" distance. The code instead computes the minimum one-to-one matching distance. The Hungarian algorithm (introduced in PR #1) gives exactly the same minimum as the brute-force permutation search it replaced.

## Step 3: Merge entanglements with the same number of crossings

For each pair of entanglements left after Step 2B, merge them if **all** of the following are true:

1. They have the same number of crossings.
2. Their signs match in the same order.
3. Their loops overlap (same test as Step 2B).
4. Every index-matched crossing is close: `|rA,k − rB,k| ≤ 20` for each `k`, with no search over other matchings.

When they merge:

- The entanglement with the **shorter loop** is kept. If the loops are the same length, the second one in the pair is kept.
- The kept entanglement takes on both `n_raw` counts and both sets of contacts.
- The longer-loop entanglement is recorded as removed, and any later pair involving it is skipped. This matches ALG1's rule to skip a pair if one of its entanglements has already been merged.

## Step 4: Cluster by primary structure within groups

1. **Group** entanglements by `(ID, number of crossings, ordered chirality list)`, for example `Q9BZ72_2_['-', '-']`.
   Entanglements in different groups are never clustered together in this step.
2. **Distance** between two entanglements is the Euclidean distance over `(i, j, |r1|, ..., |rn|)`.
3. **Greedy neighbor assignment:**
   - Go through the pairs `(X, Y)` in the order they were generated.
   - If the distance is ≤ the cutoff and neither X nor Y has already been assigned as a neighbor, Y becomes a neighbor of X and is marked as assigned.
   - Each X that has neighbors becomes a cluster made of X plus its neighbors.
   - The clusters are visited by neighbor count, with random tie-breaks, but this does not change which clusters form.
4. **Singletons:** any entanglement not placed in a cluster becomes its own cluster.
5. **Representative for clusters with more than one member:**
   - Take the geometric median of the crossings (`geom_median.numpy.compute_geometric_median`).
   - Pick the members whose crossings are closest to that median, then the one among them with the smallest loop.
   - Its `n_raw` becomes the cluster total, and its contacts string becomes all of the members' contacts.

**Deviation from ALG1:**
- ALG1 describes density-based clustering with silhouette scoring. The code uses this fixed-cutoff greedy method with the organism cutoffs listed above.
- The geometric median uses **only the first crossing** of each member (`r_value[0]`), not all `n` crossings as ALG1 describes. This only affects clusters whose entanglements have two or more crossings.
- The tie-break `random.choice(np.where(...))[0]` always returns the first tied candidate, so it is not actually random.

## Step 5: Write the output

The output is `<outdir>/<ID>.csv` (pipe-delimited) with columns
`ID|chain|i|j|crossingsN|crossingsC|gn|gc|GLNn|GLNc|TLNn|TLNc|num_contacts|contacts|CCBond`.

For each representative entanglement:

- `i, j` come from the representative loop.
- Its crossings are split into `crossingsN` (residue < `i`) and `crossingsC` (residue > `j`).
- `gn, gc, GLNn, GLNc, TLNn, TLNc` come from the representative's own raw row.
- `num_contacts` is the cluster's total `n_raw`, and `contacts` is the merged contacts string.
- `CCBond = True` if any contact in the cluster is a disulfide bond.

## Worked example: Q9BZ72 (Human, AlphaFold)

Input: `DataStore/Structures/Alphafold_structures/Human/NCLE_out_mode3/Q9BZ72/HQ_NCLE/Q9BZ72.csv`

After Step 2A, the N-terminal region contains these representatives (among others):

| Representative loop | Crossings | Raw NCLEs |
|---------------------|-----------|-----------|
| (57, 243)  | -9, -21 | 8  |
| (40, 253)  | -9, -22 | 4  |
| (180, 228) | -9      | 48 |
| (166, 232) | -10     | 2  |

- **Step 2B:** `(180,228)[-9]` and `(166,232)[-10]` each meet all five conditions against `(57,243)[-9,-21]`:
  - They share a same-chirality crossing within 3 residues (−9 vs −9, and −10 vs −9).
  - The crossing counts differ (1 vs 2).
  - Their loops lie inside 57–243.
  - No crossing falls inside 57–243.
  - The matched crossing distance is 0 or 1.

  Both single-crossing entanglements are therefore absorbed into `(57,243)[-9,-21]`.
- **Step 3:** `(40,253)[-9,-22]` has the same crossing count and signs as `(57,243)[-9,-21]`, and its crossings are within 20 residues index-by-index. Because `(57,243)` has the shorter loop, it is kept.

The result is one output row, `57|243|-9,-21`, with `num_contacts = 62`. This includes the single-crossing loops 139-217 … 180-228.

This is the intended behavior of Step 2B and of ALG1 Case 1. Note (v) of ALG1 Case 1 uses this exact situation, one entanglement with two crossings and one with a single crossing, as its example. Only Step 4 requires entanglements to have the same number of crossings and the same chirality.
