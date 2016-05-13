---
Feature name: tabular-msa
Start date: 2015-08-14
Pull request: 1
Authors: ["@jairideout", "@ebolyen"]
Contributors:
- "@gregcaporaso"
- "@wasade"
- "@rob-knight"
- "@shiffer1"
---
# Summary

Replace the current experimental `skbio.alignment.Alignment` class with a
`TabularMSA` class, providing an API for manipulating tabular (row-column)
multiple sequence alignments.

# Motivation

The current multiple sequence alignment (MSA) class assumes that it is the only
way to represent an alignment. There are multiple ways of representing MSAs,
including [partial-order MSAs](http://bioinformatics.oxfordjournals.org/content/18/3/452.short)
and [A-Bruijn MSAs](http://www.ncbi.nlm.nih.gov/pmc/articles/PMC525693/). The
reason these cannot be the same object is that a row-column MSA exists as a
total order (which means that positions exist), while partial-order MSAs do not
inherently have positions or gaps (because the MSA is represented as a directed
acyclic graph, or DAG). Thus, while some operations are shared among these MSA
structures, there are intrinsic differences that necessitate different APIs.

In scikit-bio 0.4.0 the sequence classes were completely redesigned and the
current MSA class (`Alignment`) makes invalid assumptions about the API of the
sequences it stores (for example, see
[#1036](https://github.com/biocore/scikit-bio/issues/1036)). This gives us an
opportunity to rework the MSA API to better fit the new sequence classes.

After these changes, we will have an obvious way to introduce new MSA
structures (e.g., partial-order MSA), better interaction with the new sequence
API, pandas-style vectorization, and flexible, arbitrary metadata that can
represent file formats such as
[Stockholm](http://sonnhammer.sbc.su.se/Stockholm.html).

# Detailed design

A `TabularMSA` object is a mutable, monomorphic, ordered collection of
`IUPACSequence` objects that are of equal length. As its name implies,
`TabularMSA` represents a table-like structure where the first dimension is
sequences and the second dimension is positions in the alignment. Thus, each
row corresponds to a sequence and each column corresponds to a position. The
first dimension is optionally indexed by an arbitrary key. Both dimensions are
indexable via ordinals.

The axes are named `sequence` (`0`) and `position` (`1`). Both the string and
integer form should be usable anywhere the axis is specified.

One of the advantages of `TabularMSA` over the existing `Alignment` class is
vectorization support: boolean and fancy indexing. Currently, to filter out
alignment positions that have a Shannon entropy greater than 0.5:

```python
aln.subalignment(positions_to_keep=[
    i for i, v in enumerate(aln.position_entropies()) if v <= 0.5])
```

With `TabularMSA`:

```python
aln[:, aln.entropies() <= 0.5]
```

Another example is filtering out completely-gapped positions. With the current
API, this will (probably) work:

```python
aln.omit_gap_positions(1.0 - np.finfo(float).eps)
```

With `TabularMSA`:

```python
aln[:, aln.gap_frequencies(relative=True) < 1.0]
```

Thus with boolean vector support, `TabularMSA` is significantly more flexible
in its filtering capabilities.

## Constructors
- `TabularMSA(iterable, key=None)`

   Construct a `TabularMSA` from an iterable of `IUPACSequence` objects.
   `key` is a callable or a lookup key in `IUPACSequence.metadata`. When
   `key=None` the MSA will not have any keys and certain properties may not
   be defined and will raise exceptions.

- `from_dict(dict)`

   Return a new TabularMSA from the keys and sequences of a dict. The values
   must all be an appropriate `IUPACSequence` class of the same length.

## Properties
- `dtype -> IUPACSequence`

   Type of sequence stored (e.g., `DNA`)

- `shape -> Shape(int, int)`

   Implemented as a namedtuple: `namedtuple('Shape', ['sequence', 'position'])`
   Named tuple storing number of rows (sequences) and columns (positions)

- `loc -> indexer`

   Allows indexing via keys (like pandas). Will raise exception if used when
   the MSA has no keys.

- `iloc -> indexer`

   Allows indexing via indices (like pandas)

- `metadata -> dict`

   Arbitrary metadata that applies to the alignment as a whole.

- `positional_metadata -> pd.DataFrame`

   Arbitrary metadata that applies to the positions of the alignment as a
   whole.

- `keys -> np.ndarray`

   Vector of keys in the order of the MSA. Will raise exception if used when
   the MSA has no keys.

## Methods
- `__nonzero__() -> bool`

   Return `False` if the MSA has no sequences OR if all sequences are empty.

- `__len__() -> int`

   Return number of sequences in MSA. This is equivalent to `aln.shape[0]`

- `__iter__() -> iterator`

   Iterator of `IUPACSequence` objects in the order specified by `aln.keys`

- `__reversed__() -> iterator`

   Iterator of `IUPACSequence` objects in reverse order of `aln.keys`

- `__getitem__(indexer) -> TabularMSA | Sequence`

   Accepts 1 or 2 axes separated by a comma. The rules for axes are applied
   from left to right. The least significant axis can only restrict the type
   returned, it cannot take what would have been `DNA` and return a
   `TabularMSA`. This restriction makes it easier to mentally model the
   indexing behavior. So the following must always hold:

   ```python
   aln[X, Y] == aln[X][:, Y] # or [Y] if the result is 1D
   ```

   #### First axis:
   If the index is a single key or an ordinal, return an `IUPACSequence` of
   type `dtype`:

   ```python
   aln['a'] -> DNA
   aln[0] -> DNA
   ```

   If the index is a slice or other complex structure of keys or ordinals,
   return `TabularMSA`:

   ```python
   aln[2:4] -> TabularMSA
   aln['b':'e':2] -> TabularMSA
   aln[['a', 'c']] -> TabularMSA
   aln[[1]] -> TabularMSA
   aln[[True, False, False, True]] -> TabularMSA
   ```

   In the event that the keys are integers, any index composed of them will be
   assumed to be in relation to keys and not indices
   (like `pandas.DataFrame.ix`). It is advisable to use `iloc` and `loc` in
   these situations. If there are no keys for the MSA, then the index will
   always be in relation to the indices.

   #### Second axis:
   The second axis will accept only column indices (e.g ordinals, slices,
   etc.). If the first axis has already narrowed the type returned to a
   `IUPACSequence` object, then the second axis is simply treated as an
   additional filter to that object:

   ```python
   aln[2, 4:7] -> DNA
   aln[0, 0] -> DNA
   ```

   Otherise:

   ```python
   aln[2:3, 4:7] -> TabularMSA
   ```

   If the type so far is a `TabularMSA` and the second axis is an ordinal, then
   *exactly* the type of `Sequence` will be returned. The exact `dtype` is not
   returned because the sequence does not represent a real biological sequence
   but rather a column in the alignment. A `Sequence` object is returned
   because it is convenient to reuse the existing numpythonic-string API that
   `Sequence` classes posses.

   ```python
   aln[:, 3] -> Sequence
   aln['a':'d', 3] -> Sequence
   ```

- `to_dict() -> dict`

   Returns a dictionary that maps the keys of the MSA to the aligned sequences.

- `iter_positions(reverse=False) -> iterator`

   Iterator of *exactly* `Sequence` objects in the order of positions in the
   alignment. If `reverse=True`, yields positions in reverse order.

- `append(IUPACSequence) -> None`

   Add a sequence to the end of the alignment. Must be the same `dtype` and
   length must match second dimension. It should be documented that this
   operation, if successful, is not necessarily meaningful.

- `extend(iterable | TabularMSA) -> None`

   Add sequences to the end of the alignment. Each sequence must be the same
   `dtype` and length must match second dimension. It should be documented
   that this operation, if successful, is not necessarily meaningful.

- `sort(key=None, reverse=False) -> None`

   Sort sequences in-place. If `key=None` the existing keys will be used to
   sort. Otherwise `key` matches the constructor. If there are no keys for the
   MSA and `key=None`, an exception should be raised.

- `join(TabularMSA, how='strict') -> TabularMSA`

   Joins two `TabularMSA` by sequence (horizontally). `how` may be
   `'strict'`, in which case both MSAs must have the same keys (or be the same
   length if there are no keys). `how` may also be one of the common SQL joins:
   `'inner'`, `'outer'`, `'left'` (outer), or `'right'` (outer); these require
   keys or they will raise an exception. `how` is a superset of the interface
   provided by pandas. Outer based modes will insert gaps for missing keys.
   This should be documented as *not* preforming an alignment.

- `consensus() -> IUPACSequence`

   Return majority consensus sequence matching alignment's `dtype`. Other
   types of consensus may be added later.

- `gap_frequencies(axis='position', relative=False) -> np.ndarray`

   Return vector of gap frequencies along an axis.

- `conservation(string | function) -> np.ndarray`

   Return vector of positional conservation scores based on some stat from
   `skbio.alignment.stats`. For example, one possible conservation score is
   Shannon entropy.

- `reindex(key=None) -> None`

   Assign new keys to sequences in the alignment. `key` matches constructor.

## Misc

By design, `TabularMSA` will not be able to store generic `Sequence` objects
as they lack fundamental properties needed for an alignment
(such as knowledge of gaps). We considered moving the notion of gaps to
`Sequence`, but that made the entire sequence hierarchy feel like a leaky
abstraction. Instead we propose creating a factory to make it easy to construct
arbitrary `IUPACSequence` subclasses which will require the user to specify
what gap characters are. This has the nice side-effect of forcing a real type
on the `TabularMSA` when it is read from `io` instead of having a default with
no validation and a (slightly) crippled API.

In order to better reflect the non-biological nature of `Sequence` we also
recommend changing the `repr` to be different from `IUPACSequence`.

Adding [`frequencies()`](https://github.com/biocore/scikit-bio/issues/1074) as
a method to `Sequence` will improve the value of returning `Sequence` object
for the column.

Addition of `skbio.alignment.stats` subpackage to store different conservation
scoring algorithms (such as Shannon entropy).

# Drawbacks

- scikit-bio users will have to update their code to use a new API.

- The new class name `TabularMSA` may not be as obvious as to its behavior as
the current class name `Alignment`.

- It is unclear if it is advisable to use a tabular (row-column) MSA as it
forces an ordering where none exists. See
[Lee et al. 2002](http://bioinformatics.oxfordjournals.org/content/18/3/452.short).

# Alternatives

We considered representing an MSA internally as a DAG in support of
partial-order MSAs. However, representing positions in such an MSA defeats the
entire purpose of using a DAG at all.

We considered *only* supporting partial-order MSAs. However, since scikit-bio
is a bioinformatics data munging library, we need to support row-column MSAs
because most bioinformatics tools use or produce these structures.

We considered the following alternative names for the new class:

- `LinearMSA`: it is not immediately obvious that this is a row-column
alignment
- `TotalOrderMSA`: even less obvious and it sounds "better" than partial-order
if users aren't clear on the definitions
- `RowColMSA` / `RowColumnMSA`: painful to type and hard to remember (was it
"col" or "column"?)
- `StandardAlignment`: implies "best" or preferred alignment structure
- `Alignment`: too general as it doesn't describe an MSA, also implies a
default alignment structure
- `MatrixMSA`: the name is ambiguous and could be confused with other alignment
structures such as an alignment scoring matrix. It could also be confused with
MSAs represented by graphs because graphs can be stored as matrices (e.g.,
adjacency, incidence). Does not indicate a labeled structure. May imply too
much about the internal implementation.

# Unresolved questions
