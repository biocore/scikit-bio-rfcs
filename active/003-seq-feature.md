---
Feature name: sequence-positional-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors:
- "@RNAer"
- "@mortonjt"
---

# Summary

Modify the genbank io to incorporate the new pandas `IntervalIndex`

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes 1TB
memory to read in a E. coli genome with 260K features. And the `pandas` sparse dataframe
can alleviate the mem usage, but it is still too memory hogging. The dicussion is mainly
in this [thread](https://github.com/biocore/scikit-bio/issues/1159).

Incorporating an IntervalIndex would also enable fast query times for fetch intervals, 
which is not supported by the current prototype of the `positional_metadata` at the moment.


# Detailed design

Use `IntervalIndex` in `pandas`.  The IntervalIndex is fully functional for our purposes, and we don't have to even maintain it.
The `positional_metadata` class doesn't need to be altered at all.  All we need to do is pass in a `pd.DataFrame` with a different index.

The suggested main differences to the genbank io reader would be as follows

1. The `positional metadata` dataframe can have an IntervalIndex.
2. The feature name will be a column in the `positional_metadata`.  This is to be able to link the interval back to the feature metadata contained in the hashable object of column names. Discussed here https://github.com/biocore/scikit-bio/pull/1157/files
3. The hashable object of column names will contain an interval to link to `positional_metadata`.

# Drawbacks
The `pd.IntervalIndex` is not stable.  There are some methods that don't aren't implemented yet (i.e. intersection between intervals and indeces).

# Alternatives
Stumped here.  Any ideas?

# Unresolved questions

- The `IntervalIndex` is not merged into `pandas` yet.
- Not sure if distance between intervals will be incorporated into `pandas`.
- How large of an issue will memory requirements be?  One nice thing about the `IntervalIndex` is not every base pair position needs to be documented, so it allows for row compression to some extent.  On the other hand, there is a possibility that there could be many metadata columns in the `positional_metadata` dataframe.  Then we may want to explore using a sparse data frame with the `IntervalIndex`.  Will definitely need to do more memory profiling.
- Not sure how intersections with intervals will work in `pandas` since it hasn't been implemented. In other words, will the following hold?
```python
from pandas.core.interval import Interval
gene1 = Interval(0, 100)
gene2 = Interval(50, 200)
# we want this
overlap = gene1.intersect(gene2)
assert overlap == Interval(50, 100)
```