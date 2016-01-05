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

Modify the `Sequence` class to incorporate the new pandas `IntervalIndex` into `positional_metadata`.

# Motivation

The `IntervalIndex` is not merged into pandas yet. To test the code below,
you will need the specific branch of `pandas` from github:

```
git clone git@github.com:shoyer/pandas.git
cd pandas
git pull origin IntervalIndex
pip install -e .
```

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB
memory to read in a E. coli genome with 260K features. And the `pandas` sparse dataframe
can alleviate the mem usage, but it is still too memory hogging. The dicussion is mainly
in this [thread](https://github.com/biocore/scikit-bio/issues/1159).

Incorporating an IntervalIndex would also enable fast query times for fetch intervals,
which is not supported by the current prototype of the `positional_metadata` at the moment.
For example:

```python
from pandas.core.interval import Interval, IntervalIndex
intervals = [(190, 255),
             (190, 255),
             (250, 2799),
             (337, 2799),
             (2801, 3733)]
iid = IntervalIndex.from_tuples(intervals)
df = pd.DataFrame(data={'feature': ['gene1', 'gene2', 'gene3', 'gene4', 'gene5']},
                  index=iid)
df.loc[254, :]
```
Because its implementation of `IntervalTree`, it will return all the features containing
coordinate of 254 very fast:
```
            feature
(190, 255]    gene1
(190, 255]    gene2
(250, 2799]   gene3
```

# Detailed design

Use `IntervalIndex` in `pandas`.  The IntervalIndex is functional for our purposes, and we don't have to even maintain it.

The `positional_metadata` dataframe has an IntervalIndex as its row index and has feature metadata as column(s) so that each row is sequence region and a feature can have multiple rows/regions. The `positional_metadata` can have duplicate row index to accommandate the scenario that a single region belongs to multiple features.

The dataframe should have two columns, one for the name of the feature so we can fetch the feature by name of `str` easily, the other for other metadata of the features. It would be argued that the feature medatadata column contains a hashable, dict-like object. Making it hashable will be convenient for some `pandas` operations. (The hashable implementation is [here](https://github.com/biocore/scikit-bio/pull/1157/files).


# Drawbacks

1. The `IntervalIndex` is still in active development and not merged into offical branch.
There are some methods that don't aren't implemented yet
(i.e. `IntervalIndex.is_non_overlapping_monotonic`). However, this is in high priority of `pandas`
milestone and will eventually be merged. We probably set up a dev branch and code based on that.
This branch can be merged back to master of scikit-bio once `IntervalIndex` gets into pandas.

2. Need to change the API of `positional_metadata`. The proposed data frame with `IntervalIndex` is  different from current implementation, which has feature in column
and position in row.

3. Need to change the slice implementation. we can update `__getitem__` so that we can fetch the sub sequence by providing an interval or a list of intervals, to make the following code work:

  ```python
  intervals = seq.positional_metadata.query('feature == gene1').index
  seq[intervals]
  ```

  In fact, we recommend to further slim down the operation to have this (even the current slice operation in `scikit-bio` is a little verbose):

  ```python
  gene1 = seq['gene1']
  ```

  This is certainly achievable by modifying `__getitem__`.

4. May need to some engineering to combine the per position metadata (such as quality scores) with the interval metadata (such as gene locations). The question is:

   do we really have other per position metadata besides quality score?

   If no, we propose to strip out the quality score from the `positional_metadata` to have another numpy array to store it.

   If yes, we propose that internally we have 2 objects, one storing the per position metadata and the other interval metadata, while keeping the similar, if not identical, API of a single `positional_metadata`. But this means a set of extra code to write into `Sequence` class.

# Alternatives

Stumped here. Any ideas?

# Unresolved questions

1. Not sure how intersections with intervals will work in `pandas` since it hasn't been implemented. In other words, will the following hold? Supposedly we can request this functionality in `pandas` github.


```python
from pandas.core.interval import Interval
gene1 = Interval(0, 100)
gene2 = Interval(50, 200)
gene3 = Interval(120, 150)
# we want something like this
# overlaped
overlap = gene1.intersect(gene2)
assert overlap == Interval(50, 100)
# not overlaped
overlap = gene1.intersect(gene3)
assert overlap == Interval(100, 120)
```
