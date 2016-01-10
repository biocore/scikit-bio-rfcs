---
Feature name: sequence-positional-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@RNAer", "@mortonjt"]
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

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).

Incorporating an IntervalIndex would also enable fast query times for fetch intervals,which is not supported by the current prototype of the `positional_metadata` at the moment. For example:

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

Because its implementation of `IntervalTree`, it will return all the features containing coordinate of 254 very fast:
```
            feature
(190, 255]    gene1
(190, 255]    gene2
(250, 2799]   gene3
```

# Detailed design

Use `IntervalIndex` in `pandas`. It is still in active development and not merged into offical branch.There are some methods that don't aren't implemented yet(i.e. `IntervalIndex.is_non_overlapping_monotonic`). However, this is in high priority of `pandas` milestone and will eventually be merged.

1. We will have per position data frame (`positional_metadata`) and interval data frame (`interval_metadata`) to store the sequence metadata.

2. The `interval_metadata` dataframe has an `IntervalIndex` as its row index, feature type, strand, left open, right open, and other feature metadata as a dict.

3. Each row in `interval_metadata` represents contiguous sequence region. A feature can span multiple rows. An interval can duplicate in multiple rows for different features. The map from feature to row/interval is multiple-to-multiple.

## Types of queries
1. fetch a interval.

2. fetch a feature spanning multiple intervals, such as an gene of multiple exons. This will be queried by passing a keyword, like "gene='ABC transporter X'" to get all the intervals of the same gene. This works similarly as well to fetch all features with a common attributes (such as all genes of ABC transporters).

# Drawbacks

2. Need to change the slice implementation. we can update `__getitem__` so that we can fetch the sub sequence by providing an interval or a list of intervals, to make the following code work:
  ```python
  intervals = seq.positional_metadata.query('feature == gene1').index
  seq[intervals]
  ```

  In fact, we recommend to further slim down the operation to have this (even the current slice operation in `scikit-bio` is a little verbose):

  ```python
  gene1 = seq['gene1']
  ```

  This is certainly achievable by modifying `__getitem__`.


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
