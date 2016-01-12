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

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).


# Detailed design

We propose to use `IntervalIndex` in `pandas` to store the interval features. It is still in active development and not merged into offical branch. There are some methods that don't aren't implemented yet(i.e. `IntervalIndex.is_non_overlapping_monotonic`). However, this is in high priority of `pandas` milestone and will eventually be merged.

To run some of the code below, you will need the specific branch of `pandas` from github:

```
git clone git@github.com:shoyer/pandas.git
cd pandas
git pull origin IntervalIndex
pip install -e .
```

1. We will have per position data frame (`positional_metadata`) and interval data frame (`interval_metadata`) to store the sequence metadata.

2. The `interval_metadata` dataframe has following index or columns:
   1) an `IntervalIndex` as its row index (named `INTERVAL`).

   2) feature type (`TYPE`). Type of the features, like `gene`, `CDS`, etc.

   3) strand (`STRAND`). Controled str volcabulary of `+`, `-`). `+` means the current strand and `-` means its anti-sense strand.
   4) left open (`LEFT` with `True` or `False` value), indicating if the exact lower boundary point of a feature is unknown. It would be `True` for the GenBank case of `<345..500`, as an example.
   5) right open (`RIGHT` with boolean values as `LEFT`)
   6) final column - all the other feature metadata as a dict (`ATTRIBUTES`).

   We'd like `TYPE`, `STRAND`, `LEFT`, and `RIGHT` as indiviual columns instead of storing them all in `ATTRIBUTES` dict. The idea is that we want to strip all the essential attributes for a seq feature from the `ATTRIBUTES` dict into columns of `interval_metadata` so that `ATTRIBUTES` is less cluttered.

   If any data in the columns is missing, leave it as empty in the data frame.

3. Each row in `interval_metadata` represents contiguous sequence region. A feature can span multiple rows. An interval can duplicate in multiple rows for different features. The map from feature to row/interval is multiple-to-multiple. In the following example, the `cds1` feature contains 0-30 and 70-100 intervales while both `gen3` and `gen4` are located in the 120-150 interval.

   ```python
>>> from pandas.core.interval import Interval, IntervalIndex
>>> intervals = [(0, 100), (0, 30), (70, 100), (50, 200), (120, 150), (120, 150)]
>>> iix = IntervalIndex.from_tuples(intervals, closed='left')
>>> iix
IntervalIndex(left=[0, 0, 70, 50, 120, 120],
              right=[100, 30, 100, 200, 150, 150],
              closed='left')
>>> gene1 = dict(name='gen1',
...                 gene_synonym="ECK0007; JW0006")
>>> cds1 = dict(name='cds1')
>>> gene2 = dict(name='gen2',
...                 db_xref=("EcoGene:EG11555", "GeneID:944757"))
>>> gene3 = dict(name='gen3',
...                 db_xref="EcoGene:EG11555")
>>> gene4 = dict(name='gen4')
>>> interval_md = pd.DataFrame(data=np.array([['gene', 'CDS', 'CDS', 'gene', 'gene', 'gene'],
...                                  ['+', '+', '+', '+', '+', '+'],
...                                  [True, True, True, True, True, True],
...                                  [True, True, True, True, True, False],
...                                  [gene1, cds1, cds1, gene2, gene3, gene4]]).T,
...                            copy=True,
...                            columns=['TYPE', 'STRAND', 'LEFT', 'RIGHT', 'ATTRIBUTES'],
...                            index=iix)
>>> interval_md
            TYPE STRAND  LEFT  RIGHT  \
[0, 100)    gene      +  True   True
[0, 30)      CDS      +  True   True
[70, 100)    CDS      +  True   True
[50, 200)   gene      +  True   True
[120, 150)  gene      +  True   True
[120, 150)  gene      +  True  False

                                                   ATTRIBUTES
[0, 100)    {'name': 'gen1', 'gene_synonym': 'ECK0007; JW0...
[0, 30)                                      {'name': 'cds1'}
[70, 100)                                    {'name': 'cds1'}
[50, 200)   {'db_xref': ('EcoGene:EG11555', 'GeneID:944757...
[120, 150)     {'db_xref': 'EcoGene:EG11555', 'name': 'gen3'}
[120, 150)                                   {'name': 'gen4'}
   ```

## Types of feature queries
* query by location. `IntervalIndex` module implements `IntervalTree` for fast query by location. To get features covering a unit or a region, we can do this:

``` python
>>> interval_md.loc[30, ]
          TYPE STRAND  LEFT RIGHT  \
[0, 100)  gene      +  True  True

                                                  ATTIBUTES
[0, 100)  {'name': 'gen1', 'gene_synonym': 'ECK0007; JW0...
>>> interval_md.loc[Interval(10, 30), ]
          TYPE STRAND  LEFT RIGHT  \
[0, 100)  gene      +  True  True
[0, 30)    CDS      +  True  True

                                                 ATTRIBUTES
[0, 100)  {'name': 'gen1', 'gene_synonym': 'ECK0007; JW0...
[0, 30)                                    {'name': 'cds1'}
```

* query by attributes. It can be used to fetch multiple features or a single feature. For example:

``` python
>>> interval_md[interval_md.apply(lambda row: row['ATTRIBUTES']['name']=='cds1', axis=1)]
          TYPE STRAND  LEFT RIGHT        ATTRIBUTES
[0, 30)    CDS      +  True  True  {'name': 'cds1'}
[70, 100)  CDS      +  True  True  {'name': 'cds1'}

```

## Operations that need modification

There are some methods in `Sequence` and its child classes that needs to operate on `interval_metadata` accordingly. For all the situations below, in case users only want to rc the sequence and do not intend to use `interval_metadata`, we should provide an option to skip changes and only return the sequence without `interval_metadata`.

* Reverse operation. `INTERVAL` also needs to be updated:
``` python
# original interval for a 9 nt gene
>>> a = Interval(2, 11, closed='left')
# after reverse complemenary of a seq of 12 nt, the interval should be
>>> a_rc = 12 - a
>>> a_rc
Interval(1, 10, closed='left')
```

* Complementary operation. The `STRAND` needs to be switched.

* Concatenation. The intervals of all the sequences, except for the 1st one, needs to be updated.

``` python
# interval of a feature in seq1 of 10 nt
>>> a1 = Interval(2, 9, closed='left')
# interval of a feature in seq2
>>> a2 = Interval(1, 3, closed='left')
# after concat, the two features will have the intervals of:
>>> a1  # still the same
Interval(2, 9, closed='left')
>>> a2 + 10
Interval(11, 13, closed='left')
```

* Slicing. The intervals needs to reflect the change of the sequence. Only keep the features that overlaps with sliced sequence, and need to update the coordiate as well. One caveat is that if a feature is cut in the middle of the slicing with missing left or right ends, the `LEFT` or `RIGHT` needs to be changed to `TRUE`

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
