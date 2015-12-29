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
2. The `FEATURE` list in the `metadata` attribute will become a `FEATURE` dictionary where the feature names are keys.  This is to ensure that features can be directly looked up via the feature name.
3. The feature name will be a column in the `positional_metadata`.  This is to be able to link the interval back to the feature metadata contained in `FEATURE` field within `metadata` dictionary.
4. The `FEATURE` dictionary will have a `pd.Interval` in place of the `LOCATION` value.  This is to ensure easy look up in the `positional_metadata` dataframe.

# Drawbacks
None that I can think of.  The api for `metadata` or `positional_metadata` doesn't need to be altered at all.

# Alternatives
Stumped here.  Any ideas?

# Unresolved questions

The `IntervalIndex` is not merged into `pandas` yet.
Not sure how intersections with intervals will work in `pandas` since it hasn't been implemented.
Not sure if distance between intervals will be incorporated into `pandas`.

