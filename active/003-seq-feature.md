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

Replace the current implementation of `positional_metadata` in the `Sequence` class.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes 1TB
memory to read in a E. coli genome with 260K features. And the `pandas` sparse dataframe
can alleviate the mem usage, but it is still too memory hogging. The dicussion is mainly
in this [thread](https://github.com/biocore/scikit-bio/issues/1159).


# Detailed design

Using `IntervalIndex` in `pandas` can be the ideal solution. It will not add burden for
scikit-bio to maintain the extra code for this.

# Drawbacks



# Alternatives



# Unresolved questions

The `IntervalIndex` is not merged into `pandas` yet.

How the current API will change if we use `IntervalIndex`?
