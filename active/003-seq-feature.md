---
Feature name: interval-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@RNAer", "@mortonjt", "@ebolyen", "@jairideout", "@wasade", "@gregcaporaso", "@rob-knight"]
---

# Summary

Introduce an `interval_metadata` object to allow for the storage and querying of intervals.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).  Storing intervals instead of `positional_metadata` would substantially save on memory.

Additional benefits of having a `interval_metadata` objects include the ability to rapidly query features by interval.  IntervalTrees allow for both querying by position and ranges to retrieve all overlapping intervals.  This can be particularly helpful when dealing with layers of features (i.e. transcripts made up of exons, or operons made up of genes).

# Detailed design

We propose to use the `IntervalTree` object in `bx-python.`  Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.

## `IntervalMetadata` object structure
features : dictionary.  keys are `skbio.Feature` objects.  values are a list of tuples, where each tuple corresponds to an interval.
intervals : bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `skbio.Feature` object.   

These objects are arranged such that features can be queried by both feature attributes and intervals.

## `Feature` object structure
The feature object is more a less a immutable, hashable dictionary.  This store arbiturary information about features.  For instance, strand information, CDS, ...
Since there is no standard in genbank and related formats for keywords, arbituary keywords can be encoded
```python
gene = Feature(name='sagA', type='CDS', strand='+', function='toxin')
```
The above code would encode for a gene, whose name is 'sagA', that is a coding region on the positive strand whose function is a toxin.

## `IntervalMetadata` functions
`constructor(features=None)`
- `features` : dictionary.  keys are `skbio.Feature` objects.  values are a list of tuples, where each tuple corresponds to an interval.
- Called to initialize the `IntervalMetadata` object

`reverse_complement(length)`
- `length` : int.  The length of the genome.  This is required when swapping coordinates
- This performs a reverse complement on all the coordinates with in IntervalTree.

```python
   interval_metadata = IntervalMetadata()
   interval_metadata.add(Feature(gene='sagB', cds=True), (3, 5))
   iv = interval_metadata.reverse_complement(length=10)
```

`update(old_feature, new_feature)`
- `old_feature` : `skbio.Feature`.  feature to be replaced
- `new_feature` : `skbio.Feature`.  updated feature to replace the original feature
- Since `skbio.Feature` objects are immutable, if individual features in the `IntervalMetadata` are to be updated, the need to be completely replaced.

```python
    interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
        interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
        interval_metadata.update(Feature(gene='sagB', cds=True),
                                 Feature(gene='sagB', cds=False))
```


`add(feature, *intervals)`
- `feature` : `skbio.Feature`. new feature object add into the `IntervalMetadata` object
- `*intervals` : a list of tuples corresponding to intervals where the feature object exists.
- This allows for a single feature (include those that have non-continugous features) to be added into the `IntervalMetadata` object.

```python
   interval_metadata = IntervalMetadata()
   interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
   interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
```

`query(*args, **kwargs)`
- `*args` : a list of tuples to search for features
- `**kwargs` : keyword arguments to query features by their attributes.
- This allows for features to be searched by both intervals and keywords (i.e. CDS vs exon, ...)
- Note that for a specific interval query, multiple features can be returned.
 
```python
   interval_metadata = IntervalMetadata(features={
                                        Feature(gene='sagA', location=0): [(0, 2), (4, 7)],
                                        Feature(gene='sagB', location=0): [(3, 5)]})
 
   feats = interval_metadata.query((1, 2))
   feats = interval_metadata.query(gene='sagB')
```

This is the current implementation in micronota and fits all of our use-cases at the moment.  We don't feel that IntervalIndex objects as suggested in the recent pandas poll request can easily handle all of our use cases (i.e. non-contiguous intervals isn't straightforward).

##Drawbacks
- The `query` method does a linear lookup when searching for attributes over features.  This problem can blow up (in terms of development) if we want faster querying.  Which is why we opted for slow, linear search.  Again suggestions here are welcome.
- features and intervals are loosely coupled.

#Unresolved Questions
- The `update` method may not be the most straightforward for users.  The reason why we opted for this design was because the `skbio.Feature` object can only be indexed uniquely by hashing the entire feature.  Suggestions on how to improve this will be welcome.
- Is there a better way to associate intervals to their corresponding features?
- How will features be deleted?
- How will intervals for a given feature be updated?
- How exactly will slicing be handled?  May want to be able to pass intervals/features into `__getitem__` within the `skbio.Sequence` class
- How exactly should features be added to the interval_metadata?

## Other use cases
Getting interval trees is actually fairly crucial for metadata querying.  Imagine if a user wants to query studies by a range of pH.
