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

We propose to store `interval_metadata` as `IntervalMetadata` object and associate it to `Sequence` (and its child classes) .

## `IntervalMetadata` object structure

This object will have 2 attributes:

1. `features`
   dict-like.  keys are `skbio.Feature` objects.  values are a list of tuples, where each tuple corresponds to an interval.

2. `intervals`
   bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `skbio.Feature` object. 

   `IntervalTree` object in `bx-python.`  Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.

These objects are arranged such that features can be queried by both feature attributes and intervals.

## `skbio.Feature` object structure
The feature object is more a less a immutable, hashable dictionary that inherits from `collections.Mapping`.  This store arbitrary information about features.  For instance, strand information, gene name, ...
Since there is no standard in genbank and related formats for keywords, arbitrary keywords can be encoded
```python
gene = Feature(name='sagA', type='CDS', strand='+', function='toxin')
```
Since the `Feature` object is hashable, every attribute within the `Feature` object contributes to the hash.
For example
```python
In [1]: from skbio.metadata import Feature

In [2]: hash(Feature(name='sagA', strand='-'))
Out[2]: -8103701482769186978

In [3]: hash(Feature(name='sagA', strand='+'))
Out[3]: 3837639023530960356
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
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(Feature(gene='sagB', cds=True), (3, 5))
   >>> iv = interval_metadata.reverse_complement(length=10)
   >>> iv.query(gene='sagB')
   >>> iv.features[f[0]]
   [(5, 7)]
```

`update(old_feature, new_feature)`
- `old_feature` : `skbio.Feature`.  feature to be replaced
- `new_feature` : `skbio.Feature`.  updated feature to replace the original feature
- Since `skbio.Feature` objects are immutable, if individual features in the `IntervalMetadata` are to be updated, the need to be completely replaced.

```python
    interval_metadata = IntervalMetadata()
        interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
        interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
        interval_metadata.update(Feature(gene='sagB', location=0),
                                 Feature(gene='sagB', location=1))
```


`add(feature, *intervals)`
- `feature` : `skbio.Feature`. new feature object add into the `IntervalMetadata` object
- `*intervals` : an iterable of tuples corresponding to intervals where the feature object exists.
- This allows for a single feature (include those that have non-continugous features) to be added into the `IntervalMetadata` object.

```python
   interval_metadata = IntervalMetadata()
   interval_metadata.add(Feature(gene='sagA', location=0), 1, (4, 7))
   interval_metadata.add(Feature(gene='sagB', location=0), (3, 5))
```

`query(*args, **kwargs)`
- `*args` : an iterable of tuples to search for features
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

#Unresolved Questions (arranged in priority)
- The `update` method may not be the most straightforward for users.  The reason why we opted for this design was because the `skbio.Feature` object can only be indexed uniquely by hashing the entire feature.  Suggestions on how to improve this will be welcome.
- Is there a better way to associate intervals to their corresponding features?
- How will features be deleted?
- How will intervals for a given feature be updated?
- How exactly will slicing be handled?  May want to be able to pass intervals/features into `__getitem__` within the `skbio.Sequence` class
- How exactly should features be added to the interval_metadata?
- What IntervalTree data structure should we stick with? 
- Should we re-invent the wheel and roll out our own IntervalTree object.


# Alternative design
Another alternative is to make `BoundFeature` objects that are tightly coupled with the `IntervalTree`.
This will ease the process of updating and deleting interval/features from the `IntervalMetadata` object.

## `BoundFeature` object structure
This object is a *mutable* object, that contains attributes (i.e. `name`, `function`, `strand`, ...) as well as a list of intervals.  This object would have a [weakref](https://docs.python.org/3/library/weakref.html) to the corresponding interval object in the `IntervalTree`.  if the intervals within `BoundFeature` are updated, the intervals within the `IntervalTree` are updated. We will probably want to have a reserved keyword for `intervals` in the `BoundFeature` object.  For instance, if we wanted to grab all of the intervals associated with a features, we should be able to run something as follows
```python
f = BoundFeature([1, (4, 7)], gene='sagA', function='toxin')
f.intervals
```
This means, that the user should not be able to pass in a `intervals` keyword argument.


`update(*args, **kwargs)`
- `*args`: list of intervals that need to be swapped in
- `**kwargs`: List of attributes that will be modified for the given `BoundFeature` object.
- Updates the corresponding entries in the `IntervalTree`, if intervals are modified, through the weakref

```python
f = BoundFeature(name='sagA', function='transport')
f.update(function='toxin')
f.update((1, 2))
```

`__del__`
- Deletes a `BoundFeature`.  The corresponding interval(s) in the `IntervalTree` will be removed when this is called.


## `IntervalMetadata` object structure
### Attributes
* `features`: a list of `BoundFeature` object
* `_update_flag`: boolean. Whenever any intervals of any member feature is modified, this is set to True, indicating `_intervaltree` needs to be updated.
* `_intervaltree` (implemented as property): `IntervalTree` object created from all the intervales in the `features`. When this attributes is accessed by any methods (like `query` below), it checks `_update_flag` to decide whether to re-create from all intervals before it returns the interval tree. This lazy update saves computation.

### Methods
`constructor(features=None)`
- `features` : list of `BoundFeature` objects.
- Called to initialize the `IntervalMetadata` object
- The weakrefs will be updated here to link the `BoundFeature` object to the current `IntervalMetadata` object

`reverse_complement(length)`
- `length` : int.  The length of the genome.  This is required when swapping coordinates
- This performs a reverse complement on all the coordinates with in IntervalTree.
- Same as original design

`add(feature)`
- `feature` : `skbio.BoundFeature`. new feature object add into the `IntervalMetadata` object
- This allows for a single feature (include those that have non-continugous features) to be added into the `IntervalMetadata` object.
- The weakref of the added `skbio.BoundFeature` will be updated

```python
   interval_metadata = IntervalMetadata()
   interval_metadata.add(BoundFeature([1, (4, 7)], gene='sagA', function='toxin'))
   interval_metadata.add(BoundFeature([(3, 5)], gene='sagB', function='toxin'), )
```

`query(*args, **kwargs)`
- `*args` : an iterable of tuples to search for features
- `**kwargs` : keyword arguments to query features by their attributes.
- This allows for features to be searched by both intervals and keywords (i.e. CDS vs exon, ...)
- Note that for a specific interval query, multiple features can be returned.
- Same as the original design

The mutability of the `BoundFeature` enable us to operate directly on the feature. For example:
```python
feature_list = interval_metadata.query(gene='ark')
for gene in feture_list:
    gene.GO = 'GO0003243'
# Now all ark genes in the `interval_metadata` are updated. we don't need to interject the updated genes back into `interval_metadata` like we do for the imutable implementation.
```
