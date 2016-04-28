---
Feature name: interval-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@RNAer", "@mortonjt", "@ebolyen", "@jairideout", "@wasade", "@gregcaporaso", "@rob-knight"]
---

# Summary

Introduce an `interval_metadata` object to allow for the storage, modification, and querying of features covering a region of a biological sequence. A feature is defined as a sub-region of a seuqence that is a functional entity, e.g., a gene, a riboswitch, an exon, etc.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).  Storing intervals instead of `positional_metadata` would substantially save on memory.

Additional benefits of having a `interval_metadata` objects include the ability to rapidly query features by coordinates.  Interval trees allow for both querying by position and ranges to retrieve all overlapping intervals.  This can be particularly helpful when dealing with layers of features (i.e. transcripts made up of exons, or operons made up of genes).

# Detailed design

We propose 2 new public data objects: `BoundFeature` and `IntervalMetadata`. `BoundFeature` stores all the attributes of a feature. `IntervalMetadata` stores all the `BoundFeatures` objects of a sequence and add it as `interval_metadata` attribute (as similar to `positional_metadata`) in `Sequence` (and its child classes).

## `BoundFeature` object
This object is a *mutable* object, that contains atribtrary attributes (i.e. `intervals`, `gene_name`, `function`, `strand`, ...) to store the all the info of a sequence feature. This object would also have a [weakref](https://docs.python.org/3/library/weakref.html) to the corresponding `IntervalMetadata` object. If the intervals of a `BoundFeature` are updated, the interval tree within the `IntervalMetadata` object is also auto updated. The weak reference enables the tight coupling between `BoundFeature` and `IntervalMetadata`. The mutability of `BoundFeature` enables us to directly modify a `BoundFeature` object from a `IntervalMetadata` object.

The `intervals`, `strand`, `wref` are enforced attributes to store coordinates, strand, and a weak reference. The construction would be like:
```python
>>> f = BoundFeature(intervals=[1, (4, 7)], strand='+', wref=None, gene='sagA', function='toxin')
>>> f.intervals  # get coordinates
>>> f.wref   # get weak ref to the interval metadata
>>> f.function   # get the feature info
```

### methods
`update(**kwargs)`
- `**kwargs`: List of attributes that will be modified for the given `BoundFeature` object.
- Updates the corresponding entries in the `IntervalMetadata`, if intervals are modified, through the weakref

```python
f = BoundFeature(name='sagA', function='transport')
f.update(function='toxin')
f.update(intervals=[(1, 2)])
```


## `IntervalMetadata` object
### Attributes
* `features`: a list of `BoundFeature` objects
* `_staled`: boolean. Whenever any intervals of any member feature is modified, this is set to True, indicating `_intervaltree` needs to be updated.
* `_intervaltree` (implemented as property): `IntervalTree` object created from all the intervales in the `features`. When this attributes is accessed by any methods (like `query` below), it checks `_staled` to decide whether to re-create from all intervals before it returns the interval tree. This lazy update saves computation.

   This is implemented as bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `BoundFeature` object. 

   Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.

### Methods
`__init__(self, features=None)`
- `features` : list of `BoundFeature` objects.
- Called to initialize the `IntervalMetadata` object
- The weakrefs will be updated here to link each `BoundFeature` object to the current `IntervalMetadata` object

`reverse_complement(length)`
- `length` : int.  The length of the whole sequence.  This is required when swapping coordinates
- This performs a reverse complement on all the coordinates with in IntervalTree.
- Set `_staled` to True

```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(Feature(gene='sagB', cds=True, intervals=[(3, 5)])
   >>> iv = interval_metadata.reverse_complement(length=10)
   >>> f = iv.query(gene='sagB')
   >>> f.intervals
   [(5, 7)]
```

`add(feature)`
- `feature` : `BoundFeature`. Add a new feature into the `IntervalMetadata` object
- This allows for a single feature (include those that have non-continugous features) to be added into the `IntervalMetadata` object.
- The weakref of the added `BoundFeature` will be updated
- set `_staled` to True.
```python
   interval_metadata = IntervalMetadata()
   interval_metadata.add(BoundFeature([1, (4, 7)], gene='sagA', function='toxin'))
   interval_metadata.add(BoundFeature([(3, 5)], gene='sagB', function='toxin'), )
```

`drop(feature)`
- Deletes a `BoundFeature`.
- Set `_staled` to True

`query(*args, **kwargs)`
- `*args` : an iterable of interval tuples to search for features
- `**kwargs` : keyword arguments to query features by their attributes.
- This allows for features to be searched by both intervals and attributes (i.e. CDS vs exon, ...)
- Note that for a specific interval query, multiple features can be returned.

```python
   interval_metadata = IntervalMetadata(features={
                                        BoundFeature(gene='sagA', intervals=[(0, 2), (4, 7)]),
                                        BoundFeature(gene='sagB', intervals=[(3, 5)]}))
 
   feats = interval_metadata.query((1, 2))
   feats = interval_metadata.query(gene='sagB')
```

The mutability of the `BoundFeature` enable us to operate directly on the feature. For example:
```python
for gene in feats = interval_metadata.query(gene='sagB'):
    gene.GO = 'GO0003243'
# Now all sagB genes in the `interval_metadata` are updated. we don't need to interject the updated genes back into `interval_metadata` like we do for the previouly proposed imutable implementation.
```


## Drawbacks
- The `query` method does a linear lookup when searching for attributes over features. (The query by intervals are still super fast). This problem can blow up (in terms of development) if we want faster querying by attributes. We don't think there is a need for tons of queries by attributes. Thus, we sacrifice this for much better API design by implementing `BoundFeature` as mutable rather than immutable.
