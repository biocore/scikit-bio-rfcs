---
Feature name: interval-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@RNAer", "@mortonjt", "@ebolyen", "@jairideout", "@wasade", "@gregcaporaso", "@rob-knight"]
---

# Summary

Introduce an `interval_metadata` object to allow for the storage, modification, and querying of features covering a region of a biological sequence. A feature is defined as a sub-region of a sequence that is a functional entity, e.g., a gene, a riboswitch, an exon, etc.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).  Storing intervals instead of `positional_metadata` would substantially save on memory.

Additional benefits of having a `interval_metadata` objects include the ability to rapidly query features by coordinates.  Interval trees allow for both querying by position and ranges to retrieve all overlapping intervals.  This can be particularly helpful when dealing with layers of features (i.e. transcripts made up of exons, or operons made up of genes).  In addition, this metadata could be used in other objects such as TabularMSA to store coding information within each alignment.

# Detailed design

We propose 2 new public data objects: `BoundFeature` and `IntervalMetadata`. `BoundFeature` stores all the attributes of a feature. `IntervalMetadata` stores all the `BoundFeatures` objects of a sequence and add it as `interval_metadata` attribute (as similar to `positional_metadata`) in `Sequence` (and its child classes).

## `BoundFeature` object
This object is a *mutable* object, that contains arbitrary attributes (i.e. `gene_name`, `product`, ...) to store the all the info of a sequence feature, with the exception of a few required attributes. This object would also have a reference to the corresponding `IntervalMetadata` object. If the intervals of a `BoundFeature` are updated, the interval tree within the `IntervalMetadata` object is also auto updated. The reference enables the tight coupling between `BoundFeature` and `IntervalMetadata`. The mutability of `BoundFeature` enables us to directly modify a `BoundFeature` object from a `IntervalMetadata` object.  Note that there will be private reference, to `IntervalMetadata` object.  This is a private attribute, but not exposed to the public api

### required arguments
Besides user arbitrarily given attributes, we enforce the following attributes:
* `intervals`: store intervals of coordinates.  This is represented as a list of tuples of a pair of ints
* `boundaries`: records the openness of each interval.  So a boundary of `(True, False)` would it indicate that the exact right boundary is unknown, corresponding to the examples [here](ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3).
* `data`: Dictionary of attributes


### methods
`__init__(interval, boundaries, data)`
The construction would be like:
```python
>>> im = IntervalMetadata()
>>> im.add(intervals=[(1,2), (4,7)], boundaries=None, {'foo': 'bar', 'abc': [1,2,3], 'intervals': 'yeah!!'}) # this is adding a single dict, could support an iterable of dicts
>>> bound_feature = next(im.query('foo')) # bound_feature of type BoundFeature
>>> bound_feature.intervals # returns [(1,2), (4,7)]
>>> bound_feature['foo'] # returns 'bar'
>>> bound_feature['abc'] # returns [1,2,3]
>>> bound_feature['intervals'] # returns 'yeah!!'
# there should also be a __setitem__
>>> bound_feature['new-thing'] = 42
```

Note, in the above example, the interval `1` is shorthand for `(1, 2)`.

`update(**kwargs)`
- `**kwargs`: List of attributes that will be modified for the given `BoundFeature` object.
- Updates the corresponding entries in the `IntervalMetadata`, through the private reference


```python
>>> f = BoundFeature(name='sagA', function='transport')
>>> f.update(function='toxin')
>>> f.update(intervals=[(1, 2)])
```

`__getitem__(kwd)`
- `kwd`: str, the keyword argument
- Retrieves the value of a keyword attribute

```python
>>> f = BoundFeature(intervals=[(1,2), (4,7)], boundaries=None, {'name':'sagA', 'function':'transport'})
>>> f['name']
'sagA'
```

`__setitem__(kwd, val)`
- `kwd`: str, the keyword argument
- `val`: the value of the keyword value to update to.
- Retrieves the value of a keyword attribute

```python
>>> f = BoundFeature(name='sagA', function='transport')
>>> f['name'] = 'sagB'
```

## `IntervalMetadata` object
### Attributes
* `features`: a list of unique `BoundFeature` objects, stored in an unordered fashion.
* `_is_stale_tree`: boolean. Whenever any intervals of any member feature is modified, this is set to `True`, indicating `_intervaltree` needs to be updated.
* `_intervaltree` (implemented as property): `IntervalTree` object created from all the intervales in the `features`. When this attributes is accessed by any methods (like `query` below), it checks `_is_stale_tree` to decide whether to re-create from all intervals before it returns the interval tree. This lazy update saves computation.

   This is implemented as bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `BoundFeature` object. 

   Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.

### Methods
`__init__(self, features=None)`
- `features` : list of `BoundFeature` objects.
- Called to initialize the `IntervalMetadata` object
- The references will be updated here to link each `BoundFeature` object to the current `IntervalMetadata` object


`add(intervals, features, boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `boundaries` : an iterable of boundary tuples to create features
- `features` : an iterable of dictionaries objects.
- Create the `BoundFeature` objects inside
- Inserts a list of `BoundFeature` objects into the `IntervalTree`
- This allows for multiple features (including those that have non-contiguous intervals) to be added into the `IntervalMetadata` object.
- The references of the added `BoundFeature` objects will be updated
- set `_is_stale_tree` to True.
```python
   >>> len(iv) = 10  # from the parent class
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(
   >>> BoundFeature(intervals=[(3, 5)], boundaries=None, name='sagB', cds=True, function='transport')
   >>> interval_metadata.add(BoundFeature(gene:'sagB', cds=True, intervals=[(3, 5)])
```

`drop(intervals=None, keywords=None, how='intersect', boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `keywords` : a dictionary to query features by their attributes.
- `how`: `str` specify `intersect`, or `union` of queries results
- `boundaries` : an iterable of boundary tuples to search for features
- Drops all `BoundFeature` objects that matches the query.
- Set `_is_stale_tree` to True

```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], boundaries=None, {'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], boundaries=None, {'gene':'sagB'})
   >>> interval_metadata.drop(interval=(1, 2))
   >>> feats = interval_metadata.query(intervals=(1, 2)) # returns None
   
```

`query(intervals=None, keywords=None, how='intersect', boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `keywords` : a dictionary to query features by their attributes.
- `boundaries` : an iterable of boundary tuples to search for features
- `how`: `str` specify `intersection`, or `union` of queries results
- This allows for features to be searched by both intervals and attributes (i.e. CDS vs exon, ...)
- Note that for a specific interval query, multiple features can be returned.
- Returns a generator of `BoundFeature` objects that satisfy the search criteria

```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], boundaries=None, {'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], boundaries=None, {'gene':'sagB'})
   >>> feats = interval_metadata.query(intervals=(1, 2))
   >>> feats = interval_metadata.query(keywords={'gene_name':'sagB'})
```

`reverse()`
- This performs a reverse complement on all the coordinates with in IntervalTree.
- Relies on the length method to be implemented in the parent class 
- Set `_is_stale_tree` to True

```python
   >>> len(iv) = 10  # from the parent class
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(3, 5)], boundaries=None, name='sagB', cds=True, function='transport')
   >>> iv = interval_metadata.reverse()
   >>> f = iv.query(gene='sagB')
   >>> f.intervals
   [(5, 7)]
```

The mutability of the `BoundFeature` enable us to operate directly on the feature. For example:
```python
>>> for gene in feats = interval_metadata.query(gene='sagB'):
>>>     gene['GO'] = 'GO0003243'
# Now all sagB genes in the `interval_metadata` are updated. we don't need to interject the updated genes back into `interval_metadata` like we do for the previouly proposed imutable implementation.
```

## Drawbacks
- The `query` method does a linear lookup when searching for attributes over features. (The query by intervals are still super fast). This problem can blow up (in terms of development) if we want faster querying by attributes. We don't think there is a need for tons of queries by attributes. Thus, we sacrifice this for much better API design by implementing `BoundFeature` as mutable rather than immutable.


