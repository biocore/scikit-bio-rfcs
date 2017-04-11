---
Feature name: interval-metadata
Start date: <2015-12-18>
Pull request: 1
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@RNAer", "@mortonjt", "@ebolyen", "@jairideout", "@wasade", "@gregcaporaso", "@rob-knight"]
---

# Summary

Introduce an `IntervalMetadata` object to allow for the storage, modification, and querying of features covering a region of a biological sequence. A feature is defined as a sub-region of a sequence that is a functional entity, e.g., a gene, a riboswitch, an exon, etc.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse data frame can alleviate the mem usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).  Storing intervals instead of `positional_metadata` would substantially save on memory.

Additional benefits of having a `IntervalMetadata` objects include the ability to rapidly query features by coordinates.  Interval trees allow for both querying by position and ranges to retrieve all overlapping intervals.  This can be particularly helpful when dealing with layers of features (i.e. transcripts made up of exons, or operons made up of genes).  In addition, this metadata could be used in other objects such as TabularMSA to store coding information within each alignment.

# Detailed design

We propose 2 new public data objects: `BoundFeature` and `IntervalMetadata`. `BoundFeature` stores all the information of a feature. `IntervalMetadata` stores all the `BoundFeatures` objects of a sequence and add it as `interval_metadata` attribute (as similar to `positional_metadata`) in `Sequence` (and its child classes).

## `BoundFeature` object
This object is a *mutable*, dict-like object to store the all the info of a sequence feature. This object would also have a reference to the corresponding `IntervalMetadata` object. The reference enables the tight coupling between `BoundFeature` and `IntervalMetadata`. The mutability of `BoundFeature` enables us to directly modify a `BoundFeature` object from an `IntervalMetadata` object. Note that this reference to `IntervalMetadata` object is a private attribute to `BoundFeature`, but not exposed to the public api.

### attributes:
* `intervals`: store intervals of coordinates. This is represented as a list of tuples of a pair of ints
* `boundaries`: records the openness of each interval.  So a boundary of `(True, False)` would it indicate that the exact right boundary is unknown, corresponding to the examples [here](ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3).
* `_data`: Dictionary of attributes, this is private object, storing other information of the feature like strand, gene_name, product
* `_ref`: reference to the `IntervalMetadata` object.



### methods
`__init__(data, interval, boundaries=None)`
The construction would be like:
```python
>>> f = BoundFeature(intervals=[(1,2), (4,7)], boundaries=None, data={'name':'sagA', 'function':'transport'})
>>> f['name']
'sagA'
```

`update(data)`
- `data`: Dictionary of keywords and attributes that will be modified for the given `BoundFeature` object.
```python
# note the “intervals” lives in `BoundFeature._data`, it is not the same with `BoundFeature.intervals`.
>>> f = BoundFeature(intervals=[(1, 2)], data={‘name’:'sagA', ‘function’:'transport', ‘intervals’:’some_info’})
>>> f.update({‘function’:’toxin’, ‘intervals’:’other_info’})
```


`__getitem__(key)`
- `key`: str
- Retrieves the value of a key
```python
>>> f = BoundFeature(intervals=[(1,2), (4,7)], boundaries=None, data={'name':'sagA', 'function':'transport'})
>>> f['name']
'sagA'
```

`__setitem__(key, val)`
- `key`: str
- `val`: the value of the keyword value to update to.
- set the value of a key
```python
>>> f = BoundFeature(name='sagA', function='transport')
>>> f['name'] = 'sagB'
```

`__eq__(other)`
`other` : `BoundFeature`
Tests for equality over all attributes and intervals


`@property`
intervals
Getter property for intervals.

```python
>>> f = BoundFeature(intervals=[(1, 2)], data={‘name’:'sagA', ‘function’:'transport', ‘intervals’:’some_info’})
>>> f.intervals
[(1, 2)]
```

`@intervals.setter(intervals)`
Setter property for intervals.
Sets the _is_stale_tree to true within the corresponding interval metadata

```python
>>> f = BoundFeature(intervals=[(1, 2)], data={‘name’:'sagA', ‘function’:'transport', ‘intervals’:’some_info’})
>>> f.intervals = [(10, 12)]
>>> f.intervals
[(10, 20)]
```

## `IntervalMetadata` object
### Attributes
* `features`: a list of unique `BoundFeature` objects, stored in an unordered fashion.
* `_is_stale_tree`: boolean. Whenever any intervals of any member feature is modified, this is set to `True`, indicating `_intervaltree` needs to be updated.  This is set False whenever the interval tree is retrieved through its getter.  This is where the IntervalTree is rebuilt.
* `_intervaltree` (implemented as property): `IntervalTree` object created from all the intervals in the `features`. When this attributes is accessed by any methods (like `query` below), it checks `_is_stale_tree` to decide whether to re-create from all intervals before it returns the interval tree and set `_is_stale_tree` to False. This lazy update saves computation.

   This is implemented as bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `BoundFeature` object.

   Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.


### Methods
`__init__(self, features=None)`
- `features` : list of `BoundFeature` objects.
- Called to initialize the `IntervalMetadata` object
- The references will be updated here to link each `BoundFeature` object to the current `IntervalMetadata` object

`add(data, intervals, boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `boundaries` : an iterable of boundary tuples to create features
- `data` : a single dictionary
- This has the same API of `BoundFeature` initializer and creates the `BoundFeature` objects inside
- Inserts a list of `BoundFeature` objects into the `IntervalTree`
- set `_is_stale_tree` to True.
```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(3, 5)], boundaries=None, data={‘name’:'sagB', ‘cds’:True, ‘function’:'transport')
   >>> interval_metadata.add(intervals=[(3, 5)], data={gene:'sagB', cds=True})
```
`drop(intervals=None, criteria=None, how='intersect', boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `criteria` : a dictionary to query features by their attributes.
- `how`: `str` specify `intersect`, or `union` of queries results
- `boundaries` : an iterable of boundary tuples to search for features
- Drops all `BoundFeature` objects that matches the query.
- Set `_is_stale_tree` to True
```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], boundaries=None, data={'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(40, 70)], boundaries=None, data={'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], boundaries=None, data={'gene':'sagB'})
   >>> interval_metadata.drop(criteria={‘gene’:’sagA’})
   >>> feats = interval_metadata.query(intervals=(1, 2))
   >>> feats
   None
```
Alternative `drop(features)`:
- In the 1st `drop()` design, we are concerned that users may drop additional features that they are not intended to and unaware of this dropping, if those features also have the same data. What do you guys think?
- `features`: an iterable of `BoundFeature` objects
- Compare every object in `features` to all the `BoundFeatures` within self.features and drops those that are exactly equal
- Set `_is_stale_tree` to True
```python
   >>> interval_metadata = IntervalMetadata(features)
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], boundaries=None, data={'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], boundaries=None, data={'gene':'sagB'})
   >>> feats = interval_metadata.query(intervals=(0, 2), data={‘gene’:’sagA’})
   >>> interval_metadata.drop(feats)
   >>> feats = interval_metadata.query(intervals=(0, 2), data={‘gene’:’sagA’})
   >>> feats
   None
```

`query(intervals=None, criteria=None, how='intersect', boundaries=None)`
- `intervals` : an iterable of interval tuples to search for features
- `criteria` : a dictionary to query features by their attributes.
- `boundaries` : an iterable of boundary tuples to search for features
- `how`: `str` specify `intersection`, or `union` of queries results
- This allows for features to be searched by both intervals and attributes (i.e. CDS vs exon, ...)
- Note that for a specific query, multiple features can be returned.
- Returns a generator of `BoundFeature` objects that satisfy the search criteria
- If _is_stale_tree is False, the IntervalTree is rebuilt and _is_stale_tree is set to True

```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], data={'gene':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], data={'gene':'sagB'})
   >>> feats = interval_metadata.query(intervals=(1, 2))
   >>> feats = interval_metadata.query(criteria={'gene_name':'sagB'})
   >>> feats
   BoundFeature(intervals=[(3, 4)], data={‘gene’:’sagB’})
```

`reverse()`
- This performs a reverse on all the coordinates within `IntervalTree`. Often associated with reverse complement.
- Relies on the length method to be implemented in the parent class
- Set `_is_stale_tree` to True

```python
   >>> iv = IntervalMetadata()
   >>> assert len(iv) == 10  # length inherited from the class (i.e. Sequence)
   True
   >>> iv.add(intervals=[(3, 5)], data={‘name’:'sagB', ‘cds’:True, ‘function’:'transport'})
   >>> iv = iv.reverse()
   >>> f = iv.query(criteria={‘gene’:'sagB'})
   >>> f.intervals
   [(5, 7)]
```

The mutability of the `BoundFeature` enable us to operate directly on the feature. For example:
```python
>>> for gene in feats = interval_metadata.query(criteria={gene:'sagB'}):
...     gene['GO'] = 'GO0003243'
# Now all sagB genes in the `interval_metadata` are updated. we don't need to interject the updated genes back into `interval_metadata` like we do for the previously proposed immutable implementation.
```
