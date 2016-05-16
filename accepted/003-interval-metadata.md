---
Feature name: interval-metadata
Start date: 2015-12-18
Pull request: 9
Authors: ["@RNAer", "@mortonjt"]
Contributors: ["@ebolyen", "@jairideout", "@wasade", "@gregcaporaso", "@rob-knight", "@josenavas"]
---

# Summary

Introduce an `IntervalMetadata` object to allow for the storage, modification, and querying of features covering a region of a biological sequence. A feature is defined as a sub-region of a sequence that is a functional entity, e.g., a gene, a riboswitch, an exon, etc.

# Motivation

The current implementation of `positional_metadata` in the `Sequence` class takes > 1TB memory to read in a E. coli genome with 260K features. And the `pandas` sparse DataFrame can alleviate the memory usage, but it is still too memory hogging. The discussion is mainly in this [thread](https://github.com/biocore/scikit-bio/issues/1159).  Storing intervals instead of `positional_metadata` would substantially save on memory.

Additional benefits of having a `IntervalMetadata` objects include the ability to rapidly query features by coordinates.  Interval trees allow for both querying by position and ranges to retrieve all overlapping intervals.  This can be particularly helpful when dealing with layers of features (i.e. transcripts made up of exons, or operons made up of genes).  In addition, this metadata could be used in other objects such as TabularMSA to store coding information within each alignment.

# Detailed design

We propose 2 new public objects: `Interval` and `IntervalMetadata`. `Interval` stores all the information of a feature. `IntervalMetadata` stores all the `Intervals` objects of a sequence and add it as `interval_metadata` attribute (as similar to `positional_metadata`) in `Sequence` (and its child classes).

## `Interval` object
This object is a *mutable*, dict-like object to store the all the info of a sequence feature. This object would also have a reference to the corresponding `IntervalMetadata` object. The reference enables the tight coupling between `Interval` and `IntervalMetadata`. The mutability of `Interval` enables us to directly modify a `Interval` object from an `IntervalMetadata` object. Note that this reference to `IntervalMetadata` object is a private attribute to `Interval`, but not exposed to the public api.

### Attributes
* `intervals`: store intervals of coordinates. This is represented as a list of tuples of a pair of ints.
* `boundaries`: records the openness of each interval.  So a boundary of `(True, False)` would it indicate that the exact right boundary is unknown, corresponding to the examples [here](ftp://ftp.ebi.ac.uk/pub/databases/embl/doc/FT_current.html#3.4.3).
* `metadata`: Dictionary of attributes storing information of the feature like strand, gene_name, product.
* `_interval_metadata`: reference to the `IntervalMetadata` object.



### Methods
`__init__(interval_metadata, intervals, boundaries=None, metadata=None)`
Construction should only be performed by IntervalMetadata (or in test code).
```python
>>> f = Interval(interval_metadata=IntervalMetadata(), intervals=[(1,2), (4,7)], boundaries=None, metadata={'name':'sagA', 'function':'transport'})
```

`drop()`
sets _interval_metadata to None and set intervals and boundaries to empty list. If _interval_metadata is None, this is a no-op.

`__getitem__(key)`
- `key`: hashable
- Retrieves the value of a key in metadata
```python
>>> f = Interval(interval_metadata=IntervalMetadata(), intervals=[(1,2), (4,7)], metadata={'name':'sagA', 'function':'transport'})
>>> f['name']
'sagA'
```

`__setitem__(key, val)`
- `key`: hashable
- `val`: the value to set metadata[key] to
- set the value of a key
```python
>>> f = Interval(interval_metadata=IntervalMetadata(), intervals=[(1,2), (4,7)], metadata={'name':'sagA', 'function':'transport'})
>>> f['name'] = 'sagB'
```

`__eq__(other)`
- `other` : `Interval`
- Confirms intervals, boundaries and metadata are equal.

`@intervals.setter(intervals)`
- Setter property for intervals.
- Sets the ``_is_stale_tree`` to true within the corresponding interval metadata

```python
>>> f = Interval(interval_metadata=IntervalMetadata(), intervals=[(1, 2)], metadata={'name':'sagA', 'function':'transport'})
>>> f.intervals = [(10, 12)]
>>> f.intervals
[(10, 12)]
```

## `IntervalMetadata` object
### Attributes
* `_is_stale_tree`: boolean. Whenever any intervals of any member feature is modified, this is set to `True`, indicating `_intervaltree` needs to be updated.  This is set False whenever the interval tree is retrieved through its getter.  This is where the IntervalTree is rebuilt.
* `_intervaltree` (implemented as property): `IntervalTree` object created from all the intervals in the `features`. When this attributes is accessed by any methods (like `query` below), it checks `_is_stale_tree` to decide whether to re-create from all intervals before it returns the interval tree and set `_is_stale_tree` to False. This lazy update saves computation.

   This is implemented as bx-python `IntervalTree` object.  The keys correspond to intervals and the values correspond to a single `Interval` object.

   Some benchmarks have been made for various IntervalTree data structures - in this [benchmark](https://gist.github.com/shoyer/c939325f509d7c027949), bx-python was determined to have one of the fastest query times.


### Methods
`__init__(self)`
- Called to initialize the `IntervalMetadata` object. Intervals are subsequently added with IntervalMetadata.add.

`add(intervals, boundaries=None, metadata=None)`
- See Interval constructor
- set `_is_stale_tree` to True.
```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(3, 5)], metadata={'name':'sagB', 'function':'transport')
```
`drop(intervals=None, boundaries=None, metadata=None)`
- `intervals` : an iterable of interval tuples to search for features
- `metadata` : a dictionary to query features by their metadata.
- `boundaries` : an iterable of boundary tuples to search for features
- Drops all `Interval` objects that matches the query (this is an intersection).
- Set `_is_stale_tree` to True
```python
   >>> interval_metadata = IntervalMetadata()
   >>> interval_metadata.add(intervals=[(0, 2), (4, 7)], boundaries=None, metadata={'name':'sagA'})
   >>> interval_metadata.add(intervals=[(40, 70)], boundaries=None, metadata={'name':'sagA'})
   >>> interval_metadata.add(intervals=[(3, 4)], boundaries=None, metadata={'name':'sagB'})
   >>> interval_metadata.drop(metadata={'name':'sagA'})
   >>> feats = list(interval_metadata.query(intervals=[(1, 2)]))
   >>> feats
   []
```

`query(intervals=None, boundaries=None, metadata=None)`
- `intervals` : an iterable of interval tuples to search for features
- `metadata` : a dictionary to query features by their metadata.
- `boundaries` : an iterable of boundary tuples to search for features
- This allows for features to be searched by both intervals and metadata (i.e. CDS vs exon, ...)
- Note that for a specific query, multiple features can be returned.
- Returns a generator of `Interval` objects that satisfy the search criteria
- If `_is_stale_tree` is False, the IntervalTree is rebuilt and `_is_stale_tree` is set to True

```python
>>> interval_metadata = IntervalMetadata()
>>> interval_metadata.add(intervals=[(0, 2), (4, 7)], metadata={'name':'sagA'})
>>> interval_metadata.add(intervals=[(3, 4)], metadata={'name':'sagB'})
>>> feats = list(interval_metadata.query(metadata={'gene_name':'sagB'}))
>>> feats
[Interval(intervals=[(3, 4)], metadata={'name':'sagB'})]
```

`_reverse(length)`
This is private, and only called from the mixed-in class reverse/reverse complement methods (which can correctly provide their length).
- This performs a reverse on all the coordinates within `IntervalTree`. Often associated with reverse complement.
- Set `_is_stale_tree` to True

## Editorial comments

`BoundFeature` is now called `Interval` - we think this is a better name (now that we stripped away Feature last week) as `IntervalMetdata` is a collection of intervals and their associated metadata.  

`_data` (private) is now `metadata` (public) for consistency with other types of metadata and ease of use. Related: we dropped `Interval.update` in favor of just using `interval.metadata.update`.

`IntervalMetadata.features` was removed since all of the information is accessible through query (and a blank query returns all of the intervals).

Include the standard "dunder" (double under) methods like ``__eq__``, ``__ne__``, ``__copy__`` - those will go into the code, but do not need to be added to this RFC.

For now you should **not** define ``__len__``, ``__iter__``, and ``__getitem__`` as these will open up a lot of new design questions, and they do not sound necessary for the current functionality.

``IntervalMetadata.__init__`` does not take features, those are added with IntervalMetadata.add (you cannot construct the individual Interval objects without the reference to the IntervalMetadata object).

Renamed ``criteria`` in ``drop`` and ``query`` to ``metadata``. This way you have parameters that correspond directly to the three ``Interval`` attributes.

Dropped ``how`` from ``drop``, because it is not clear if this applies to intervals, metadata, or boundaries. We think you should just implement what you think is the most important use case, and we can work this out when we need alternatives.

``drop(features)`` alternative: we added an ``Interval.drop``, which can be called on results of ``IntervalMetadata.query``, so we think you should keep the query-like drop API.

``reverse`` should become private. Since ``__len__`` would be ambiguous here (and open a can of worms), skip it for now, and just pass the length from the "parent object" to ``_reverse`` when you call it from ``reverse_complement``.
