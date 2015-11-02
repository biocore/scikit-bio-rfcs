---
Feature name: py3-only
Start date: 2015-09-23
Pull request: 6
Authors: ["@jairideout", "@gregcaporaso"]
Contributors:
- "@wasade"
- "@anderspitman"
- "@ElDeveloper"
- "@jni"
- "@Jorge-C"
- "@josenavas"
- "@ebolyen"
- "@JWDebelius"
---

# Summary

Remove Python 2 (Legacy Python) support from scikit-bio in favor of Python 3
support only.

# Motivation

We propose removing Python 2 (Legacy Python) support from scikit-bio for the
following reasons:

1. The scientific Python community is in favor of dropping support for Python 2
as soon as possible (e.g., see [recent discussion at DS4DS](https://bids.hackpad.com/Planning-an-early-death-for-Python-2-Legacy-Python.-6aB7agi5XCr)).
Dropping Python 2 support in scikit-bio will force Python 2 and Python 2/3
packages depending on scikit-bio to also drop Python 2 support if they continue
depending on scikit-bio, encouraging adoption of this critical change.
2. Dropping Python 2 support will reduce development and maintenance burden on
scikit-bio developers. For example, we can remove the hacks in `skbio.io`
necessary for Python 2 compatibility. These hacks were particularly
time-consuming to put in place in order to support Python 2 and 3. For example,
[this hack](https://github.com/biocore/scikit-bio/blob/0.4.0/skbio/io/_fileobject.py#L31-L80)
was non-trivial and required reading CPython's source code.
3. Supporting Python 3 only will allow us to start using new features of the
language, ultimately making scikit-bio a better package. Python 3 contains new
and useful features which are not backported to Python 2, including function
annotations, type-hinting, the matrix multiplication operator, and asynchronous
programming support.
4. Given that scikit-bio is still in beta and has relatively few users (e.g.,
compared to a package like [numpy](http://www.numpy.org/)), making this switch
now would be better than waiting and could set a precedent for other scientific
Python packages. Switching sooner is also less of a developer burden because
the codebase is relatively small.

# Detailed design

## Timeline
We will drop Python 2 support in scikit-bio 0.5.0
([milestone](https://github.com/biocore/scikit-bio/milestones/0.5.0:%20Python%203%20support%20only)).
We will have at least
[one more release](https://github.com/biocore/scikit-bio/milestones/0.4.1) in
the 0.4.x series that includes Python 2 and 3 support.

## Implementation
Changes to support Python 3-only include:

1. Remove Python 2 tests from Travis-CI build matrix.
2. Remove `future` and `six` dependencies from scikit-bio.
3. Update documentation to indicate Python 3-only, including (but not limited
to) `README.md` and `CONTRIBUTING.md`. Make it very clear on the scikit-bio
front page that scikit-bio only supports Python 3 and link to resources
describing why.
4. Modify `doc/source/development/py3.rst` to be a guide for helping Python 2
developers become proficient Python 3 developers (e.g., documenting common
gotchas or new language features that apply to scikit-bio).
5. Update PyPI classifiers in `setup.py` to indicate Python 3-only.
6. Search codebase for references to Python 2 and remove/update accordingly.
This includes comments and actual code modifications.
7. Removes hacks in `skbio.io` necessary for Python 2 compatibility.
8. Test against latest release and nightly build of Python 3 on Travis-CI.
Allow failure against nightly Python 3 build.

Note: scikit-bio's doctests are already Python 3-only, so no updates necessary
there! :shipit:

# Drawbacks

Packages/tools supporting Python 2 (including those supporting Python 2
**and** 3) will not be able to depend on scikit-bio.

Excludes Python 2 user and developer communities.

# Alternatives

We could continue supporting Python 2 and 3 but this requires substantial
developer effort (see Motivation 2 above). Given our limited developer
bandwidth and free (as in :beer:) software model, supporting Python 2 and 3 is
not feasible.

Python 3 has been available for 7 years. This transition **must** happen at
some point because Python 3 is here to stay and Python 2 will eventually not
be supported by other systems and tools. By not making the switch now (or
soonish), we are delaying the inevitable and will affect more users as our
userbase increases in size. It would be good form to make this change now
while we are still in beta.

# Unresolved questions
