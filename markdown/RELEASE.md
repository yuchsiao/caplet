**1.1.0** - 2013-08-12

+ New: Improved algorithm for merging projection rectangles. Ver1.1 of `mergeProjection()` takes projection source distances into account.
+ New: Add function `removeBadProjection()` after `mergeProjection()` to remove bad projections that may increase the condition number.
+ New: Add options for `projectionDistance` and `projectionMergeDistance` for `caplet_geo` and `caplet_geo_cli`.
+ New: Add a button to `caplet_geo` for loading reference Cmat. Reference Cmat now are autoloaded when opening a `geo` file if there exists.
+ New: Add a compilation flag in `caplet_debug.h` for printing out system matrices and right-hand sides for `caplet_solver`.
+ New: Add an option for using double-precision LAPACK with instantiable basis functions.
* Fix: Fixed bugs in `ROBUST_INTEGRAL_CHECK`, `computeAdjacency()`, and `instantiableBasisFunction()`. Now fewer basis functions give the same level of accuracy. 

_________
**1.0.5** - 2013-08-07

+ New: Added default flag `ROBUST_INTEGRAL_CHECK` and modified default value of zero to 1e-12 in `caplet_solver/include/caplet_parameter.h` to improve integral quality.
* Fix: Improved integral quality to avoid sometimes nan or inf values.

_________
**1.0.4** - 2013-08-03

* Fix: Updated Makefile to be compatible with gcc that has linker flag '--as-needed' by default.
* Fix: Fixed bugs in `extend()` and `poly2rect()` of `caplet_geo/geoloader.cpp` that cause segmentation fault for certain cases.

_________
**1.0.3** - 2013-07-23

* Fix: Ensured output of `caplet_gds2geo.py` in integers so that `caplet_geo` reads correctly.

_________
**1.0.2** - 2013-07-20

+ New: Hosted a repo on GitHub.
+ New: Prepared this README.md file.
+ New: Tested installation with various Qt versions.
* Fix: Improved installation experience.
* Fix: Cleaned up several compilation warnings.

_________
**1.0.1** - 2013-07-12

* Fix: Improved installation by utilizing qmake for caplet_geo.
* Fix: Updated Makefile in caplet_solver that originally generates linking warnings.
* Fix: Updated caplet_gds2geo/install.sh so that gdsii library can be downloaded through proxies or behind firewalls.

__________
[**1.0.0** - 2013-02-15](http://sourceforge.net/projects/caplet/files/?source=navbar)

+ New: Repackaged caplet into three tools.
+ New: Enabled standard BEM solver in GUI (piecewise constant basis functions with collocation testing).
+ New: Simplified installation and uninstallation process.
+ New: Supported Command Line Interface for `caplet_geo` (CLI): `caplet_geo_cli`

__________
[**0.9.0** - 2013-01-31](http://sourceforge.net/projects/caplet/files/?source=navbar)
    
* New: Debut
