[Caplet]
========

[Caplet] is a capacitance extraction and visualization toolkit for VLSI interconnects that

- Visualizes GDS2 layout files
- Generates basis functions automatically
- Extracts parasitic capacitance matrices in a parallel manner 
- Provides field-solver accurate solutions

**If you find this software useful in your personal or commercial purpose, please encourage the author by mentioning your project to the [author].**

Copyright (C) 2013 [Yu-Chung Hsiao]


Features
--------
1. Capable of transforming GDS2 layout files
2. Aimed for ultra-fast extraction of small-to-medium structures
3. Accurate within 5% of reference solutions (using standard boundary element method on finely discretized geometries)
4. Faster than [FASTCAP] at single-core execution for the same accuracy (5X faster for a NAND gate)
5. Efficiently parallelized (90% efficiency for 24x24 buses at eight-core execution)
6. Able to automatically generate piecewise constant and instantiable basis functions 
7. GUI and OpenGL visualization


Screenshots
-----------
![GUI](http://www.mit.edu/~yuchsiao/caplet/img/gui_piecewise_constant_basis_function.png "GUI and Piecewise constant basis functions")

Version
-------
1.0.2


License
-------
GNU Lesser General Public License Version 3 (GPLv3)

Caplet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


Prerequisites
-------------

* g++ 4.4 or higher

`caplet_gds2geo`

* Python 2.6 or higher
* [GDSII for Python](http://gdspy.sourceforge.net/): automatically downloaded with `install.sh`

`caplet_geo`

* Qt 4.5 or higher
* OpenGL

`caplet_solver`

* openmpi, openmp
* gfortran, lapack, blas


Installation
------------

To install, execute `install.sh`. Make sure `qmake` is in your system path.

To uninstall, execute `uninstall.sh`.

If you encounter any problem, go to the individual folder and run `make`.
If problems happen at `caplet_geo`, run `qmake` then `make`.


Quickstarts
-----------

`caplet` consists of three parts: **`caplet_gds2geo`**, **`caplet_geo`**, and **`caplet_solver`**.

**`caplet_gds2geo`** transforms binary GDSii layout files into ascii geometry definitions, including squares and polygons. The generated geometry files end with `.geo` extension. This program needs an additional file to specify the elevation of each metal layer and connection relationship between layers and vias. The usage is as the following:

```
python caplet_gds2geo.py -l LAYER_FILE GDS2_FILE
```

**Example** (under folder `caplet_gds2geo`):

```
python caplet_gds2geo.py -l sample.tech cap_inverter.gds
```

**`caplet_geo`** decomposes 2D polygons into non-overlapping 3D rectangles, and generate piecewise constant (PWC) basis functions or instantiable basis functions of your choice. The usage should be straightforward: open a .geo file, select the type of basis function type and parameters for your purpose, and click on **Extract** to extract the capacitance matrix using `caplet_solver`. `caplet_geo` also provides iterative schemes for calculating the finely discreted PWC reference capacitance matrices for accuracy comparison.

Similar to `caplet_geo`, The Command Line Interface (CLI) version `caplet_geo_cli` also generates either type of basis functions but does not provide visualization. The command line usage is the following:

```
caplet_geo_cli [--type pwc or ins] [--unit n or u or m or 1] [--size value] filename.geo
```

`--type` is followed by either `pwc` for piecewise constant basis functions or by `ins` for instantiable basis functions. The default is `--type ins`.

`--unit` specifies the unit for parameters following `--size`. `n` for nanometer, `u` for micrometer, `m` for millimeter, and `1` for meter. The default is `--unit n`.

`--size` sets up the panel size (one size) for constructing piecewise constant basis functions or the arch length for instantiable basis functions. The default value for `--type pwc` is 50 and for `--type ins` is 300.

When `--type pwc`, the output file format and extension is `.qui`, which is compatible with [FASTCAP]. When `--type ins`, the output file ends with extension `.caplet`.

**Example** (under folder `caplet_geo`)

```
./caplet_geo_cli example/cap_inverter.geo 
```

generates instantiable basis functions with 300nm arch length.

```
./caplet_geo_cli --type pwc --unit n --size 100 example/cap_inverter.geo
```

generates piecewise constant basis functions with panel size 100nm.

**`caplet_solver`** extracts capacitance matrices from `.qui` files which list PWC basis functions or from `.caplet` file which list instantiable basis functions for all conductors.



Examples
--------



Tested Environments
-------------------

* platform: lmde-64bit (linux mint debian edition), mint15-64bit (ubuntu13.04), ubuntu10.04-64bit
* g++: 4.7.3, 4.4.3
* python: 2.7.4, 2.6.5 
* Qt: 5.1.0, 5.0.2, 4.8.5, 4.7.4, 4.6.4, 4.5.3, 4.4.3, 4.3.5 (not compatible)


Release Notes
-------------
**1.0.2** - 2013-07-16

* Improve installation experience.
* Host a repo on GitHub.
* Clean up several compilation warnings.
* Prepare this README file.

______
[**1.0.0** - 2013-02-15](http://sourceforge.net/projects/caplet/files/?source=navbar)

* Repackage caplet into three tools.
* Enable standard BEM solver in GUI (piecewise constant basis functions with collocation testing).
* Simplify installation and uninstallation process.
* Support Command Line Interface for `caplet_geo` (CLI): `caplet_geo_cli`

______
[**0.9.0** - 2013-01-31](http://sourceforge.net/projects/caplet/files/?source=navbar)
    
* Debut




[Caplet]: http://www.rle.mit.edu/cpg/codes/caplet/
[Yu-Chung Hsiao]: yuchsiao@gmail.com
[author]: yuchsiao@gmail.com
[FASTCAP]: http://www.rle.mit.edu/cpg/research_codes.htm


