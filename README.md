[Caplet]
========

[Caplet] is a capacitance extraction and visualization toolkit for VLSI interconnects that

- Visualizes GDS2 layout files
- Generates basis functions automatically
- Extracts parasitic capacitance matrices in a parallel manner 
- Provides field-solver accurate solutions

**If you find this software useful in your personal or commercial purpose, please encourage the author of this work by mentioning your project to <project.caplet@gmail.com>.**


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
![GUI](http://www.feature.space/caplet/img/gui_color_by_layer.png "GUI and Piecewise constant basis functions")


License
-------
GNU Lesser General Public License Version 3 (GPLv3)

Caplet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.


Requirements
------------

#### General

- gcc 4.4 or higher

#### `caplet_gds2geo`

- Python 2.6 or higher
- [GDSII for Python](http://gdspy.sourceforge.net/): automatically downloaded with `install.sh`

#### `caplet_geo`

- Qt 4.4 or higher
- OpenGL

#### `caplet_solver`

- openmpi, openmp
- gfortran, lapack, blas


Installation
------------

To install with default settings, 

1. Execute `install.sh`. Make sure `qmake` is in your system path.
2. Put a symbolic link or binary of FASTCAP in the same folder of `caplet_geo` for computing reference capacitance matrices.

To modify settings, see the file `caplet_solver/include/caplet_parameter.h` for `caplet_solver`.

To uninstall, execute `uninstall.sh`.

If you encounter any problem, go to the individual folder and then `make`.
If problems happen at `caplet_geo`, run `qmake` first and then `make`.


Tested Environments
-------------------

- **platform:** lmde-64bit (linux mint debian edition), mint15-64bit (ubuntu13.04), ubuntu12.04-64bit, ubuntu10.04-64bit
- **g++:** 4.7.3, 4.4.3
- **python:** 2.7.4, 2.6.5 
- **Qt:** 5.1.0, 5.0.2, 4.8.5, 4.7.4, 4.6.4, 4.5.3, 4.4.3, 4.3.5 (not compatible)
- **OpenMPI:** 1.4.5, 1.4.1
- **python-gdsii:** 0.2.1


Quickstarts
-----------

#### `caplet`
`caplet` consists of three parts: `caplet_gds2geo`, `caplet_geo`, and `caplet_solver`, each of which, respectively, processes layout `.gds` files into geometry `.geo` files, generates basis function files `.qui` or `.caplet`, and extracts and prints performance and capacitance matrices on screen. Each program takes the output file of the previous stage as input.

#### `caplelt_gds2geo`
`caplet_gds2geo` transforms binary GDSii layout files into ascii geometry definitions, including squares and polygons. The generated geometry files end with `.geo` extension. This program needs an additional file to specify the elevation of each metal layer and connection relationship between layers and vias. The usage is as the following:

```
python caplet_gds2geo.py -l LAYER_FILE GDS2_FILE
```

**Example** (under folder `caplet_gds2geo`)

```
python caplet_gds2geo.py -l sample.tech cap_inverter.gds
```
 
#### `caplet_geo`
`caplet_geo` decomposes 2D polygons into non-overlapping 3D rectangles, and generate piecewise constant (PWC) basis functions or instantiable basis functions of your choice. The usage should be straightforward: open a .geo file, select the type of basis function type and parameters for your purpose, and click on **Extract** to extract the capacitance matrix using `caplet_solver`. `caplet_geo` also provides iterative schemes for calculating the finely discreted PWC reference capacitance matrices for accuracy comparison.

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

####`caplet_solver`
`caplet_solver` extracts capacitance matrices from `.qui` files which list PWC basis functions or from `.caplet` files which list instantiable basis functions for all conductors. Two binary executables `capletMPI` and `capletOpenMP` are generated after compilation. As suggested by their names, `capletMPI` is the capacitance extraction solver parallelized by MPI, and `capletOpenMP` is parallelized by OpenMP. The usage of `capletOpenMP` is as the following:

For single thread extraction, simply run `capletMPI` as common executables:

```
capletMPI filename.ext
```

When the extension is `.caplet`, `capletMPI` extracts capacitance matrices using instantiable basis functions. When the extension is `.qui`, `capletMPI` solves the problem by the standard boundary element method with piecewise constant basis functions.

The number of processors `N` can be specified through `mpirun` (assuming that `mpirun` is in the system path:

```
mpirun -np N capletMPI filename.caplet
```

Such parallelization is only implemented for instantiable basis functions. There is no parallelization effect for `.qui` files.

The usage of `capletOpenMP` is similar:

```
capletOpenMP filename.ext
```

The number of threads in `capletOpenMP` is fixed after compilation. The parameter is defined as `CAPLET_OPENMP_NUM_THREADS` in `caplet_solver/include/caplet_parameter.h`. Once `caplet_parameter.h` is modified, recompilation is required.

**Example** (under folder `caplet_solver`)
Use four cores to extract capacitance out of instantiable basis functions and save the result in `result` (given `mpirun` is in the system path)

```
mpirun -np 4 bin/capletMPI example/cap_inverter_300nm.caplet > result
```

Use a single core to extract capacitance out of piecewise constant basis functions and print the result on screen

```
bin/capletMPI example/cap_inverter_200nm.qui
```


FORMAT
------

####`.geo`

Geometry files `.geo` define metal layer elevations (z-direction), via elevations and which metal layers are connected from above and below, and the polygon descriptions for each layer. Currently, only boundary and path types of GDS2 are processed into polygons. All other types of GDS2 are ignored. 

The format is as the following


    #metal_layers
    layer_index, z_min, z_max 
    .
    .
    #vias
    layer_index, z_min, z_max, bottom_layer_index, top_layer_index
    layer_index
    {#polygons_in_a_layer
    [#vertices_in_a_polygon
    vertex_1
    vertex_2
    .
    .
    vertex_last]}
    .
    .
    

It is worth noting that metal layers and vias are indexed together, staring from 0. The square bracket part [] describes x- and y-coordinates of vertices for a polygon. There are total `#polygons_in_a_layer` such square brackets for a layer. Once all polygons in a single layer are described, we move to the next layer description enclosed by curly brackets {}. The file ends when all layers, including metal layers and vias, and all polygons in each layer are listed. The brackets are not written in output files explicitly. We put them here only for explanation purposes.

All numbers are integers with an implicit unit nm. Check if your GDS2 files are defined in different units.

Examples can be found under folder `caplet_geo/example`.



####`.caplet`

Instantiable basis function files `.caplet` describe fundamental shapes that are used to constitute instantiable basis functions. The format is as the following

    #conductors
    #shapes_for_conductor_1 #shapes_for_conductor_2 ...
    total_number_of_shapes
    [ST Inc XL XU YL YU ZL ZU Dir SDir Decay Ignore]
    .
    .

The square bracket records information for one shape. The file terminates when all shapes of `total_number_of_shapes` are listed. Each element is explained below:

`ST`: Shape Type is either 'A' for arch or 'F' for flat.

`Inc`: Increment value is either '1' or '0'. A shape with increment '0' belongs to the same basis function of the closest previous shape with '1'.

`XL` to `ZU`: Lower and Upper x-, y-, and z-coordinates.

`Dir`: Surface normal direction.

`SDir`: Shape varying direction.

`Decay`: Positive sign for shape decaying toward the positive direction, and negative sign for decaying toward the negative direction.

`Ignore`: Reserved and currently not used in `caplet_solver`.

Examples can be found under folder `caplet_solver/example`.



####`.qui`

Files are originally used for FASTCAP. `caplet_geo` can also load and visualize them in 3D. The format is as the following:

    0
    shape conductor_label x1 y1 z1 x2 y2 z2 ...
    .
    .

The first line led by 0 is always ignored. From the second line to the last, each line defines a flat shape or a piecewise constant basis function. 

`shape`: `Q` for quadrilateral and `T` for triangular shapes.
`conductor_label`: arbitrary single word label. The common choice is integer indices.
`xn yn zn`: x-, y-, and z-coordinates for a single vertex. For `Q`, four such triplets follow `conductor_label`. For `T`, there are three triplets instead. The vertices should be recorded either clockwise or counter-clockwise.
In Manhattan geometries, we only use `Q` for piecewise constant basis functions.

More information can be found in [FASTCAP](http://www.rle.mit.edu/cpg/research_codes.htm) manual (`ug.tex` under `doc`).


Limitation
----------

#### `caplet_gds2geo`

- GDS files are assumed flattened.

#### `caplet_geo`

- Metal layers and vias are disjoint in the vertical direction.
- Sublayers are not taken into account.

#### `caplet_solver`

- Relative dielectric constant is fixed to be 1.
- Dielectric material is assumed to be uniform.


Release Notes
-------------

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

+ New: Added default flag ROBUST_INTEGRAL_CHECK and modified default value of zero to 1e-12 in `caplet_solver/include/caplet_parameter.h` to improve integral quality.
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



[Caplet]: http://www.rle.mit.edu/cpg/codes/caplet/
[FASTCAP]: http://www.rle.mit.edu/cpg/research_codes.htm

