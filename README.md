[Caplet]
========

[Caplet] is a capacitance extraction and visualization toolkit for VLSI interconnects that

- Visualizes GDS layout files
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


Dependences
-----------


Installation
------------


Quickstarts
-----------


Tested Environments
-------------------



Release Notes
-------------
1.0.2 -2013-07-16

* Improve installation experience.
* Host a repo on GitHub.
* Clean up several compilation warnings.
* Prepare this README file.
______
[1.0.0 - 2013-02-15](http://sourceforge.net/projects/caplet/files/?source=navbar)

* Repackage caplet into three tools.
* Enable standard BEM solver in GUI (piecewise constant basis functions with collocation testing).
* Simplify installation and uninstallation process.
* Support Command Line Interface for `caplet_geo` (CLI): `caplet_geo_cli`

______
[0.9.0 - 2013-01-31](http://sourceforge.net/projects/caplet/files/?source=navbar)
    
* Debut




[Caplet]: http://www.rle.mit.edu/cpg/codes/caplet/
[Yu-Chung Hsiao]: yuchsiao@gmail.com
[author]: yuchsiao@gmail.com
[FASTCAP]: http://www.rle.mit.edu/cpg/codes/


  