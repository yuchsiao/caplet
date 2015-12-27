INTRODUCTION='[Caplet]
========

[Caplet] is a capacitance extraction and visualization toolkit for VLSI interconnects that

- Visualizes GDS2 layout files
- Generates basis functions automatically
- Extracts parasitic capacitance matrices in a parallel manner 
- Provides field-solver accurate solutions

**If you find this software useful in your personal or commercial purpose, please encourage the author by mentioning your project at <project.caplet@gmail.com>.**
'

FEATURE='Features
--------
1. Capable of transforming GDS2 layout files
2. Aimed for ultra-fast extraction of small-to-medium structures
3. Accurate within 5% of reference solutions (using standard boundary element method on finely discretized geometries)
4. Faster than [FASTCAP] at single-core execution for the same accuracy (5X faster for a NAND gate)
5. Efficiently parallelized (90% efficiency for 24x24 buses at eight-core execution)
6. Able to automatically generate piecewise constant and instantiable basis functions 
7. GUI and OpenGL visualization
'

SCREENSHOT='Screenshots
-----------
![GUI](http://www.mit.edu/~yuchsiao/caplet/img/gui_piecewise_constant_basis_function.png "GUI and Piecewise constant basis functions")
'

LICENSE='License
-------
GNU Lesser General Public License Version 3 (GPLv3)

Caplet is distributed in the hope that it will be useful, but WITHOUT ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License for more details.
'

LINK='
[Caplet]: http://www.rle.mit.edu/cpg/codes/caplet/
[FASTCAP]: http://www.rle.mit.edu/cpg/research_codes.htm
'

REQUIREMENT='Requirements
------------
'

INSTALLATION='Installation
------------
'

QUICKSTART='Quickstarts
-----------
'

TESTENV='Tested Environments
-------------------
'

RELEASE='Release Notes
-------------
'

FORMAT='FORMAT
------
'

LIMITATION='Limitation
----------
'

echo "$INTRODUCTION"
echo

echo "$FEATURE"
echo

echo "$SCREENSHOT"
echo 

echo "$LICENSE"
echo 

echo "$REQUIREMENT"
cat REQUIREMENT.md
echo 
echo

echo "$INSTALLATION"
cat INSTALLATION.md
echo
echo

echo "$TESTENV"
cat TESTENV.md
echo
echo

echo "$QUICKSTART"
cat QUICKSTART.md
echo
echo

echo "$FORMAT"
cat FORMAT.md
echo
echo

echo "$LIMITATION"
cat LIMITATION.md
echo
echo

echo "$RELEASE"
cat RELEASE.md
echo
echo

echo "$LINK"

