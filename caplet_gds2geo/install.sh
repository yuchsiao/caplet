#!/bin/bash

export GDSII=python-gdsii-0.2.1
export GDSII_EXT=.tar.gz

echo "caplet_gds2geo: Clean up gdsii folder and fetch python-gdsii ..."
rm -r gdsii
wget --no-check-certificate "https://pypi.python.org/packages/source/p/python-gdsii/$GDSII$GDSII_EXT"
echo "caplet_gds2geo: Decompress python-gdsii ..."
tar -xf $GDSII$GDSII_EXT
echo "caplet_gds2geo: Move gdsii folder and delete downloaded files ..."
mv $GDSII/gdsii ./
rm -r $GDSII $GDSII$GDSII_EXT
echo "caplet_gds2geo: Done."
