#!/bin/bash

echo -e "\e[1;36mCAPLET:\e[m \e[1mClean and remove caplet_gds2geo ...\e[m"
cd caplet_gds2geo
make clean
rm -f caplet_gds2geo
cd ..

echo -e "\e[1;36mCAPLET:\e[m \e[1mClean and remove caplet_geo and caplet_geo_cli ...\e[m"
cd caplet_geo
qmake
make clean
make -f MakefileCLI clean
rm -f Makefile
rm -f caplet_geo caplet_geo_cli
cd ..

echo -e "\e[1;36mCAPLET:\e[m \e[1mClean and remove caplet_solver ...\e[m"
cd caplet_solver
make clean
rm -f bin/capletMPI bin/capletOpenMP
cd ..
echo -e "\e[1;36mCAPLET:\e[m \e[1;35mDone uninstallation.\e[m"
