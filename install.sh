#!/bin/bash

echo -e "\e[1;36mCAPLET:\e[m \e[1mMaking caplet_gds2geo ...\e[m"
cd caplet_gds2geo
make

echo -e "\e[1;36mCAPLET:\e[m \e[1mMaking caplet_geo and caplet_geo_cli ...\e[m"
cd ../caplet_geo
make
make -f MakefileCLI

echo -e "\e[1;36mCAPLET:\e[m \e[1mMaking caplet_solver ...\e[m"
cd ../caplet_solver
make
cd ..
echo -e "\e[1;36mCAPLET:\e[m \e[1;35mDone\e[1m installation.\e[m"
