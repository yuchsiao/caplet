/*
CREATED : Feb 15, 2013
MODIFIED:
AUTHOR  : Yu-Chung Hsiao
EMAIL   : project.caplet@gmail.com

This file is part of CAPLET.

CAPLET is free software: you can redistribute it and/or modify
it under the terms of the GNU Lesser General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

CAPLET is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU Lesser General Public License for more details.

You should have received a copy of the Lesser GNU General Public License
along with CAPLET.  If not, see <http://www.gnu.org/licenses/>.
*/

#include "geoloader.h"

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <cstdlib>
using namespace std;

enum  BasisFunctionType{ PWC_BASIS, INSTANTIABLE_BASIS };
const string geoExt     = "geo";
const string fastcapExt = "qui";
const string capletExt  = "caplet";

int main(int argc, char *argv[])
{
    if (argc<2){
        cout << "CAPLET_GEO_CLI: Construct basis functions of specified type " << endl;
        cout << "Usage  : " << argv[0] << " --options filename.geo" << endl;
        cout << "Output : filename.qui    for piecewise constant basis functions" << endl
             << "         filename.caplet for instantiable basis functions" << endl;
        cout << "Options: " << endl
             << "       Basis Function Type:" << endl
             << "       --type pwc  : piecewise constant basis functions" << endl
             << "       --type ins  : instantaible basis functions  (default)" << endl
             << endl
             << "       Grid Unit for GDS and Size:" << endl
             << "       --unit n    : unit in nm (default)" << endl
             << "       --unit u    : unit in um" << endl
             << "       --unit m    : unit in mm" << endl
             << "       --unit 1    : unit in  m" << endl
             << endl
             << "       Basis Function Parameter:" << endl
             << "       --size value: panel size for piecewise constant basis functions " << endl
             << "                     (default:  50, value>0)" << endl
             << "                     arch length for instantiable basis functions" << endl
             << "                     (default: 300, value>0: normal, " << endl
             << "                                    value=0: no arch, " << endl
             << "                                    value<0: no flat shapes" << endl
             << endl;
        return 0;
    }

    //* Convert argv to list
    list<string> argvList;
    for (int i=1; i<argc; ++i){ //* from 1: no command name
        argvList.push_back(argv[i]);
    }

    //* Read options
    BasisFunctionType basisFunctionType = INSTANTIABLE_BASIS;
    double unit = 1e-9;
    double size = 300;
    bool   isSizeInput = false;
    for ( list<string>::iterator each=argvList.begin();
          each!=argvList.end(); ){

        //* --type
        if (each->compare("--type")==0){
            each = argvList.erase(each);
            string type = *each;
            if (type.compare("pwc")==0){
                basisFunctionType = PWC_BASIS;
            }
            else if(type.compare("ins")==0){
                basisFunctionType = INSTANTIABLE_BASIS;
            }
            else{ //* unexpected basis function type
                cout << "CAPLET_GEO: Unknown basis function type (" << type << ")" << endl;
                exit(0);
            }
            each = argvList.erase(each);
            continue;
        }

        //* --unit
        if (each->compare("--unit")==0){
            each = argvList.erase(each);
            switch (*each->begin()){
            case 'u':
                unit = 1e-6;
                break;
            case 'n':
                unit = 1e-9;
                break;
            case 'm':
                unit = 1e-3;
                break;
            case '1':
                unit = 1;
                break;
            default:
                unit = 1e-9;
            }
            each = argvList.erase(each);
            continue;
        }

        //* --size
        if (each->compare("--size")==0){
            each = argvList.erase(each);
            isSizeInput = true;
            string sizeString = *each;
            istringstream sizeSS(sizeString);
            sizeSS >> size;
            each = argvList.erase(each);
            continue;
        }

        //* increment
        ++each;
    }

    //* Set size to default value if size is not specified.
    if (isSizeInput==false){
        switch(basisFunctionType){
        case PWC_BASIS:
            size = 50;
            break;
        case INSTANTIABLE_BASIS:
            size = 300;
            break;
        default:
            ;
        }
    }
    else {
        if (basisFunctionType == PWC_BASIS && size <= 0){
            cout << "CAPLET_GEO: PWC panel size has to be positive." << endl;
            exit(0);
        }
    }

    //* Read file name
    string fileName;
    string folderPath = ".";
    string fileBaseName;
    string fileExtName;

    //* Get folder path
    fileName = argvList.front();
    unsigned index = fileName.find_last_of('/');
    if (index<fileName.size()){
        //* '/' found
        folderPath = fileName.substr(0, index);
        fileName = fileName.substr(index+1);
    }

    //* Get file name and extension
    index = fileName.find_last_of('.');
    if (index<fileName.size() && index!=0){
        //* '.' found and not for denoting a hidden file
        fileBaseName = fileName.substr(0, index);
        fileExtName = fileName.substr(index+1);
    }

    //* Check ext name
    if ( fileExtName.compare(geoExt)!=0 ){
        //* unexpected file extension
        cout << "CAPLET_GEO: Not supported file type. (" << fileExtName << ")" << endl;
        exit(0);
    }


    //* Setup GeoLoader
    GeoLoader geoloader;
    try{
        geoloader.loadGeo(folderPath + "/" + fileName);
    }
    catch (FileNotFoundError e){
        cerr << "ERROR: File not found. (" << (folderPath+"/"+fileName) << ")" << endl;
        exit(1);
    }
    catch (GeometryNotManhattanError e){
        cerr << "ERROR: " << e.what() << endl;
        exit(1);
    }

    //* Construct basis functions
    string outputFileName = fileBaseName+".";
    switch(basisFunctionType){
    case PWC_BASIS:{
        outputFileName += fastcapExt;
        try{
            writeFastcapFile(fileBaseName, geoloader.getPWCBasisFunction(unit, size*unit));
        }
        catch (FileNotFoundError e){
            cerr << "ERROR: Cannot write file. (" << outputFileName << ")" << endl;
            exit(1);
        }
        const ConductorFPList   condList = geoloader.getPWCBasisFunction(unit, size*unit);
        cout << condList.front().size() << endl;
    }break;
    case INSTANTIABLE_BASIS:{
        outputFileName += capletExt;
        try{
            writeCapletFile(fileBaseName, geoloader.getInstantiableBasisFunction(unit, size*unit));
        }
        catch (FileNotFoundError e){
            cerr << "ERROR: Cannot write file. (" << outputFileName << ")" << endl;
            exit(1);
        }
    }break;
    default:
        cerr << "ERROR: Unknown basis function type." << endl;
        exit(1);
    }
    cout << "CAPLET_GEO: Basis functions constructed! (" << outputFileName << ")" << endl;

    return 0;
}

