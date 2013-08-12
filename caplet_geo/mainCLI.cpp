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

void printUsage(const string &command){
    cout << "Usage  : " << command << " [OPTIONS] filename.geo" << endl;
    cout << "Output : filename.qui    for piecewise constant basis functions" << endl
         << "         filename.caplet for instantiable basis functions" << endl;
    cout << "Options: " << endl
         << "       Basis Function Type:" << endl
         << "       -t,--type pwc  : piecewise constant basis functions" << endl
         << "       -t,--type ins  : instantaible basis functions  (default)" << endl
         << endl
         // << "       Grid Unit for GDS and Size:" << endl
         // << "       --unit n    : unit in nm (default)" << endl
         // << "       --unit u    : unit in um" << endl
         // << "       --unit m    : unit in mm" << endl
         // << "       --unit 1    : unit in  m" << endl
         // << endl
         << "       Basis Function Parameter:" << endl
         << "       -s,--size        value: panel size for piecewise constant basis functions " << endl
         << "                               (default:  50e-9, value>0)" << endl
         << "                               arch length for instantiable basis functions" << endl
         << "                               (default: 300e-9, value>0: normal, " << endl
         << "                                                 value=0: no arch, " << endl
         << "                                                 value<0: no flat shapes" << endl
         << "       -p,--proj-dist   value: projection distance (default: 2e-6)" << endl
         << "       -m,--merge-dist  value: projection merge distance (default: 1e-7)" << endl
         << endl;    
}

int main(int argc, char *argv[])
{
    if (argc<2){
        printUsage(argv[0]);
        return 0;
    }

    //* Convert argv to list
    list<string> argvList;
    for (int i=1; i<argc; ++i){ //* from 1: no command name
        argvList.push_back(argv[i]);
    }

    //* Read options
    BasisFunctionType basisFunctionType = INSTANTIABLE_BASIS;
    float unit = 1e-9;
    float size = 300;
    bool   isSizeInput = false;

    float projDist  = 2000 *unit;
    float mergeDist =   10 *unit;

    for ( list<string>::iterator each=argvList.begin();
          each!=argvList.end(); ){

        //* -t,--type
        if (each->compare("--type")==0 || each->compare("-t")==0){
            if (argvList.empty()==true){
                printUsage(argv[0]);
                return 0;
            }
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

        //* -u,--unit
        // if (each->compare("--unit")==0 || each->compare("-u")==0){
        //     each = argvList.erase(each);
        //     switch (*each->begin()){
        //     case 'u':
        //         unit = 1e-6;
        //         break;
        //     case 'n':
        //         unit = 1e-9;
        //         break;
        //     case 'm':
        //         unit = 1e-3;
        //         break;
        //     case '1':
        //         unit = 1;
        //         break;
        //     default:
        //         unit = 1e-9;
        //     }
        //     each = argvList.erase(each);
        //     continue;
        // }

        //* -s,--size
        if (each->compare("--size")==0 || each->compare("-s")==0){
            if (argvList.empty()==true) {
                printUsage(argv[0]);
                return 0;
            }
            each = argvList.erase(each);
            isSizeInput = true;
            string sizeString = *each;
            istringstream sizeSS(sizeString);
            sizeSS >> size;
            each = argvList.erase(each);
            continue;
        }

        //* -p,--proj-dist
        if (each->compare("--proj-dist")==0 || each->compare("-p")==0){
            if (argvList.empty()==true) {
                printUsage(argv[0]);
                return 0;
            }
            each = argvList.erase(each);
            string projDistString = *each;
            istringstream projDistSS(projDistString);
            projDistSS >> projDist;
            each = argvList.erase(each);
            continue;
        }

        //* -m,--merge-dist
        if (each->compare("--merge-dist")==0 || each->compare("-m")==0){
            if (argvList.empty()==true) {
                printUsage(argv[0]);
                return 0;
            }
            each = argvList.erase(each);
            string mergeDistString = *each;
            istringstream mergeDistSS(mergeDistString);
            mergeDistSS >> mergeDist;
            each = argvList.erase(each);
            continue;
        }

        //* increment
        ++each;
    }


    //* Check if nonknown options
    for ( list<string>::iterator each=argvList.begin();
          each!=argvList.end(); ++each){

        if ( (*each)[0] =='-'){
            cout << "CAPLET_GEO_CLI: Unknown option (" << *each << ")." << endl;
            return 0;
        }
    }


    //* If not input file specified
    if ( argvList.empty()==true ){
        cout << "CAPLET_GEO_CLI: No input file specified." << endl;
        printUsage(argv[0]);
        return 0;
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
            writeCapletFile(fileBaseName, geoloader.getInstantiableBasisFunction(unit, size*unit, projDist, mergeDist));
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
    cout << "CAPLET_GEO: Done basis functions construction. (" << outputFileName << ")" << endl;

    return 0;
}

