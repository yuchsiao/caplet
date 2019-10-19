/*
Created : Jul 28, 2010
Author  : Yu-Chung Hsiao
Email   : project.caplet@gmail.com
*/

/*
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

#include "caplet.h"

#include "mpi.h" 

#include <iostream>
#include <string>
#include <sstream>
#include <list>
#include <cstdlib>
using namespace std;

 
void printUsage(char* command){
    cout << "CAPLET: Extract capacitance matrix " << endl
         << "        from instantiable basis functions (.caplet)" << endl
         << "        or   piecewise constant basis functions (.qui)" << endl
         << "Usage  : " << command << " [OPTION] INPUT.caplet [-o OUTPUT]" << endl
         << "   or  : " << command << " [OPTION] INPUT.qui    [-o OUTPUT]" << endl
         << "Option : " << endl
         << "  -d, --double              use double-precision LAPACK" << endl
         << "  -v, --version             print version info" << endl;
} 

void printVersion(){
    cout << "CAPLET Version 1.1" << endl;
}

int main(int argc, char *argv[]){
    using namespace caplet;

    string fileName;
    string folderPath = ".";
    string fileBaseName;
    string fileExtName;
    const string capletExt  = "caplet";
    const string fastcapExt = "qui";
    string fileNameCmat  = "";

    bool flagDouble = false; //* single precision fast solution

    list<string> argvList;
    for (int i=1; i<argc; ++i){ //* skip command name
        argvList.push_back(argv[i]);
    }

    //* Read options
    for ( list<string>::iterator each=argvList.begin();
          each!=argvList.end(); ){

        //* Read Cmat file name
        if (each->compare("-o")==0){
            each = argvList.erase(each);
            if (each == argvList.end()){
                printUsage(argv[0]);
                return 0;
            }
            fileNameCmat = *each;
            each = argvList.erase(each);
        }

        //* Flag -f for single-precision fast solution
        else if (each->compare("-d")==0 || each->compare("--double")==0 ){
            flagDouble = true;
            each = argvList.erase(each);
        }

        //* Flag -v --version
        else if (each->compare("-v")==0 || each->compare("--version")==0 ){
            printVersion();
            return 0;
        }

        else{
          ++each;
        }
    }

    //* Check if nonknown options
    for ( list<string>::iterator each=argvList.begin();
          each!=argvList.end(); ++each){

        if ( (*each)[0] =='-'){
            cout << "CAPLET: Unknown option (" << *each << ")." << endl;
            return 0;
        }
    }


    //* If not input file specified
    if ( argvList.empty()==true ){
        printUsage(argv[0]);
        return 0;
    }

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

    //* Call Caplet
    Caplet caplet;


    if ( fileExtName.compare(capletExt)==0 ){
        caplet.loadCapletFile(folderPath+"/"+fileName);
        if (flagDouble==true){
            caplet.extractC( Caplet::DOUBLE_GALERKIN );
        }
        else{
            caplet.extractC( Caplet::FAST_GALERKIN );
        }
    }
    else if ( fileExtName.compare(fastcapExt)==0 ){
        caplet.loadFastcapFile(folderPath+"/"+fileName);
        caplet.extractC( Caplet::DOUBLE_COLLOCATION );
    }
    else{
        //* unexpected file extension
        cout << "CAPLET: File type not supported (" << fileExtName << ")." << endl;
        return 0;
    }

    //* Save Cmat
    if (fileNameCmat.empty()==false){
        caplet.saveCmat(fileNameCmat);
    }

    return 0; 
} 

