/*
Created : Jul 28, 2010
Modified: Feb 14, 2013
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

 
int main(int argc, char *argv[]){
	using namespace caplet;

	if (argc<2){
        cout << "CAPLET: Extract capacitance matrix " << endl
             << "        from instantiable basis functions (.caplet)" << endl
             << "        or   piecewise constant basis functions (.qui)" << endl;
        cout << "Usage : " << argv[0] << " filename.caplet" << endl;
        cout << "   or : " << argv[0] << " filename.fastcap" << endl;
		return 0;
	}
    string fileName;
    string folderPath = ".";
    string fileBaseName;
    string fileExtName;
    const string capletExt  = "caplet";
    const string fastcapExt = "qui";
    string fileNameCmat  = "";

	list<string> argvList;
	for (int i=1; i<argc; ++i){
		argvList.push_back(argv[i]);
	}
	
    //* Read the output file name
	for ( list<string>::iterator each=argvList.begin();
		  each!=argvList.end(); ){
		if (each->compare("-o")==0){
			each = argvList.erase(each);
			fileNameCmat = *each;
			each = argvList.erase(each);
			break;
		}
		++each;
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
        caplet.extractC( Caplet::FAST_GALERKIN );
    }
    else if ( fileExtName.compare(fastcapExt)==0 ){
        caplet.loadFastcapFile(folderPath+"/"+fileName);
        caplet.extractC( Caplet::DOUBLE_COLLOCATION );
    }
    else{
        //* unexpected file extension
        cout << "CAPLET: Not supported input file type (" << fileExtName << ")." << endl;
        exit(0);
    }

    //* Save Cmat
    if (fileNameCmat.empty()==false){
        caplet.saveCmat(fileNameCmat);
    }

    return 0; 
} 

