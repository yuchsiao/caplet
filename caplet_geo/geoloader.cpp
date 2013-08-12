/*
CREATED : Jan 31, 2013
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

#include "debug.h"
#include <list>
#include <string>
#include <vector>
#include <utility>
#include <fstream>
#include <algorithm>
#include <limits>
#include <map>
#include <cmath>
#include <cstdlib>

#include <fstream>
#include <iomanip>

#include <ctime>

using namespace std;

//**
//* debugging utitlity functions
template<class T>
void printMatrix(T** matrix, int nrow, int ncol){
    using namespace std;
    for (int i=0; i<nrow; i++){
        cout << i << " \t[ ";
        for (int j=0; j<ncol; j++){
            cout << matrix[i][j] << ", \t";
        }
        cout << "] " << endl;
    }
}



//****
//*
//* GeoLoader
//*
//*

//**
//* GeoLoader constructor
GeoLoader::GeoLoader()
    :isLoaded(false){
}

//**
//* GeoLoader destructor
GeoLoader::~GeoLoader(){
    clear();
}


//**
//* loadGeo
//* - Test if geoFile exists first before clear things up
void GeoLoader::loadGeo(const string &geoFile) throw (FileNotFoundError, GeometryNotManhattanError){

    //* read geomery definitions from geoFile
    LayeredPolygonList metalLayeredPolygonList;
    LayeredPolygonList viaLayeredPolygonList;
    try{
        readGeo(geoFile, metalLayeredPolygonList, viaLayeredPolygonList); // may throw FileNotFoundError
    }
    catch (FileNotFoundError &e){
        throw;
    }

    //* check if Manhattan geometries
    for (unsigned int i=0; i<metalLayeredPolygonList.size(); ++i){
        for ( PolygonList::const_iterator eachPolyIt = metalLayeredPolygonList[i].begin();
              eachPolyIt != metalLayeredPolygonList[i].end(); ++eachPolyIt){
            if ( eachPolyIt->isManhattan() == false ){
                throw GeometryNotManhattanError();
            }
        }
    }
    //* poly to rect for metal
    LayeredRectangleList metalLayeredRectangleList(nMetal);
    for ( int i=0; i<nMetal; ++i ){
        poly2rect(metalLayeredPolygonList[i], metalLayeredRectangleList[i]);
        metalLayeredPolygonList[i].clear();
    }
    //* poly to rect for via
    viaLayeredRectangleList.resize(nVia);
    for ( int i=0; i<nVia; ++i ){
        poly2rect(viaLayeredPolygonList[i], viaLayeredRectangleList[i]);
        viaLayeredPolygonList[i].clear();
    }

    //* Group connected 2D x-y plane rects into a list (metal only)
    LayeredConnectedRectangleList metalLayeredConnectedRectangleList(nMetal);
    for ( int i=0; i<nMetal; ++i ){
        generateConnectedRects(metalLayeredRectangleList[i], metalLayeredConnectedRectangleList[i]);
        metalLayeredRectangleList[i].clear();
    }


    //* decompose x-y plane rects in each connected conductor
    //* and then generate 3d rects
    for ( int i=0; i<nMetal; ++i ){
        //* each metal layer
        for ( ConnectedRectangleList::iterator eachListIt
              = metalLayeredConnectedRectangleList[i].begin();
              eachListIt != metalLayeredConnectedRectangleList[i].end();
              ++eachListIt )
        {
            //* each connected 2D rectangle list in a metal layer

            //* decompose x-y plane rects if they overlap
            eachListIt->decompose();
            //* combine x-y plane rects if adjacent ones can form a single rect
            eachListIt->merge();
            //* compute adjacency for each 2D rect
            DirAdjacencyListOfRectangleList adjacency;
            DirAdjacencyListOfRectangleList compAdjacency;
            computeAdjacency(*eachListIt, adjacency, compAdjacency);
            //* construct 3D rects
            metalConductorList.push_back(Conductor(nMetal, nVia));
            generate3dRects(*eachListIt, compAdjacency, metalDef[i][0], metalDef[i][1], i, metalConductorList.back());
            //* merge 3D rects
            for (int dirIndex=0; dirIndex<Conductor::nDir; ++dirIndex){
                metalConductorList.back().layer[i][dirIndex].merge();
            }

            eachListIt->clear();
        }
    }


    //* UNCOMMENT to generate matlab structure output
    //printConductorListMatlab(conductorList);
    isLoaded = true;
}

//**
//* GeoLoader::getGeometryConductorList
//* -
const ConductorFPList &GeoLoader::getGeometryConductorList(const float unit)
{
    generateConductorList(geometryConductorList, false);
    geometryConductorFPList.constructFrom(geometryConductorList, unit);
    geometryConductorList.clear();
    return geometryConductorFPList;
}

const ConductorFPList &GeoLoader::getPWCBasisFunction(const float unit, const float suggestedPanelSize)
{
    clock_t tBefore = clock();
    generateConductorList(geometryConductorList, true);
    pwcConductorFPList.constructFrom(geometryConductorList, unit);
    geometryConductorList.clear();
    discretizeDisjointSurface(pwcConductorFPList, suggestedPanelSize);
    clock_t tAfter = clock();

    tPWCConstruction = difftime(tAfter, tBefore)/CLOCKS_PER_SEC;

    return pwcConductorFPList;
}

const ConductorFPList &GeoLoader::getPWCBasisFunction() const
{
    return pwcConductorFPList;
}


const ConductorFPList &GeoLoader::getInstantiableBasisFunction(const float unit, const float archLength,
                                                               const float projectionDistance, const float projectionMergeDistance){

    clock_t tBefore = clock();

    generateConductorList(geometryConductorList, false);
    instantiableConductorFPList.constructFrom(geometryConductorList, unit);
    geometryConductorList.clear();
    instantiateBasisFunction(instantiableConductorFPList, archLength, projectionDistance, projectionMergeDistance);
    clock_t tAfter = clock();
    tInstantiableConstruction = difftime(tAfter, tBefore)/CLOCKS_PER_SEC;

    return instantiableConductorFPList;
}

void GeoLoader::loadQui(const string &inputFileName) throw (FileNotFoundError)
{
    map<int, RectangleGLList> rectListMap;

    ifstream fin(inputFileName.c_str());

    if(!fin.is_open()){
        fin.close();
        throw FileNotFoundError(inputFileName);
    }

    string line;
    //* skip the first line
    getline(fin, line);
    int lineNumber = 1;
    while(getline(fin, line)){
        ++lineNumber;
        stringstream ss(line);
        char initial;
        ss >> initial;
        if (initial=='*'){
            //* skip comment
        }
        else if (initial != 'Q'){
            cerr << "ERROR: only support Q (line " << lineNumber << endl;
        }
        int condNum;
        float x1,y1,z1,x2,y2,z2,x3,y3,z3;
        ss >> condNum
           >> x1 >> y1 >> z1 >> x2 >> y2 >> z2 >> x3 >> y3 >> z3;
        RectangleGL rect;
        rect.x1 = min(x1, x3);
        rect.x2 = max(x1, x3);
        rect.y1 = min(y1, y3);
        rect.y2 = max(y1, y3);
        rect.z1 = min(z1, z3);
        rect.z2 = max(z1, z3);
        rect.xn = (rect.x1==rect.x2) ? 1 : 0;
        rect.yn = (rect.y1==rect.y2) ? 1 : 0;
        rect.zn = (rect.z1==rect.z2) ? 1 : 0;

        rectListMap[condNum].push_back(rect);
    }
    fin.close();

    //* copy to pwcConductorFPList
    pwcConductorFPList.clear();
    const int nLayer = 1;     // single layer for visualization purpose
    const int layerIndex = 0; // single layer
    const int dirIndex = 0;   // single dir
    int condIndex = 1;
    for ( map<int,RectangleGLList>::const_iterator each = rectListMap.begin();
          each != rectListMap.end(); ++each, ++condIndex){
        pwcConductorFPList.push_back(ConductorFP(nLayer, 0));
        RectangleGLList &rectList = pwcConductorFPList.back().layer[layerIndex][dirIndex];
        rectList.insert(rectList.begin(), each->second.begin(), each->second.end());
    }
}


//**
//* GeoLoader::runCapletQui
//* - write pwcConductorFPList to .qui in path
//* - run caplet
//* - write output to .stdsolver_output
ExtractionInfo &GeoLoader::runCapletQui(const std::string &pathFileBaseName)
        throw (FileNotFoundError)
{
    int coreNum = 1;

    const string program = "../caplet_solver/bin/capletMPI";
    const string suffix  = "stdsolver_output";
    const string resultSuffix = "stdsolver_result";
    const string outputFileName = pathFileBaseName+"."+suffix;
    writeFastcapFile(pathFileBaseName, pwcConductorFPList);

    stringstream ssCommand;
    ssCommand << "/usr/bin/mpirun -np " << coreNum << " " << program
              << " -o " << pathFileBaseName << ".cmat "
              << pathFileBaseName << ".qui | tee " << outputFileName;
    cout << ssCommand.str() << endl;
    int systemReturn = system(ssCommand.str().c_str());
    if (systemReturn != 0){
        cerr << "ERROR: Command returns nonzero value. (" << ssCommand.str() << ")" << endl;
    }

    ifstream fin(outputFileName.c_str());
    if (!fin.is_open()){
        fin.close();
        throw FileNotFoundError(outputFileName);
    }

    bool flagNCondFound    = false;
    bool flagNPanelFound   = false;
    bool flagTTotalFound   = false;
    bool flagTSetupFound   = false;
    bool flagMatrixFound   = false;

    extractionInfoList.push_back(ExtractionInfo());
    ExtractionInfo &result = extractionInfoList.back();

    string line;
    string temp;
    while(getline(fin,line)){
        if (line.empty()==true){
            continue;
        }
        if (flagNCondFound ==false && line.find("Number of conductors")  !=string::npos){
            flagNCondFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.nConductor;
        }
        if (flagNPanelFound ==false && line.find("Number of basis functions")  !=string::npos){
            flagNPanelFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.nBasisFunction;
        }
        if (flagTTotalFound ==false && line.find("Total extraction time")  !=string::npos){
            flagTTotalFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.tTotal;
        }
        if (flagTSetupFound ==false && line.find("Setup time")  !=string::npos){
            flagTSetupFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.tSetup;
        }
        if (flagMatrixFound==false && line.find("Cmat")!=string::npos){
            flagMatrixFound = true;

            //* read matrix
            result.capacitanceMatrix.resize(result.nConductor, vector<float>(result.nConductor));
            for (int i=0; i<result.nConductor; ++i){
                getline(fin, line);
                stringstream ss(line);
                for (int j=0; j<result.nConductor; ++j){
                    ss >> result.capacitanceMatrix[i][j];
                }
            }
        }

    }
    fin.close();

    result.tSolving = result.tTotal-result.tSetup;
    result.tTotal += tInstantiableConstruction;
    result.tBasis = tInstantiableConstruction;

    const string resultFileName = pathFileBaseName+"."+resultSuffix;

    ofstream fout(resultFileName.c_str());
    if (!fout.is_open()){
        fout.close();
        throw FileNotFoundError(resultFileName);
    }
    result.print(fout);
    fout.close();

    return result;
}


//**
//* GeoLoader::runFastcap
//* - write pwcConductorFPList to .qui in path
//* - run fastcap with flag
//* - write output to .fastcap_output
ExtractionInfo& GeoLoader::runFastcap(const string &pathFileBaseName, const string &option)
        throw (FileNotFoundError)
{
    const string program = "./fastcap";
    const string suffix  = "fastcap_output";
    const string resultSuffix = "fastcap_result";
    const string outputFileName = pathFileBaseName+"."+suffix;

    writeFastcapFile(pathFileBaseName, pwcConductorFPList);
    const string command = program + " " + option + " " + pathFileBaseName + ".qui | tee " + outputFileName;
    int systemReturn = system(command.c_str());
    if (systemReturn != 0){
        cerr << "ERROR: Command returns abnormally. (" << command << ")" << endl;
    }

    ifstream fin(outputFileName.c_str());
    if (!fin.is_open()){
        fin.close();
        throw FileNotFoundError(outputFileName);
    }

    bool flagNCondFound    = false;
    bool flagNPanelFound   = false;
    bool flagTTotalFound   = false;
    bool flagTSetupFound   = false;
    bool flagMatrixFound   = false;

    extractionInfoList.push_back(ExtractionInfo());
    ExtractionInfo &result = extractionInfoList.back();
    float unit = 1.0f;

    string line;
    string temp;
    while(getline(fin,line)){
        if (line.empty()==true){
            continue;
        }
        if (flagNCondFound ==false && line.find("Number of conductors:")  !=string::npos){
            flagNCondFound = true;
            stringstream ss(line);
            ss >> temp >>temp >> temp >> result.nConductor;
        }
        if (flagNPanelFound==false && line.find("Total number of panels:")!=string::npos){
            flagNPanelFound = true;
            stringstream ss(line);
            ss >> temp >> temp >>temp >> temp >> result.nBasisFunction;
        }
        if (flagTTotalFound==false && line.find("Total time:")!=string::npos){
            flagTTotalFound = true;
            stringstream ss(line);
            ss >> temp >> temp >> result.tTotal;
        }
        if (flagTSetupFound==false && line.find("Total setup time:")!=string::npos){
            flagTSetupFound = true;
            stringstream ss(line);
            ss >> temp >> temp >> temp >> result.tSetup;
        }
        if (flagMatrixFound==false && line.find("CAPACITANCE MATRIX,")!=string::npos){
            flagMatrixFound = true;

            //* read unit
            if (line.find("attofarads")!=string::npos){
                //* unit is attofarads
                unit = 1e-18;
            }
            else if (line.find("CAPACITANCE MATRIX, every unit is")!=string::npos){
                //* unit is represented in floating point
                stringstream ss(line);
                ss >> temp >> temp >> temp >> temp >> temp >> unit;
            }

            //* read matrix
            result.capacitanceMatrix.resize(result.nConductor, vector<float>(result.nConductor));
            getline(fin, line); // skip table header
            for (int i=0; i<result.nConductor; ++i){
                getline(fin, line);
                stringstream ss(line);
                ss >> temp >> temp;
                for (int j=0; j<result.nConductor; ++j){
                    ss >> result.capacitanceMatrix[i][j];
                    result.capacitanceMatrix[i][j] *= unit;
                }
            }
        }


    }
    fin.close();

    result.tSolving = result.tTotal-result.tSetup;
    result.tTotal += tPWCConstruction;
    result.tBasis = tPWCConstruction;

    const string resultFileName = pathFileBaseName+"."+resultSuffix;

    ofstream fout(resultFileName.c_str());
    if (!fout.is_open()){
        fout.close();
        throw FileNotFoundError(resultFileName);
    }
    result.print(fout);
    fout.close();

    return result;
}

ExtractionInfo &GeoLoader::runCaplet(const string &pathFileBaseName, const unsigned coreNum)
        throw (FileNotFoundError)
{
    const string program = "../caplet_solver/bin/capletMPI";
    const string suffix  = "caplet_output";
    const string resultSuffix = "caplet_result";
    const string outputFileName = pathFileBaseName+"."+suffix;

    writeCapletFile(pathFileBaseName, instantiableConductorFPList);
    stringstream ssCommand;
    ssCommand << "/usr/bin/mpirun -np " << coreNum << " " << program
              << " -o " << pathFileBaseName << ".cmat "
              << pathFileBaseName << ".caplet | tee " << outputFileName;
    cout << ssCommand.str() << endl;
    int systemReturn = system(ssCommand.str().c_str());
    if (systemReturn != 0){
        cerr << "ERROR: Command returns abnormally. (" << ssCommand.str() << ")" << endl;
    }

    ifstream fin(outputFileName.c_str());
    if (!fin.is_open()){
        fin.close();
        throw FileNotFoundError(outputFileName);
    }

    bool flagNCondFound    = false;
    bool flagNPanelFound   = false;
    bool flagTTotalFound   = false;
    bool flagTSetupFound   = false;
    bool flagMatrixFound   = false;

    extractionInfoList.push_back(ExtractionInfo());
    ExtractionInfo &result = extractionInfoList.back();

    string line;
    string temp;
    while(getline(fin,line)){
        if (line.empty()==true){
            continue;
        }
        if (flagNCondFound ==false && line.find("Number of conductors")  !=string::npos){
            flagNCondFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.nConductor;
        }
        if (flagNPanelFound ==false && line.find("Number of basis functions")  !=string::npos){
            flagNPanelFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.nBasisFunction;
        }
        if (flagTTotalFound ==false && line.find("Total extraction time")  !=string::npos){
            flagTTotalFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.tTotal;
        }
        if (flagTSetupFound ==false && line.find("Setup time")  !=string::npos){
            flagTSetupFound = true;
            size_t pos = line.find_last_of(':');
            line = line.substr(pos);
            stringstream ss(line);
            ss >> temp >> result.tSetup;
        }
        if (flagMatrixFound==false && line.find("Cmat")!=string::npos){
            flagMatrixFound = true;

            //* read matrix
            result.capacitanceMatrix.resize(result.nConductor, vector<float>(result.nConductor));
            for (int i=0; i<result.nConductor; ++i){
                getline(fin, line);
                stringstream ss(line);
                for (int j=0; j<result.nConductor; ++j){
                    ss >> result.capacitanceMatrix[i][j];
                }
            }
        }

    }
    fin.close();

    result.tSolving = result.tTotal-result.tSetup;
    result.tTotal += tInstantiableConstruction;
    result.tBasis = tInstantiableConstruction;

    const string resultFileName = pathFileBaseName+"."+resultSuffix;

    ofstream fout(resultFileName.c_str());
    if (!fout.is_open()){
        fout.close();
        throw FileNotFoundError(resultFileName);
    }
    result.print(fout);
    fout.close();

    return result;
}

const ExtractionInfo &GeoLoader::getLastResult() const
{
    return extractionInfoList.back();
}

const ExtractionInfo &GeoLoader::getReferenceResult() const
{
    return referenceResult;
}

void GeoLoader::loadReferenceResult(const string &filename)
{
    size_t size = geometryConductorFPList.size();

    ifstream fin(filename.c_str());

    //* read matrix
    referenceResult.capacitanceMatrix.resize(size, vector<float>(size));
    string line;
    for (size_t i=0; i<size; ++i){
        getline(fin, line);
        stringstream ss(line);
        for (size_t j=0; j<size; ++j){
            ss >> referenceResult.capacitanceMatrix[i][j];
            referenceResult.capacitanceMatrix[i][j];
        }
    }
    fin.close();
}

const ExtractionInfo &GeoLoader::storeLastAsReference(const string &pathFileNameCmat="")
{
    referenceResult = extractionInfoList.back();
    if (pathFileNameCmat.empty()==false){
        ofstream fout(pathFileNameCmat.c_str());
        referenceResult.printMatrix(fout);
        fout.close();
    }
    return referenceResult;
}

list<ExtractionInfo> GeoLoader::compareAllAgainstReference() const
{
    list<ExtractionInfo> result;
    for ( list<ExtractionInfo>::const_iterator each = extractionInfoList.begin();
          each != extractionInfoList.end(); ++each){
        result.push_back(ExtractionInfo());
        ExtractionInfo &current = result.back();
        current.nConductor = each->nConductor;
        current.nBasisFunction = each->nBasisFunction;
        current.tTotal = each->tTotal;
        current.tBasis = each->tBasis;
        current.tSetup = each->tSetup;
        current.tSolving = each->tSolving;
        current.error = each->compare(referenceResult);
    }
    return result;
}

void GeoLoader::clearResult()
{
    extractionInfoList.clear();
}

//__________________________________________________________
//*
//* Getters and setters
//*
size_t GeoLoader::getNumberOfConductor() const {
    return geometryConductorFPList.size();
}


//**
//* GeoLoader::generateConductorList
//* - Init: Copy metalConductorList to conductorList
ConductorList& GeoLoader::generateConductorList(ConductorList &conductorList, bool flagDecomposed)
{
    //* Init:
    //* For each ci on Layer k, check all vias on Layer k+, and add overlapping vias vj to ci
    //*     For each vj, add LayeredCond cl connected to vj on Layer k+1 to ci

    conductorList.clear();

    if (isLoaded==false){
        return conductorList;
    }

    //* Copy metalConductorList to conductorList
    conductorList.insert(conductorList.begin(), metalConductorList.begin(), metalConductorList.end());

    //* Construct 3D vias and put together connected conductors
    for ( int viaIndex = 0; viaIndex < nVia; ++viaIndex ){
        int lowerMetalIndex = viaConnect[viaIndex][0];
        int upperMetalIndex = viaConnect[viaIndex][1];

        for ( RectangleList::const_iterator eachViaIt = viaLayeredRectangleList[viaIndex].begin();
              eachViaIt != viaLayeredRectangleList[viaIndex].end(); ++eachViaIt )
        {
            list< Conductor >::iterator eachBottomCondIt;
            list< Conductor >::iterator eachTopCondIt;

            //* flags that indicate whether the metals connected to via bottom and via top found, respectively.
            bool flagBottom = false;
            bool flagTop = false;

            //* Search for metal conductors that connect to the bottom and the top of the via
            for ( list< Conductor >::iterator eachCondIt = conductorList.begin();
                    eachCondIt != conductorList.end(); ++eachCondIt ){
                if ( flagBottom == false && eachCondIt->isContaining(*eachViaIt, lowerMetalIndex) == true ){
                    flagBottom = true;
                    eachBottomCondIt = eachCondIt;
                }
                if ( flagTop == false && eachCondIt->isContaining(*eachViaIt, upperMetalIndex) == true ){
                    flagTop = true;
                    eachTopCondIt = eachCondIt;
                }
                //* Both metal connected to via bottom and to via top are found. Break loop.
                if ( flagBottom == true && flagTop == true ){
                    break;
                }
            }

            //* if the via is connected from both the bottom and the top
            if ( flagBottom == true && flagTop == true){
                //* If the conds are distinct, combine them
                if ( eachBottomCondIt != eachTopCondIt ){
                    *eachBottomCondIt += *eachTopCondIt;
                    conductorList.erase(eachTopCondIt);
                }
                //* generate a 3D via
                eachBottomCondIt->generateVia(*eachViaIt, viaIndex, viaDef, viaConnect, flagDecomposed); //[viaIndex][0], viaDef[viaIndex][1]);
            }else{
                if ( flagBottom == true ){
                    eachBottomCondIt->generateVia(*eachViaIt, viaIndex, viaDef, viaConnect, flagDecomposed); //[viaIndex][0], viaDef[viaIndex][1]);
                }
                if ( flagTop == true ){
                    eachTopCondIt->generateVia(*eachViaIt, viaIndex, viaDef, viaConnect, flagDecomposed); //[viaIndex][0], viaDef[viaIndex][1]);
                }
            }
        }

    }


    #ifdef DEBUG_PANEL_DISCRETIZATION
    //* check if any rect overlap each other (which is not correct)
    //* check if any rect has zero area (which is not correct)
    for (list<Conductor>::iterator eachCond = conductorList.begin();
            eachCond != conductorList.end(); ++eachCond ){
        if ( flagDecomposed == true && eachCond->checkSelfOverlapping(viaConnect) == true ){
            cerr << "ERROR: overlapping rectangles" << endl;
        }
        if ( eachCond->checkZeroAreaRectangle() == true ) {
            cerr << "ERROR: rectangles with zero area" << endl;
        }
    }
    #endif

    return conductorList;
}



//**
//* GeoLoader::readGeo
//* - aux function of loadGeo()
void GeoLoader::readGeo(
        const std::string   &geoFileName,
        LayeredPolygonList  &metalLayeredPolygonList,
        LayeredPolygonList  &viaLayeredPolygonList ) throw (FileNotFoundError)
{
    ifstream fin(geoFileName.c_str());
    if ( !fin.is_open() ){
        fin.close();
        throw FileNotFoundError(geoFileName);
    }

    clear();
    this->fileName = geoFileName;

    //* read layer info
    readLayerInfo(fin);

    //* read metal struc
    readStruc(fin, nMetal, metalLayeredPolygonList);
    //* read via struc
    readStruc(fin, nVia, viaLayeredPolygonList);

    //* close file
    fin.close();
}


//**
//* GeoLoader::readLayerInfo
//* - aux function of readGeo()
void GeoLoader::readLayerInfo(ifstream &fin)
{
    //* declare temp variables
    string line;
    int nLine = 0;
    int layerIndex = 0;
    int bottomElevation = 0;
    int topElevation = 0;
    int bottomConnect = 0;
    int topConnect = 0;

    //* get metal layer infomation
    getline(fin, line);
    sscanf(line.c_str(), "%d", &nLine);
    nMetal = nLine;
    //* allocate metalDef
    metalDef = new int*[nMetal];
    for ( int i=0; i<nMetal; i++ ){
        metalDef[i] = new int[2];
    }

    for ( int i=0; i<nMetal; i++ ){
        getline(fin, line);
        sscanf(line.c_str(), "%d, %d, %d", &layerIndex, &bottomElevation, &topElevation);
        metalDef[layerIndex][0] = bottomElevation;
        metalDef[layerIndex][1] = topElevation;
    }


    //* get via infomation
    getline(fin, line);
    sscanf(line.c_str(), "%d", &nLine);
    nVia = nLine;
    //* allocate viaDef, viaConnect
    viaDef = new int*[nVia];
    viaConnect = new int*[nVia];
    for ( int i=0; i<nVia; i++ ){
        viaDef[i] = new int[2];
        viaConnect[i] = new int[2];
    }

    for ( int i=0; i<nVia; i++ ){
        getline(fin, line);
        sscanf(line.c_str(), "%d, %d, %d, %d, %d", &layerIndex, &bottomElevation,
               &topElevation, &bottomConnect, &topConnect);
        viaDef[i][0] = bottomElevation;
        viaDef[i][1] = topElevation;
        viaConnect[i][0] = bottomConnect;
        viaConnect[i][1] = topConnect;
    }
}


//**
//* GeoLoader::readStruc
//* - aux function of readGeo()
//* - Store fin stream info to struc ( 2D array of vector of Point )
//* - input: fin, nLayer
//* - output: struc, nPoly
void GeoLoader::readStruc(ifstream &fin, int nLayer, vector<PolygonList> &struc){
    string  line;
    int     temp;
    int     nPoint;
    int     nPoly;
    int     xCoord;
    int     yCoord;

    PolygonList::reverse_iterator   eachPolygon;
    Polygon::reverse_iterator       eachPoint;

    struc.resize(nLayer);
    for ( int i=0; i<nLayer; ++i ){
        getline(fin, line);
        sscanf(line.c_str(), "%d", &temp);
        getline(fin, line);
        sscanf(line.c_str(), "%d", &nPoly);

        for ( int j=0; j<nPoly; j++ ){
            getline(fin, line);
            sscanf(line.c_str(), "%d", &nPoint);

            struc[i].push_back(Polygon());
            eachPolygon = struc[i].rbegin();

            for ( int k=0; k<nPoint-1; k++ ){
                getline(fin, line);
                sscanf(line.c_str(), "%d, %d", &xCoord, &yCoord);


                (*eachPolygon).push_back(Point());
                eachPoint = (*eachPolygon).rbegin();

                (*eachPoint).x = xCoord;
                (*eachPoint).y = yCoord;
            }
            getline(fin, line);
        }
    }
}


//void GeoLoader::layeredConductorList2ConductorList()
//{
//    conductorList.clear();
//    for ( int i=0; i<nMetal+nVia; ++i )
//    {
//        conductorList.insert(conductorList.end(), layeredMetalConductorList[i].begin(), layeredMetalConductorList[i].end());
//        layeredMetalConductorList[i].clear();
//    }
//}

//**
//* GeoLoader::printStruc
//* - Print struc to std::cout for each layer, each polygon, and each point
//* - input: nLayer, struc
void GeoLoader::printStruc(int nLayer, vector<PolygonList> &struc){
    PolygonList::iterator   eachPolygon;
    Polygon::iterator       eachPoint;

    for ( int i=0; i<nLayer; ++i ){
        cout << "Layer " << i << ":" << endl;
        int j = 0;
        for ( eachPolygon = struc[i].begin(); eachPolygon != struc[i].end(); ++eachPolygon, ++j ){
            cout << "  Polygon" << j << ":" << endl;

            for ( eachPoint = (*eachPolygon).begin(); eachPoint != (*eachPolygon).end(); ++eachPoint ){
                cout << "    ( " << (*eachPoint).x << ", " << (*eachPoint).y << " )" << endl;
            }
        }
    }
}

void GeoLoader::clear(){
    if (isLoaded){
        //* delete metalDef
        for (int i=0; i<nMetal; i++){
            delete[] metalDef[i];
        }
        delete[] metalDef;

        //* delete viaDef, viaConnect
        for (int i=0; i<nVia; i++){
            delete[] viaDef[i];
            delete[] viaConnect[i];
        }
        delete[] viaDef;
        delete[] viaConnect;

        //* clean up geo info
        viaLayeredRectangleList.clear();
        metalConductorList.clear();

        geometryConductorList.clear();
        geometryConductorFPList.clear();

        pwcConductorFPList.clear();

        //* clean up extraction info
        extractionInfoList.clear();
        referenceResult = ExtractionInfo();
    }
}

//* utility functions for poly2rect
bool xLess(const Point &p1, const Point &p2){
    return p1.x < p2.x;
}
void bring2front(Polygon &poly, Polygon::iterator &marker){
    poly.insert(poly.end(), poly.begin(), marker);
    poly.erase(poly.begin(), marker);
}
int sign(int a){
    return  ( a>0 ) ? 1 : ( (a<0) ? -1 : 0 );
}
bool lenLess(const Point &p1, const Point &p2){
    return p1.len < p2.len;
}
inline bool isConvex(int pos, int dir2){
    //* pos: the connection position of edge2 on edge1
    //* dir2: the inward direction of edge2
    return (pos*dir2)<0;
}
inline Dir getExtDir( Point &p1, Point &p2){
    Point temp = p2 - p1;

    if ( temp.x==0 ){
        //* in y-dir
        return (temp.y>0)? YP : YM;
    }else{
        //* in x-dir
        return (temp.x>0)? XP : XM;
    }
}
inline bool isCollinear(const Point &p1, const Point &p2, const Point &p3 ){
    if ( ((p1.x==p2.x)&&(p2.x==p3.x)) || ((p1.y==p2.y)&&(p2.y==p3.y)) ){
        return true;
    }else{
        return false;
    }
}

inline Polygon::iterator prevIterator(Polygon& poly, Polygon::iterator currentIterator){
    if ( currentIterator == poly.begin() ){
        return --poly.end();
    }else{
        return --currentIterator;
    }
}
inline Polygon::iterator nextIterator(Polygon& poly, Polygon::iterator currentIterator){
    ++currentIterator;
    if ( currentIterator == poly.end() ){
        return poly.begin();
    }else{
        return currentIterator;
    }
}
void updateEdgeDirLen( Point &prev, Point& current, Point &next )
{
    //* compute the dir and len of the current point
    //* based on the positions of all three points and the correct "prev.dir"

    Point vec_prev = current - prev;
    Point vec_current = next - current;

    //* determine the inward direction
    int rotation_dir = sign(vec_prev.x*vec_current.y - vec_prev.y*vec_current.x);
    if ( rotation_dir > 0 ){ //* right-hand rotation
        switch( prev.dir ){
        case XP:
            current.dir = YP; break;
        case XM:
            current.dir = YM; break;
        case YP:
            current.dir = XM; break;
        case YM:
            current.dir = XP; break;
        default:
            cerr << "ERROR: impossible direction detected!" << endl;
        }
    }else{ //* left-hand rotation
        switch( prev.dir ){
        case XP:
            current.dir = YM; break;
        case XM:
            current.dir = YP; break;
        case YP:
            current.dir = XP; break;
        case YM:
            current.dir = XM; break;
        default:
            cerr << "ERROR: impossible direction detected!" << endl;
        }
    }
    //* determine the edge length
    current.len = vec_current.vecLen2();
}


PolygonList extend(Polygon& poly, Polygon::iterator extensionPointIt,
               Polygon::iterator theOtherPointIt, Dir extDir)
{
    PolygonList cutSet;

    while (true){
        int dist = numeric_limits<int>::max();

        Polygon::iterator eachPointIt;
        Polygon::iterator closestPointIt;
        //* STEP 1. Find the closest point
        for (eachPointIt=poly.begin(); eachPointIt!=poly.end(); ++eachPointIt){
            Polygon::iterator nextPointIt = nextIterator(poly, eachPointIt);

            int tempDist = numeric_limits<int>::max();

            if ( abs(extDir) == X ){ //* extends to x-dir
                if ( ( (*eachPointIt).y-(*extensionPointIt).y )
                   * ( (*nextPointIt).y-(*extensionPointIt).y ) <= 0 )
                {
                    tempDist = ( (*eachPointIt).x - (*extensionPointIt).x )*sign(extDir);
                }
            }else{ //* extends to y-dir
                if ( ( (*eachPointIt).x-(*extensionPointIt).x )
                   * ( (*nextPointIt).x-(*extensionPointIt).x ) <= 0 )
                {
                    tempDist = ( (*eachPointIt).y - (*extensionPointIt).y )*sign(extDir);
                }
            }

            //* update the dist if closer
            if ( tempDist < dist && tempDist>0 ){
                dist = tempDist;
                closestPointIt = eachPointIt;
            }
        }
        if ( dist == numeric_limits<int>::max() ){
            //* if no closest edge is found
            return cutSet;
        }
        Polygon::iterator nextClosestPointIt = nextIterator(poly, closestPointIt);

        //* STEP 2. Extend the end point

        //* cut out the subset in either forward or backward indexing
        //* from insertedPoint to extensionPoint (keep the clock direction)
        Point insertedPoint;
        if ( abs(extDir)==X ){ //* x-dir
            insertedPoint.x = (*closestPointIt).x;
            insertedPoint.y = (*extensionPointIt).y;
        }else{ //* y-dir
            insertedPoint.x = (*extensionPointIt).x;
            insertedPoint.y = (*closestPointIt).y;
        }

        bool breakFlag = true;
        bool includeInsertedPointInCut = true;
        bool removeInsertedPointFromPoly = false;

        //* STEP 2 CASE 1:
        //* - Inserted point coincides with some point in poly
        //* - Same dir for collinear edges
        Polygon::iterator insertedPointIt;
        if ( (insertedPoint == (*nextClosestPointIt) &&
              theOtherPointIt == prevIterator(poly, extensionPointIt))
              ||
             (insertedPoint == (*closestPointIt) &&
              theOtherPointIt == nextIterator(poly, extensionPointIt)) )
        {
            //* the inserted point coincides the far end of the closest edge
            //* the next edge of the closest edge will be collinear with extensionEdge
            //* need to further extend extensionEdge and cut the U-shaped exterior polygon
            if (theOtherPointIt == prevIterator(poly, extensionPointIt)){
                insertedPointIt = nextClosestPointIt;
            }
            else{
                insertedPointIt = closestPointIt;
            }
            breakFlag = false;
            removeInsertedPointFromPoly = true;

        }
        //* STEP 2 CASE 2:
        //* - Inserted point coincides with some point in poly
        //* - Opposite dir for collinear edges
        else if ( (insertedPoint == (*closestPointIt) &&
                   theOtherPointIt == prevIterator(poly, extensionPointIt))
                   ||
                  (insertedPoint == (*nextClosestPointIt) &&
                   theOtherPointIt == nextIterator(poly, extensionPointIt)) )
        {
            //* the inserted point coincides the near end of the closest edge
            //* simply cut the exterior L-shaped polygon
            if (insertedPoint == (*closestPointIt)){
                insertedPointIt = closestPointIt;
            }
            else{
                insertedPointIt = nextClosestPointIt;
            }
            includeInsertedPointInCut = false;

        }
        //* STEP 2 OTHER CASES:
        //* - Inserted point is strictly on one edge in poly
        else{
            //* the inserted point is a new point in poly
            insertedPoint.dir = (*closestPointIt).dir;
            poly.insert(nextClosestPointIt, insertedPoint);
            insertedPointIt = nextClosestPointIt;
            insertedPointIt = prevIterator(poly, insertedPointIt);
        }
        //* STEP 3. Cut by the newly inserted extension edge
        Polygon cut;
        if ( theOtherPointIt == prevIterator(poly, extensionPointIt) ){
            //* forward cut from extensionPointIt to insertedPointIt
            for ( eachPointIt = extensionPointIt; eachPointIt != insertedPointIt; )
            {
                cut.push_back(*eachPointIt);
                eachPointIt = poly.erase(eachPointIt);
                eachPointIt = (eachPointIt==poly.end()) ? poly.begin() : eachPointIt;
            }
            if ( includeInsertedPointInCut ){
                cut.push_back(*eachPointIt);
                Polygon::reverse_iterator prevIt = ++cut.rbegin();
                Polygon::reverse_iterator prevPrevIt = prevIt; ++prevPrevIt; // cut.size()>=4. OK without boundary check
                updateEdgeDirLen( *prevIt, cut.back(), cut.front());
                updateEdgeDirLen(*prevPrevIt, *prevIt, cut.back());
            }
            if ( removeInsertedPointFromPoly ){
                insertedPointIt = nextIterator(poly, insertedPointIt);
                eachPointIt = poly.erase(eachPointIt);
                eachPointIt = (eachPointIt==poly.end()) ? poly.begin() : eachPointIt;
            }
            Polygon::iterator prevTheOtherPointIt = prevIterator(poly, theOtherPointIt);
            Polygon::iterator nextInsertedPointIt = nextIterator(poly, insertedPointIt);
            updateEdgeDirLen( *prevTheOtherPointIt, *theOtherPointIt, *insertedPointIt );
            updateEdgeDirLen( *theOtherPointIt, *insertedPointIt, *nextInsertedPointIt);
        }
        else{
            //* forward cut from insertedPointIt to extensionPointIt
            eachPointIt = insertedPointIt;
            if ( includeInsertedPointInCut ){
                cut.push_back(*eachPointIt);
            }
            if ( removeInsertedPointFromPoly ){
                insertedPointIt = prevIterator(poly, insertedPointIt);
                eachPointIt = poly.erase(eachPointIt);
                eachPointIt = (eachPointIt==poly.end()) ? poly.begin() : eachPointIt;
            }
            else{
                eachPointIt = nextIterator(poly, eachPointIt);
            }
            Polygon::iterator tempIt = extensionPointIt;
            tempIt = nextIterator(poly, tempIt);
            Polygon::iterator nextIt = nextIterator(poly, extensionPointIt);
            for ( ; eachPointIt != tempIt; )
            {
                cut.push_back(*eachPointIt);
                eachPointIt = poly.erase(eachPointIt);
                eachPointIt = (eachPointIt==poly.end()) ? poly.begin() : eachPointIt;
            }
            Polygon::reverse_iterator rprevIt = cut.rbegin();
            updateEdgeDirLen(*(++rprevIt), cut.back(), *cut.begin());
            updateEdgeDirLen(cut.back(), cut.front(), *(++cut.begin()));

            Polygon::iterator prevIt = prevIterator(poly, insertedPointIt);
            Polygon::iterator prevPrevIt = prevIterator(poly, prevIt);
            updateEdgeDirLen(*prevPrevIt, *prevIt, *insertedPointIt);
            updateEdgeDirLen(*prevIt, *insertedPointIt, *nextIt);
        }

        cutSet.push_back(cut);

        if ( breakFlag ){
            break;
        }
        else{
            //* reset extensionPointIt
            extensionPointIt = insertedPointIt;
        }
    }
    return cutSet;
}

PolygonList cut(Polygon &poly, Polygon::iterator cutPointIt){
    PolygonList cutSet;

    Polygon::iterator nextCutPointIt = nextIterator(poly, cutPointIt);
    Polygon::iterator prevCutPointIt = prevIterator(poly, cutPointIt);

    Point temp = *nextCutPointIt - *cutPointIt;
    int pos = -sign( temp.x + temp.y );

    #ifdef DEBUG_CUT
    printPolygon(poly);
    #endif

    //* Consider to extend edge on one end
    if ( !isConvex( pos, (*prevCutPointIt).dir ) ){
        //* concave. need extension
        Dir extDir = getExtDir(*nextCutPointIt, *cutPointIt);
        PolygonList tempPolygonList = extend( poly, cutPointIt, nextCutPointIt, extDir ) ;
        cutSet.insert( cutSet.end(), tempPolygonList.begin(), tempPolygonList.end() );
    }
    #ifdef DEBUG_CUT
    printPolygon(poly);
    printPolygonList(cutSet);
    #endif

    //* Consider to extend edge on the other end
    cutPointIt = prevIterator(poly, nextCutPointIt); // recompute cutPointIt
    if ( !isConvex(-pos, (*nextCutPointIt).dir ) ){
        //* concave. need extension
        Dir extDir = getExtDir(*cutPointIt, *nextCutPointIt);
        PolygonList tempPolygonList = extend( poly, nextCutPointIt, cutPointIt, extDir );
        cutSet.insert( cutSet.end(), tempPolygonList.begin(), tempPolygonList.end() );
    }
    #ifdef DEBUG_CUT
    printPolygon(poly);
    printPolygonList(cutSet);
    #endif

    return cutSet;
}


//**
//* poly2rect
void poly2rect(PolygonList &polyList, RectangleList &rectList){

    Polygon poly;

    //* STEP 1: Compute lens and dirs (inward normal) for each edge for each polygon
    //* - Each edge is represented by a point.
    //* - The other point of an edge is the next point iterator inferred by list data structure.
    for ( PolygonList::iterator eachPolygon = polyList.begin();
         eachPolygon != polyList.end(); ++eachPolygon )
    {
        //* Find the leftmost edge to get an idea of inner direction
        Polygon::iterator leftMostEdge = min_element( (*eachPolygon).begin(), (*eachPolygon).end(), xLess );
        if ( leftMostEdge == (*eachPolygon).begin() && (*leftMostEdge).x==(*eachPolygon).back().x ){
            //* fix the case when the edge connects the first and the last points
            leftMostEdge = --(*eachPolygon).end();
        }
        (*leftMostEdge).dir = XP;

        Polygon::iterator prevPointIt = leftMostEdge;
        Polygon::iterator thisPointIt = nextIterator( *eachPolygon, prevPointIt);
        Polygon::iterator nextPointIt = nextIterator( *eachPolygon, thisPointIt);
        for ( ; thisPointIt!=leftMostEdge;
              prevPointIt = thisPointIt,
              thisPointIt = nextPointIt,
              nextPointIt = nextIterator((*eachPolygon), nextPointIt) )
        {
            updateEdgeDirLen(*prevPointIt, *thisPointIt, *nextPointIt);
        }
        updateEdgeDirLen(*prevPointIt, *thisPointIt, *nextPointIt); // complete a loop of polygon edge

    }

    //* STEP 2: Decompose each polygon into a list of disjoint rectangles for a single metal layer
    //* - Pop the first polygon from polyList
    //* - Find the longest edge, sweep the edge to cover a rect area, and cut the rect from the polygon.
    //* - Push back the list of polygons that remain after the cut.
    //* - Keep popping the first polygon from polyList until polyList is empty.
    while ( !polyList.empty() )
    {
        poly = polyList.front();
        polyList.pop_front();

        //* if the poly contains only four points
        if ( poly.size() == 4 ){
            rectList.push_back( Rectangle(poly) );
            continue;
        }

        //* Find the longest edge and determine whether the vertex is convex or non-convex
        //* then cut both vertices of the longest edge if applicable
        Polygon::iterator longestEdgeIt = max_element( poly.begin(), poly.end(), lenLess );
        PolygonList tempPolygonList = cut(poly, longestEdgeIt);
        polyList.insert( polyList.end(), tempPolygonList.begin(), tempPolygonList.end() );

        //* update longestEdgeIt
        longestEdgeIt = max_element( poly.begin(), poly.end(), lenLess );
        Polygon::iterator nextLongestEdgeIt = nextIterator(poly, longestEdgeIt);

        #ifdef DEBUG_POLY2RECT
        printPolygon(poly);
        printPolygonList(tempPolygonList);
        #endif

        //* Find the cloest parallel edge
        Polygon::iterator eachPointIt;
        Polygon::iterator closestParallelEdgeIt;
        int dist = numeric_limits<int>::max();
        for ( eachPointIt = poly.begin(); eachPointIt != poly.end(); ++eachPointIt )
        {
            if ( eachPointIt == longestEdgeIt ){
                //* skip the current edge
                continue;
            }
            if ( (*eachPointIt).dir == -(*longestEdgeIt).dir )
            {
                //* opposite, inclosing parallel edge
                Polygon::iterator nextEachPointIt = nextIterator(poly, eachPointIt);
                int tempDist = numeric_limits<int>::max();

                int longestEdgeUpperLimit;
                int longestEdgeLowerLimit;
                int eachEdgeUpperLimit;
                int eachEdgeLowerLimit;
                if (abs((*eachPointIt).dir)==X){
                    //* x-dir
                    if ( (*longestEdgeIt).y > (*nextLongestEdgeIt).y ){
                        longestEdgeUpperLimit = (*longestEdgeIt).y;
                        longestEdgeLowerLimit = (*nextLongestEdgeIt).y;
                    }else{
                        longestEdgeUpperLimit = (*nextLongestEdgeIt).y;
                        longestEdgeLowerLimit = (*longestEdgeIt).y;
                    }
                    if ( (*eachPointIt).y > (*nextEachPointIt).y ){
                        eachEdgeUpperLimit = (*eachPointIt).y;
                        eachEdgeLowerLimit = (*nextEachPointIt).y;
                    }else{
                        eachEdgeUpperLimit = (*nextEachPointIt).y;
                        eachEdgeLowerLimit = (*eachPointIt).y;
                    }
                }else{
                    //* y-dir
                    if ( (*longestEdgeIt).x > (*nextLongestEdgeIt).x ){
                        longestEdgeUpperLimit = (*longestEdgeIt).x;
                        longestEdgeLowerLimit = (*nextLongestEdgeIt).x;
                    }else{
                        longestEdgeUpperLimit = (*nextLongestEdgeIt).x;
                        longestEdgeLowerLimit = (*longestEdgeIt).x;
                    }
                    if ( (*eachPointIt).x > (*nextEachPointIt).x ){
                        eachEdgeUpperLimit = (*eachPointIt).x;
                        eachEdgeLowerLimit = (*nextEachPointIt).x;
                    }else{
                        eachEdgeUpperLimit = (*nextEachPointIt).x;
                        eachEdgeLowerLimit = (*eachPointIt).x;
                    }
                }

                if ( !( eachEdgeLowerLimit >= longestEdgeUpperLimit || eachEdgeUpperLimit <= longestEdgeLowerLimit ) ){
                    // within longestEdge parallel range
                    tempDist =  (abs((*eachPointIt).dir)==X)? abs((*eachPointIt-*longestEdgeIt).x) : abs((*eachPointIt-*longestEdgeIt).y);
                }

                if ( tempDist < dist ){
                    dist = tempDist;
                    closestParallelEdgeIt = eachPointIt;
                }

            }
        }

        //* Extend the closest parallel edge and cut along the edge
        #ifdef DEBUG_POLY2RECT
        printPolygon(poly);
        #endif
        tempPolygonList = cut(poly, closestParallelEdgeIt);
        #ifdef DEBUG_POLY2RECT
        printPolygon(poly);
        printPolygonList(tempPolygonList);
        #endif

        polyList.insert( polyList.end(), tempPolygonList.begin(), tempPolygonList.end() );
        polyList.push_back( poly );
    }
}


void generateConnectedRects( RectangleList &rectList, ConnectedRectangleList &rectListList ){
    rectListList.clear();

    while( rectList.empty() == false ){
        //* create a new empty RectangleList
        rectListList.push_back(RectangleList());

        //* move the 1st item of rectList to the new rectangle list of rectListList
        rectListList.back().push_back(rectList.front());
        rectList.pop_front();

        //* loop over rectList to find any rect that is connected to the back pushed item
        for(RectangleList::iterator eachRectIt = rectListList.back().begin();
            eachRectIt != rectListList.back().end(); ++eachRectIt)
        {
            for(RectangleList::iterator originalRectIt = rectList.begin();
                originalRectIt != rectList.end(); )
            {
                if( !( eachRectIt->x2 < originalRectIt->x1 ||
                       eachRectIt->x1 > originalRectIt->x2 ||
                       eachRectIt->y2 < originalRectIt->y1 ||
                       eachRectIt->y1 > originalRectIt->y2    ) )
                {
                    //* if connected
                    //* then move *originalRectIt to rectListList.back()
                    rectListList.back().push_back(*originalRectIt);
                    rectList.erase(originalRectIt++);
                }
                else{
                    ++originalRectIt;
                }
            }
        }
    }
}


//**
//* computeAdjacency
//* - CURRENTLY DOES NOT SUPPORT SUBLAYERS
//* - Assume adjacent segments are disjoint (not work properly if not)
//* - Input a list of x-y plane rectangles
//* - Output the adjacency between rectangles and the complement
template<class T>
bool lessPairFirst(const pair<T,T> p1, const pair<T,T> p2){
    return p1.first < p2.first;
}
void computeAdjacency(const RectangleList               &rectList,
                      DirAdjacencyListOfRectangleList   &adjacency,
                      DirAdjacencyListOfRectangleList   &compAdjacency){


    enum ADJ {LEFT, RIGHT, BACK, FRONT, BOTTOM, TOP};

    adjacency.clear();
    compAdjacency.clear();

    for ( unsigned int i=0; i<rectList.size(); ++i ){
        adjacency.push_back(vector< list< pair<int,int> > >(4, list< pair<int,int> >()));
    }

    RectangleList::const_iterator               rectIit;
    RectangleList::const_iterator               rectListEndMinusOne = --rectList.end();
    DirAdjacencyListOfRectangleList::iterator   adjIit;
    for ( rectIit = rectList.begin(), adjIit = adjacency.begin();
          rectIit != rectListEndMinusOne; ++rectIit, ++adjIit )
    {

        RectangleList::const_iterator               rectJit = rectIit;
        DirAdjacencyListOfRectangleList::iterator   adjJit = adjIit;
        for ( ++rectJit, ++adjJit;
//        for ( rectJit = rectList.begin(), adjJit = adjacency.begin();
              rectJit != rectList.end(); ++rectJit, ++adjJit )
        {
//            if (rectJit == rectIit){
//                continue;
//            }

            //* I contains J in y-dir
            if ( (rectIit->y1 <= rectJit->y1 && rectJit->y1 <  rectIit->y2) ||
                 (rectIit->y1 <  rectJit->y2 && rectJit->y2 <= rectIit->y2) ||
                 (rectJit->y1 <= rectIit->y1 && rectIit->y1 <  rectJit->y2) ||
                 (rectJit->y1 <  rectIit->y2 && rectIit->y2 <= rectJit->y2)    )
            {
                //* connected from I's right
                if ( rectIit->x2 == rectJit->x1 ){
                    (*adjIit)[RIGHT].push_back(make_pair(rectJit->y1, rectJit->y2));
                    (*adjJit)[LEFT ].push_back(make_pair(rectIit->y1, rectIit->y2));
                }
                //* connected from I's left
                else if ( rectIit->x1 == rectJit->x2 ){
                    (*adjIit)[LEFT ].push_back(make_pair(rectJit->y1, rectJit->y2));
                    (*adjJit)[RIGHT].push_back(make_pair(rectIit->y1, rectIit->y2));
                }
            }
            //* I contains J in x-dir
            if ( (rectIit->x1 <= rectJit->x1 && rectJit->x1 <  rectIit->x2) ||
                 (rectIit->x1 <  rectJit->x2 && rectJit->x2 <= rectIit->x2) ||
                 (rectJit->x1 <= rectIit->x1 && rectIit->x1 <  rectJit->x2) ||
                 (rectJit->x1 <  rectIit->x2 && rectIit->x2 <= rectJit->x2)    )
            {
                //* connected from I's top
                if ( rectIit->y2 == rectJit->y1 ){
                    (*adjIit)[FRONT].push_back(make_pair(rectJit->x1, rectJit->x2));
                    (*adjJit)[BACK ].push_back(make_pair(rectIit->x1, rectIit->x2));
                }
                //* connected from I's bottom
                else if ( rectIit->y1 == rectJit->y2 ){
                    (*adjIit)[BACK ].push_back(make_pair(rectJit->x1, rectJit->x2));
                    (*adjJit)[FRONT].push_back(make_pair(rectIit->x1, rectIit->x2));
                }
            }
        }
    }


    //* sort and trim adjacent rect ranges
    const int nDir = 4; // LEFT, RIGHT, FRONT, BACK four directions
    for ( adjIit = adjacency.begin(), rectIit = rectList.begin(); adjIit!=adjacency.end(); ++adjIit, ++rectIit ){

        compAdjacency.push_back( DirAdjacencyList(nDir, AdjacencyList()));

        (*adjIit)[LEFT ].sort( lessPairFirst<int>);
        (*adjIit)[RIGHT].sort( lessPairFirst<int>);
        (*adjIit)[FRONT].sort( lessPairFirst<int>);
        (*adjIit)[BACK ].sort( lessPairFirst<int>);

        //* trim protrusive adjacent range and generate compAdjacency
        for ( int i=0; i<nDir; ++i ){
            int upper = 0;
            int lower = 0;

            if ( i==LEFT || i==RIGHT ){
                lower = rectIit->y1;
                upper = rectIit->y2;
            }else if( i==FRONT || i==BACK ){
                lower = rectIit->x1;
                upper = rectIit->x2;
            }

            //* if no adjacent rects, use full range
            if ( (*adjIit)[i].empty() ){
                compAdjacency.back()[i].push_back( make_pair(lower, upper) );
                continue;
            }

            if ( (*adjIit)[i].front().first > lower ){
                compAdjacency.back()[i].push_back( make_pair(lower, (*adjIit)[i].front().first) );
            }else{
                (*adjIit)[i].front().first = lower;
            }
            for ( list< pair<int,int> >::iterator eachRangeIt = (*adjIit)[i].begin();
                  eachRangeIt != --(*adjIit)[i].end(); ++eachRangeIt )
            {
                list< pair<int,int> >::iterator nextRangeIt = eachRangeIt;
                ++nextRangeIt;

                // if there is a gap between two adjacent rects
                if ( eachRangeIt->second != nextRangeIt->first ){
                    compAdjacency.back()[i].push_back( make_pair( eachRangeIt->second, nextRangeIt->first ) );
                }
            }
            if ( (*adjIit)[i].back().second < upper ){
                compAdjacency.back()[i].push_back( make_pair((*adjIit)[i].back().second, upper) );
            }else{
                (*adjIit)[i].back().second = upper;
            }
        }
    }

}



//**
//* generate3dRects
void generate3dRects(const RectangleList    &rect2dList,
                     const DirAdjacencyListOfRectangleList        &compAdjacency,
                     const int              elevationBottom,
                     const int              elevationTop,
                     const int              layerIndex,
                     Conductor              &cond)
{
    enum ADJ {LEFT, RIGHT, BACK, FRONT, BOTTOM, TOP};

    vector<RectangleList> &rect3dList = cond.layer[layerIndex];

    RectangleList::const_iterator each2dRectIt;
    list< vector< list< pair<int,int> > > >::const_iterator eachCompAdjIt;
    for ( each2dRectIt = rect2dList.begin(), eachCompAdjIt = compAdjacency.begin();
          each2dRectIt != rect2dList.end(); ++each2dRectIt, ++eachCompAdjIt )
    {
        //* generate top and bottom rects
        rect3dList[BOTTOM].push_back(*each2dRectIt);
        rect3dList[BOTTOM].back().normal = ZM;
        rect3dList[BOTTOM].back().z1 = elevationBottom;
        rect3dList[BOTTOM].back().z2 = elevationBottom;

        rect3dList[TOP].push_back(*each2dRectIt);
        rect3dList[TOP].back().normal = ZP;
        rect3dList[TOP].back().z1 = elevationTop;
        rect3dList[TOP].back().z2 = elevationTop;

        //* generate LEFT wall rects
        if( (*eachCompAdjIt)[LEFT].empty() == false ){
            for ( list< pair<int,int> >::const_iterator eachRangeIt = (*eachCompAdjIt)[LEFT].begin();
                  eachRangeIt != (*eachCompAdjIt)[LEFT].end(); ++eachRangeIt )
            {
                rect3dList[LEFT].push_back(Rectangle());
                rect3dList[LEFT].back().normal = XM;
                rect3dList[LEFT].back().x1 = each2dRectIt->x1;
                rect3dList[LEFT].back().x2 = each2dRectIt->x1;
                rect3dList[LEFT].back().y1 = eachRangeIt->first;
                rect3dList[LEFT].back().y2 = eachRangeIt->second;
                rect3dList[LEFT].back().z1 = elevationBottom;
                rect3dList[LEFT].back().z2 = elevationTop;
            }
        }

        //* generate RIGHT wall rects
        if( (*eachCompAdjIt)[RIGHT].empty() == false ){
            for ( list< pair<int,int> >::const_iterator eachRangeIt = (*eachCompAdjIt)[RIGHT].begin();
                  eachRangeIt != (*eachCompAdjIt)[RIGHT].end(); ++eachRangeIt )
            {
                rect3dList[RIGHT].push_back(Rectangle());
                rect3dList[RIGHT].back().normal = XP;
                rect3dList[RIGHT].back().x1 = each2dRectIt->x2;
                rect3dList[RIGHT].back().x2 = each2dRectIt->x2;
                rect3dList[RIGHT].back().y1 = eachRangeIt->first;
                rect3dList[RIGHT].back().y2 = eachRangeIt->second;
                rect3dList[RIGHT].back().z1 = elevationBottom;
                rect3dList[RIGHT].back().z2 = elevationTop;
            }
        }

        //* generate BACK wall rects
        if( (*eachCompAdjIt)[BACK].empty() == false ){
            for ( list< pair<int,int> >::const_iterator eachRangeIt = (*eachCompAdjIt)[BACK].begin();
                  eachRangeIt != (*eachCompAdjIt)[BACK].end(); ++eachRangeIt )
            {
                rect3dList[BACK].push_back(Rectangle());
                rect3dList[BACK].back().normal = YM;
                rect3dList[BACK].back().x1 = eachRangeIt->first;
                rect3dList[BACK].back().x2 = eachRangeIt->second;
                rect3dList[BACK].back().y1 = each2dRectIt->y1;
                rect3dList[BACK].back().y2 = each2dRectIt->y1;
                rect3dList[BACK].back().z1 = elevationBottom;
                rect3dList[BACK].back().z2 = elevationTop;
            }
        }

        //* generate FRONT wall rects
        if( (*eachCompAdjIt)[FRONT].empty() == false ){
            for ( list< pair<int,int> >::const_iterator eachRangeIt = (*eachCompAdjIt)[FRONT].begin();
                  eachRangeIt != (*eachCompAdjIt)[FRONT].end(); ++eachRangeIt )
            {
                rect3dList[FRONT].push_back(Rectangle());
                rect3dList[FRONT].back().normal = YP;
                rect3dList[FRONT].back().x1 = eachRangeIt->first;
                rect3dList[FRONT].back().x2 = eachRangeIt->second;
                rect3dList[FRONT].back().y1 = each2dRectIt->y2;
                rect3dList[FRONT].back().y2 = each2dRectIt->y2;
                rect3dList[FRONT].back().z1 = elevationBottom;
                rect3dList[FRONT].back().z2 = elevationTop;
            }
        }

    }
}

//**
//* discretizeDisjointSurface
//* - uniformly discretize rectangle surfaces with one-end length no larger
//*   than suggestedPanelSize
//* - used to generate PWC basis functions with disjoint rectangles as input
void discretizeXDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize);
void discretizeYDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize);
void discretizeZDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize);
void discretizeDisjointSurface(ConductorFPList &cond, const float suggestedPanelSize)
{
    for ( ConductorFPList::iterator eachCond = cond.begin();
          eachCond != cond.end(); ++eachCond){
        LayeredDirRectangleGLList &layer = eachCond->layer;
        for ( unsigned layerIndex = 0; layerIndex < layer.size(); ++layerIndex ){
            for ( unsigned dirIndex = 0; dirIndex < ConductorFP::nDir; ++dirIndex ){
                RectangleGLList &rectList = layer[layerIndex][dirIndex];

                //* push back a dummy rect as a marker (will be removed after the loop)
                rectList.push_back(RectangleGL());
                RectangleGLList::iterator rectEnd = --rectList.end();
                for ( RectangleGLList::iterator eachRectIt = rectList.begin();
                      eachRectIt != rectEnd; ){

                    if (eachRectIt->xn != 0){
                        //* x-dir
                        discretizeXDirRectangleGL(rectList, eachRectIt, suggestedPanelSize);
                    }
                    else if (eachRectIt->yn != 0){
                        //* y-dir
                        discretizeYDirRectangleGL(rectList, eachRectIt, suggestedPanelSize);
                    }
                    else if (eachRectIt->zn != 0){
                        //* z-dir
                        discretizeZDirRectangleGL(rectList, eachRectIt, suggestedPanelSize);
                    }
                    else{
                        ++eachRectIt;
                        cerr << "ERROR: impossible dir in discretizeDisjointSurface" << endl
                             << "       normal = (" << eachRectIt->xn << ","
                                                    << eachRectIt->yn << ","
                                                    << eachRectIt->zn << ")" << endl;
                    }
                }
                //* remove the empty marker
                rectList.erase(rectEnd);
            }
        }
    }
}


//**
//* discretizeZDirRectangleGL
//* - aux function of discretizeDisjointSurface
//* - discretize +z and -z RectangleGL
void discretizeZDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize)
{
    float xlen = rectIt->x2 - rectIt->x1;
    float ylen = rectIt->y2 - rectIt->y1;
    int nx = static_cast<int>(ceil(xlen/suggestedPanelSize));
    int ny = static_cast<int>(ceil(ylen/suggestedPanelSize));
    float dx = xlen/nx;
    float dy = ylen/ny;

    float y1 = rectIt->y1 - dy;
    float y2 = rectIt->y1;
    for ( int i = 0; i<ny; ++i ){
        y1 += dy;
        y2 += dy;
        float x1 = rectIt->x1 - dx;
        float x2 = rectIt->x1;

        for ( int j=0; j<nx; ++j ){
            x1 += dx;
            x2 += dx;

            rectList.push_back(*rectIt);
            RectangleGL &rect = rectList.back();

            rect.x1 = x1;
            rect.x2 = x2;
            rect.y1 = y1;
            rect.y2 = y2;
        }
    }
    rectIt = rectList.erase(rectIt);
}


//**
//* discretizeYDirRectangleGL
//* - aux function of discretizeDisjointSurface
//* - discretize +y and -y RectangleGL
void discretizeYDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize)
{
    float zlen = rectIt->z2 - rectIt->z1;
    float xlen = rectIt->x2 - rectIt->x1;

    int nz = static_cast<int>(ceil(zlen/suggestedPanelSize));
    int nx = static_cast<int>(ceil(xlen/suggestedPanelSize));

    float dz = zlen/nz;
    float dx = xlen/nx;

    float x1 = rectIt->x1 - dx;
    float x2 = rectIt->x1;
    for ( int i = 0; i<nx; ++i ){
        x1 += dx;
        x2 += dx;
        float z1 = rectIt->z1 - dz;
        float z2 = rectIt->z1;

        for ( int j=0; j<nz; ++j ){
            z1 += dz;
            z2 += dz;

            rectList.push_back(*rectIt);
            RectangleGL &rect = rectList.back();

            rect.z1 = z1;
            rect.z2 = z2;
            rect.x1 = x1;
            rect.x2 = x2;
        }
    }
    rectIt = rectList.erase(rectIt);
}


//**
//* discretizeXDirRectangleGL
//* - aux function of discretizeDisjointSurface
//* - discretize +x and -x RectangleGL
void discretizeXDirRectangleGL(
        RectangleGLList &rectList,
        RectangleGLList::iterator &rectIt,
        const float suggestedPanelSize)
{
    float ylen = rectIt->y2 - rectIt->y1;
    float zlen = rectIt->z2 - rectIt->z1;

    int ny = static_cast<int>(ceil(ylen/suggestedPanelSize));
    int nz = static_cast<int>(ceil(zlen/suggestedPanelSize));

    float dy = ylen/ny;
    float dz = zlen/nz;

    float z1 = rectIt->z1 - dz;
    float z2 = rectIt->z1;
    for ( int i = 0; i<nz; ++i ){
        z1 += dz;
        z2 += dz;
        float y1 = rectIt->y1 - dy;
        float y2 = rectIt->y1;

        for ( int j=0; j<ny; ++j ){
            y1 += dy;
            y2 += dy;

            rectList.push_back(*rectIt);
            RectangleGL &rect = rectList.back();

            rect.y1 = y1;
            rect.y2 = y2;
            rect.z1 = z1;
            rect.z2 = z2;
        }
    }
    rectIt = rectList.erase(rectIt);
}


//**
//* instantiateBasisFunction
//* - shapeShift is always set to zero
//* - do not extend the arch to the rect that is not right under the edge of
//*   the source rect (at this point).
//* - fixed normal distance.
void generateArch (RectangleGLList &rectList, const float archLength);
void instantiateBasisFunction (ConductorFPList &cond, const float archLength,
                               const float projectionDistance, const float projectionMergeDistance)
{
    float closestDistanceOfInterest = projectionDistance;

    unsigned nMetal = cond.front().nMetal;

    if (archLength<0){
        return;
    }

    for ( ConductorFPList::iterator eachCond1 = cond.begin();
          eachCond1 != cond.end(); ++eachCond1){
        for ( ConductorFPList::iterator eachCond2 = cond.begin();
              eachCond2 != cond.end(); ++eachCond2 ){
            //* if same conductor, skip
            if (eachCond1==eachCond2){
                continue;
            }
            //* determine overlapping rectangle
            //* 1. generate a list of overlapping rectangles
            //* 2. for each overlapping rectangle,
            //*    extend arches in four directions if dir1 is TOP or BOTTOM
            //*    extend arches in two  directions otherwise. (side walls)
            //* 3. for each extending arch,
            //*    find the overlapping rectangles that share the extending edge
            //* 4. If one arch is contained completed in another, erase it.
            //* The way to insert
            //* - insert them directly into the input list
            //* - return a list of iterators pointing to those inserted RectangleGL


            //* 1. generate a list of overlapping rectangles
            for( unsigned layer1=0; layer1<nMetal; ++layer1){
                for ( unsigned layer2=0; layer2<nMetal; ++layer2 ){
                    //* face-to-face
                    for (unsigned dir1=0; dir1<ConductorFP::nDir; ++dir1 ){
                        unsigned dir2 = ( dir1%2 == 0 ) ? (dir1+1) : (dir1-1);
                        RectangleGLList &rectList1 = eachCond1->layer[layer1][dir1];
                        RectangleGLList &rectList2 = eachCond2->layer[layer2][dir2];


                        RectangleGLList::IteratorList projectionRectItList;
                        for ( RectangleGLList::const_iterator eachRect2 = rectList2.begin();
                              eachRect2 != rectList2.end(); ++eachRect2){

                            //* If *eachRect2 is a projection, continue
                            if (eachRect2->shapeShift==1){
                                continue;
                            }

                            //* project rect from eachCond2 onto rects from eachCond1
                            RectangleGLList::IteratorList insertedRectItList =
                                    rectList1.insertProjectedOverlappingRectangleGL(*eachRect2, closestDistanceOfInterest);
                            projectionRectItList.insert(projectionRectItList.end(),
                                    insertedRectItList.begin(), insertedRectItList.end());
                        }

                    }
                    //* face-to-side
                    //- may not be necessary
                }
            }
        }
    }
    for ( ConductorFPList::iterator eachCond = cond.begin();
          eachCond != cond.end(); ++eachCond){

        for( unsigned layer=0; layer<nMetal; ++layer){
            for ( unsigned dir=0; dir<ConductorFP::nDir; ++dir){
                //* 1-a. remove rectangles with common support
                //* - make sure inserted projections are untouched.
                //* 1-b. remove projections with common support
                //* 2. genereate arches
                //* 3. for each extending arch,
                //*    find the overlapping rectangles that share the extending edge
                RectangleGLList &rectList = eachCond->layer[layer][dir];

                //* Ver1.0 Obsolete
                #ifdef MERGE_PROJECTION_VER1_0
                rectList.mergeProjection();
                rectList.absorbCommonSupport();
                #else
                //* Ver1.1
                rectList.mergeProjection1_1(projectionMergeDistance);
                rectList.absorbCommonSupport();
                rectList.removeBadProjection(caplet::DEFAULT_COINCIDENTAL_MARGIN);
                #endif

                if (archLength>0){
                    generateArch(rectList, archLength);
                }
            }
        }
    }
}
void intersectArch(RectangleGLList &rectList, RectangleGLList::iterator intersectEnd, RectangleGLList::iterator archPos);
void generateArch (RectangleGLList &rectList, const float archLength)
{

    //* find the head of projection flat rects
    RectangleGLList::iterator first = rectList.begin();
    for ( RectangleGLList::iterator each=rectList.begin(); each!=rectList.end(); ++each ){
        if (each->shapeShift!=0){
            first = each;
            break;
        }
    }

    if (first==rectList.begin()){
        //* there is no projection
        return;
    }
    for ( RectangleGLList::iterator each = first; each!=rectList.end(); ++each ){
        if (each->shapeType!=RectangleGL::FLAT_TYPE){
            continue;
        }
        RectangleGLList::iterator insertPos = each;
        ++insertPos;
        RectangleGLList::iterator archPos;

        if (each->xn!=0 || each->zn!=0){
            //* x-dir normal: two arches
            //* z-dir normal: bottom and top
            //* bottom
            archPos = rectList.insert(insertPos, *each);
            archPos->y1 = each->y1 - archLength;
            archPos->y2 = each->y1;
            archPos->shapeType = RectangleGL::ARCH_TYPE;
            archPos->shapeDir  = RectangleGL::Y_DECAY;
            archPos->shapeShift= 0;
            archPos->shapeNormalDistance = each->shapeNormalDistance * -1;
            intersectArch(rectList, first, archPos);
            //* top
            archPos = rectList.insert(insertPos, *each);
            archPos->y1 = each->y2;
            archPos->y2 = each->y2 + archLength;
            archPos->shapeType = RectangleGL::ARCH_TYPE;
            archPos->shapeDir  = RectangleGL::Y_DECAY;
            archPos->shapeShift= 0;
            archPos->shapeNormalDistance = each->shapeNormalDistance *  1;
            intersectArch(rectList, first, archPos);
        }
        if (each->yn!=0 || each->zn!=0){
            //* y-dir normal: two arches
            //* z-dir normal: left and right
            //* left
            archPos = rectList.insert(insertPos, *each);
            archPos->x1 = each->x1 - archLength;
            archPos->x2 = each->x1;
            archPos->shapeType = RectangleGL::ARCH_TYPE;
            archPos->shapeDir  = RectangleGL::X_DECAY;
            archPos->shapeShift= 0;
            archPos->shapeNormalDistance = each->shapeNormalDistance * -1;
            intersectArch(rectList, first, archPos);
            //* right
            archPos = rectList.insert(insertPos, *each);
            archPos->x1 = each->x2;
            archPos->x2 = each->x2 + archLength;
            archPos->shapeType = RectangleGL::ARCH_TYPE;
            archPos->shapeDir  = RectangleGL::X_DECAY;
            archPos->shapeShift= 0;
            archPos->shapeNormalDistance = each->shapeNormalDistance *  1;
            intersectArch(rectList, first, archPos);
        }
    }
}
void intersectArch(RectangleGLList &rectList, RectangleGLList::iterator intersectEnd, RectangleGLList::iterator archPos)
{
    for ( RectangleGLList::iterator eachRect = rectList.begin();
          eachRect != intersectEnd; ++eachRect){

        RectangleGL rect = archPos->intersectArchOnFlat(*eachRect);
        if (rect.isEmpty()==false){
            rectList.insert(archPos, rect);
        }
    }
    archPos = rectList.erase(archPos);
}


//****
//*
//* File writer
//*
//*
void writeFastcapFile(
        const std::string &outputFileName,
        const ConductorFPList &cond)
        throw (FileNotFoundError)
{
    string fullFileName = outputFileName + ".qui";
    ofstream fout(fullFileName.c_str());
    if (fout.is_open()==false){
        fout.close();
        throw FileNotFoundError(outputFileName);
    }

    //* first line
    fout << "0 " << outputFileName << endl;

    unsigned condIndex = 1;
    for ( ConductorFPList::const_iterator eachCond = cond.begin();
          eachCond != cond.end(); ++eachCond, ++condIndex ){
        const LayeredDirRectangleGLList &layer = eachCond->layer;
        for ( unsigned layerIndex = 0; layerIndex < layer.size(); ++layerIndex ){
            for ( unsigned dirIndex = 0; dirIndex < ConductorFP::nDir; ++dirIndex){
                const RectangleGLList *eachList = &layer[layerIndex][dirIndex];
                for ( RectangleGLList::const_iterator each = eachList->begin();
                      each != eachList->end(); ++each ){
                    fout << "Q " << condIndex << "    ";

                    if (each->zn > 0){
                        //* +z-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << endl;
                    }
                    else if (each->zn < 0){
                        //* -z-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << endl;
                    }
                    else if (each->xn > 0){
                        //* +x-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << endl;
                    }
                    else if (each->xn < 0){
                        //* -x-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y2
                             << setw(14) << setprecision(6) << each->z1 << endl;
                    }
                    else if (each->yn > 0){
                        //* +y-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << endl;
                    }
                    else if (each->yn < 0){
                        //* -y-dir
                        fout << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z1 << "    "

                             << setw(14) << setprecision(6) << each->x2
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << "    "

                             << setw(14) << setprecision(6) << each->x1
                             << setw(14) << setprecision(6) << each->y1
                             << setw(14) << setprecision(6) << each->z2 << endl;
                    }
                }
            }
        }
    }


    fout.close();
}

void writeCapletFile(
        const std::string &outputFileName,
        const ConductorFPList &cond)
        throw (FileNotFoundError)
{
    string fullFileName = outputFileName + ".caplet";
    ofstream fout(fullFileName.c_str());
    if (fout.is_open()==false){
        fout.close();
        throw FileNotFoundError(outputFileName);
    }

    fout << cond.size() << endl;
    unsigned totalCount = 0;
    for ( ConductorFPList::const_iterator eachCond = cond.begin();
          eachCond != cond.end(); ++eachCond)
    {
        unsigned thisCondSize = eachCond->size();
        fout << thisCondSize << " ";
        totalCount += thisCondSize;
    }
    fout << endl;
    fout << totalCount << endl;
    for ( ConductorFPList::const_iterator eachCond = cond.begin();
          eachCond != cond.end(); ++eachCond){
        for (unsigned layerIndex=0; layerIndex<eachCond->layer.size(); ++layerIndex){
            for (unsigned dirIndex=0; dirIndex<ConductorFP::nDir; ++dirIndex){
                const RectangleGLList &rectList = eachCond->layer[layerIndex][dirIndex];
                for ( RectangleGLList::const_iterator eachRect = rectList.begin();
                      eachRect != rectList.end(); ++eachRect){
                    eachRect->printCapletLine(fout);
                }
            }
        }
    }

    fout.close();
}



void printPolygon(Polygon &poly){
    Polygon::iterator eachPoint;
    cout << "  Polygon" << endl;
    cout << "  size = " << poly.size() << endl;
    cout << " [ ..." << endl;
    for ( eachPoint = poly.begin(); eachPoint != poly.end(); ++eachPoint ){
        cout << (*eachPoint).x << ", " << (*eachPoint).y << "..." << endl;
    }
}

void printPolygonList(PolygonList &polyList){
    PolygonList::iterator eachPolygon;
    cout << "PolygonList" << endl;
    cout << "size = " << polyList.size() << endl;
    for ( eachPolygon = polyList.begin(); eachPolygon!=polyList.end(); ++eachPolygon ){
        printPolygon(*eachPolygon);
    }
}


void printRectList(RectangleList &rectList){
    RectangleList::iterator eachRectIt;
    cout << "RectList" << endl;
    cout << "size = " << rectList.size() << endl;
    for ( eachRectIt = rectList.begin(); eachRectIt != rectList.end(); ++eachRectIt ){
        cout << "  ( " << (*eachRectIt).x1 << ", " << (*eachRectIt).x2 << ", "
             << (*eachRectIt).y1 << ", " << (*eachRectIt).y2 << ", "
             << (*eachRectIt).z1 << ", " << (*eachRectIt).z2 << " ), area: "
             << (*eachRectIt).area() << endl;
    }
}


void printRectListMatlab(RectangleList &rectList, Dir dir, ofstream& fout){
    if ( rectList.size() == 0 ){
        return;
    }

    RectangleList::iterator eachItemIt;
    // x-coord
    if ( dir == X ){
        fout << "x = [..." << endl;
        for ( eachItemIt = rectList.begin(); eachItemIt != rectList.end(); ++eachItemIt ){
            Rectangle eachRect = (*eachItemIt);
            if ( eachRect.normal == XP || eachRect.normal == XM ){
                fout << eachRect.x1 << ", " << eachRect.x2 << ", " << eachRect.x2 << ", " << eachRect.x1 << "; ..." <<endl;
            }else if( eachRect.normal == YP || eachRect.normal == YM ){
                fout << eachRect.x1 << ", " << eachRect.x2 << ", " << eachRect.x2 << ", " << eachRect.x1 << "; ..." <<endl;
            }else{
                fout << eachRect.x1 << ", " << eachRect.x2 << ", " << eachRect.x2 << ", " << eachRect.x1 << "; ..." <<endl;
            }
        }
        fout << "];" <<endl;
        fout << "x = x';" << endl;
    }

    // y-coord
    if ( dir == Y ){
        fout << "y = [..." << endl;
        for ( eachItemIt = rectList.begin(); eachItemIt != rectList.end(); ++eachItemIt ){
            Rectangle eachRect = (*eachItemIt);
            if ( eachRect.normal == XP || eachRect.normal == XM ){
                fout << eachRect.y1 << ", " << eachRect.y2 << ", " << eachRect.y2 << ", " << eachRect.y1 << "; ..." <<endl;
            }else if( eachRect.normal == YP || eachRect.normal == YM ){
                fout << eachRect.y1 << ", " << eachRect.y2 << ", " << eachRect.y2 << ", " << eachRect.y1 << "; ..." <<endl;
            }else{
                fout << eachRect.y1 << ", " << eachRect.y1 << ", " << eachRect.y2 << ", " << eachRect.y2 << "; ..." <<endl;
            }
        }
        fout << "];" <<endl;
        fout << "y = y';" << endl;
    }

    // z-coord
    if ( dir == Z ){
        fout << "z = [..." << endl;
        for ( eachItemIt = rectList.begin(); eachItemIt != rectList.end(); ++eachItemIt ){
            Rectangle eachRect = (*eachItemIt);
            if ( eachRect.normal == XP || eachRect.normal == XM ){
                fout << eachRect.z1 << ", " << eachRect.z1 << ", " << eachRect.z2 << ", " << eachRect.z2 << "; ..." <<endl;
            }else if( eachRect.normal == YP || eachRect.normal == YM ){
                fout << eachRect.z1 << ", " << eachRect.z1 << ", " << eachRect.z2 << ", " << eachRect.z2 << "; ..." <<endl;
            }else{
                fout << eachRect.z1 << ", " << eachRect.z1 << ", " << eachRect.z2 << ", " << eachRect.z2 << "; ..." <<endl;
            }
        }
        fout << "];" <<endl;
        fout << "z = z';" << endl;
    }
}

void printConductorListMatlab(ConductorList &conductorList){
    string filename("structure_output.m");

    ofstream fout(filename.c_str());

    fout << "figure; hold all;" << endl;

    char color[] = {  'r','g','b',  'c','m','y','w','k' };
    int k=0;
    int Nk = 8;

    for ( list<Conductor>::iterator eachCondIt = conductorList.begin();
          eachCondIt != conductorList.end(); ++eachCondIt )
    {
        // layer
        for ( int i=0; i<eachCondIt->nLayer; ++i ){
            // dir
            for ( int j=0; j<6; j++ ){
                if ( eachCondIt->layer[i][j].empty() == false ){
                    printRectListMatlab(eachCondIt->layer[i][j], X, fout);
                    printRectListMatlab(eachCondIt->layer[i][j], Y, fout);
                    printRectListMatlab(eachCondIt->layer[i][j], Z, fout);
                    fout<< "fill3(x,y,z, "<<"'"<<color[k] <<"')" << endl;
                }
            }
        }
        k = (k+1)%Nk;
    }


    fout.close();
}
void printRectMap(RectangleMap &rectMap){
    RectangleMap::iterator eachItemIt;
    cout << "RectMap" << endl;
    cout << "size = " << rectMap.size() << endl;

    for ( eachItemIt = rectMap.begin(); eachItemIt != rectMap.end(); ++eachItemIt ){
        Rectangle eachRect = (*eachItemIt).second;
        double area = (*eachItemIt).first;
        cout << "  ( " << eachRect.x1 << ", " << eachRect.x2 << ", "
             << eachRect.y1 << ", " << eachRect.y2 << ", "
             << eachRect.z1 << ", " << eachRect.z2 << " ), area: "
             << area << endl;
    }
}


//****
//*
//* ExtractionInfo
//*
//*

//**
//* ExtractionInfo::printBasic
void ExtractionInfo::printBasic(ostream &out) const
{
    out << "Number of conductors        : " << this->nConductor << endl;
    out << "Number of basis functions   : " << this->nBasisFunction << endl;
    out << "Total extraction time (s)   : " << this->tTotal << endl;
    out << "  Setup time (s)            : " << this->tSetup << endl;
    out << "  Solving time (s)          : " << this->tSolving << endl;
    out << "Basis Construction time (s) : " << this->tBasis << endl;
}

void ExtractionInfo::printBasicHeader(ostream &out) const
{
    out << "#conductors = " << nConductor << endl;
    out << "#basis_functions, error(%), total_time(s), basis_construction_time(s), setup_time(s), solving_time(s)";
}

void ExtractionInfo::printBasicLine(ostream &out) const
{
    out << setw(5)  << nBasisFunction << ", "
        << setw(8) << setprecision(4) << error*100 << ", "
        << setw(8) << setprecision(4) << tTotal << ", "
        << setw(8) << setprecision(4) << tBasis << ", "
        << setw(8) << setprecision(4) << tSetup << ", "
        << setw(8) << setprecision(4) << tSolving ;
}

//**
//* ExtractionInfo::printMatrix
void ExtractionInfo::printMatrix(ostream &out) const
{
    if (capacitanceMatrix.size()==0){
        return;
    }

    for ( int i=0; i<nConductor; ++i ){
        for ( int j=0; j<nConductor; ++j ){
            out << setw(14) << setprecision(6) << capacitanceMatrix[i][j];
        }
        out << endl;
    }
}

//**
//* ExtractionInfo::print
void ExtractionInfo::print(ostream &out) const
{
    printBasic(out);
    out << "Matrix: " << endl;
    printMatrix(out);
}

string ExtractionInfo::toString() const
{
    stringstream ss;
    print(ss);
    return ss.str();
}

float ExtractionInfo::compare(const ExtractionInfo &ref) const
throw (length_error)
{
    if (ref.capacitanceMatrix.size()==0){
        return numeric_limits<float>::quiet_NaN();
    }

    if (this->capacitanceMatrix.size() != ref.capacitanceMatrix.size()){
        std::cerr << capacitanceMatrix.size() << std::endl;
        std::cerr << ref.capacitanceMatrix.size() << std::endl;
        throw length_error("Matrix sizes are not compatiable");
    }

    float maxErr = numeric_limits<float>::min();
    unsigned size = this->capacitanceMatrix.size();

    for (unsigned i=0; i<size; ++i){
        for (unsigned j=0; j<size; ++j){
            float err = abs((this->capacitanceMatrix[i][j]-ref.capacitanceMatrix[i][j])/ref.capacitanceMatrix[i][i]);
            maxErr = max(maxErr, err);
        }
    }
    return maxErr;
}

vector<float> ExtractionInfo::compareDiagonal(const ExtractionInfo &ref) const
throw (length_error)
{
    if (ref.capacitanceMatrix.size()==0){
        return vector<float>(1, numeric_limits<float>::quiet_NaN());
    }

    if (this->capacitanceMatrix.size() != ref.capacitanceMatrix.size()){
        throw length_error("Matrix sizes are not compatiable");
    }

    vector<float> errorVector(nConductor, 0);

    for (int i=0; i<nConductor; ++i){
        errorVector[i] = abs((this->capacitanceMatrix[i][i]-ref.capacitanceMatrix[i][i])/ref.capacitanceMatrix[i][i]);
    }
    return errorVector;
}
























