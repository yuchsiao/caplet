/*
Created : Jul 28, 2010
Modified: Feb 13, 2013
Author  : Yu-Chung Hsiao
Email   : yuchsiao@mit.edu
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



#include "caplet_const.h"
#include "caplet_parameter.h"
#include "caplet.h"
#include "caplet_blas.h"
#include "caplet_widgets.h"
#include "caplet_elem.h"
#include "caplet_int.h"
#include "caplet_gauss.h"

#include "mpi.h"

#include <fstream>
#include <sstream>
#include <iostream>
#include <string>
#include <list>
#include <cmath>
#include <cstdlib>
#include <iomanip>


namespace caplet{

using namespace std;

//**** Private constants
float Caplet::epsilon0 	= 8.8541878176e-12f;
float Caplet::pi		= 3.1415926535f;


//**** Public functions
Caplet::Caplet()
    : isLoaded(false), isSolved(false){
}


Caplet::~Caplet(){
    this->clear();
}


void Caplet::clear(){
    if ( this->isLoaded == true){
	    this->isLoaded = false;
        this->isInstantiable = false;

        delete[] this->nWirePanels;
        delete[] this->nWireCoefs;

        delete[] this->panels;

        delete[] this->dirs;
        delete[] this->areas;
        delete[] this->indexIncrements;
        delete[] this->basisTypes;
        delete[] this->basisDirs;
        delete[] this->basisZs;
        delete[] this->basisShifts;

		if ( this->isSolved == true ){
			this->isSolved = false;

            //* single precision version
            switch(this->mode){
            case FAST_GALERKIN:
                delete[] this->P;
                delete[] this->rhs;
                delete[] this->coefs;

                delete[] this->Cmat;
                break;
            case DOUBLE_COLLOCATION:
                delete[] this->dP;
                delete[] this->drhs;
                delete[] this->dcoefs;

                delete[] this->dCmat;
                break;
            default:
                ;
            }
        }
	}
    this->mode = FAST_GALERKIN;
    this->isInstantiable = false;
}


void Caplet::loadFastcapFile(const std::string filename){

    //* Clean up
	this->clear();

    ifstream ifile(filename.c_str());
    //* Check if open successfully
	if (!ifile){
        cerr << "ERROR: Cannot open the file: " << filename << endl;
		return;
	}
    
    //* Initialize
	this->nPanels = 0;
	this->nWires = 0;

    string		stringTemp;
    string		lineTemp;
    string      currentConductorName = "";
    char        charTemp;
    list<int>	nWirePanelsList;
    int			nWirePanelTemp = 0;

	float xCoord[4];
	float yCoord[4];
	float zCoord[4];

	float xmin, xmax, ymin, ymax, zmin, zmax;

    list<float> x1;
    list<float> x2;
    list<float> y1;
    list<float> y2;
    list<float> z1;
    list<float> z2;

    //* skip the first line
	getline(ifile, lineTemp);

	while( getline(ifile, lineTemp) ){
		nWirePanelTemp++;
		this->nPanels++;

		std::stringstream	stringTokenizer(lineTemp);

		stringTokenizer >> charTemp;
		if (charTemp != 'Q'){
            cerr << "ERROR: only 'Q' is allowed in this solver!" << endl;
		}

		stringTokenizer >> stringTemp;

		stringTokenizer >> xCoord[0]
						>> yCoord[0]
						>> zCoord[0]

						>> xCoord[1]
						>> yCoord[1]
						>> zCoord[1]

						>> xCoord[2]
						>> yCoord[2]
						>> zCoord[2]

						>> xCoord[3]
						>> yCoord[3]
						>> zCoord[3];

		xmin = xCoord[0];
		xmax = xCoord[0];
		ymin = yCoord[0];
		ymax = yCoord[0];
		zmin = zCoord[0];
		zmax = zCoord[0];

		for (int i=1; i<4; i++){
			if ( xmin > xCoord[i] ) xmin = xCoord[i];
			if ( xmax < xCoord[i] ) xmax = xCoord[i];
			if ( ymin > yCoord[i] ) ymin = yCoord[i];
			if ( ymax < yCoord[i] ) ymax = yCoord[i];
			if ( zmin > zCoord[i] ) zmin = zCoord[i];
			if ( zmax < zCoord[i] ) zmax = zCoord[i];
		}

		x1.push_back(xmin);
		x2.push_back(xmax);
		y1.push_back(ymin);
		y2.push_back(ymax);
		z1.push_back(zmin);
		z2.push_back(zmax);

		if ( currentConductorName.compare(stringTemp) != 0 ){ // is different
            //* change to another conductor
			currentConductorName = stringTemp;
			this->nWires++;
			nWirePanelsList.push_back(nWirePanelTemp);
			nWirePanelTemp = 0;
		}
	}
    //* Register the last wire information
	nWirePanelsList.push_back(++nWirePanelTemp);

    //* Collect all information from the .caplet file until here
    //  fill in and allocate each array and field below

	this->nCoefs = this->nPanels;

    //* Convert nWirePanelsList to the nWirePanels array
	this->nWirePanels = new int[this->nWires];
	this->nWireCoefs  = new int[this->nWires];
	int i=0;
	for ( std::list<int>::iterator it = ++nWirePanelsList.begin();
			it != nWirePanelsList.end(); ++it){
		this->nWirePanels[i] = *it;
		this->nWireCoefs [i] = *it;
		i++;
	}

	this->panels 	= new float	[this->nPanels][3][4];
	this->dirs 		= new int	[this->nPanels];
	this->areas		= new float	[this->nPanels];
	for ( int i=0; i<this->nPanels; i++){
		this->panels[i][X][MIN] = x1.front(); x1.pop_front();
		this->panels[i][X][MAX] = x2.front(); x2.pop_front();
		this->panels[i][X][LENGTH]
			   = this->panels[i][X][MAX] - this->panels[i][X][MIN];
		this->panels[i][X][CENTER]
			   = (this->panels[i][X][MAX] + this->panels[i][X][MIN])/2;

		this->panels[i][Y][MIN] = y1.front(); y1.pop_front();
		this->panels[i][Y][MAX] = y2.front(); y2.pop_front();
		this->panels[i][Y][LENGTH]
			   = this->panels[i][Y][MAX] - this->panels[i][Y][MIN];
		this->panels[i][Y][CENTER]
			   = (this->panels[i][Y][MAX] + this->panels[i][Y][MIN])/2;

		this->panels[i][Z][MIN] = z1.front(); z1.pop_front();
		this->panels[i][Z][MAX] = z2.front(); z2.pop_front();
		this->panels[i][Z][LENGTH]
			   = this->panels[i][Z][MAX] - this->panels[i][Z][MIN];
		this->panels[i][Z][CENTER]
			   = (this->panels[i][Z][MAX] + this->panels[i][Z][MIN])/2;

		if ( std::abs(this->panels[i][X][LENGTH]) < zero ){ // in x-dir
			this->dirs[i] = X;
			this->areas[i] = this->panels[i][Y][LENGTH] * this->panels[i][Z][LENGTH];
		}else if ( std::abs(this->panels[i][Y][LENGTH]) < zero ){ // in y-dir
			this->dirs[i] = Y;
			this->areas[i] = this->panels[i][Z][LENGTH] * this->panels[i][X][LENGTH];
		}else{ // in z-dir
			this->dirs[i] = Z;
			this->areas[i] = this->panels[i][X][LENGTH] * this->panels[i][Y][LENGTH];
		}

	}

    //* For a uniform treatment of generateRHS
	this->indexIncrements 	= new int[this->nPanels];
	this->basisTypes		= new char[this->nPanels];
	for (int i=0; i<this->nPanels; i++){
		this->indexIncrements[i] = 1;
		this->basisTypes[i] = 'F';
	}


    //* Allocate dummy array for a uniform clear operation
	this->basisDirs			= new int[this->nPanels];
	for ( int i=0; i< this->nPanels; i++ ){
		this->basisDirs[i] = FLAT;
	}
	this->basisZs			= new float[1];
	this->basisShifts		= new float[1];

    //* Set flags
	this->isSolved = false;
	this->isLoaded = true;

	ifile.close();

    #ifdef DEBUG_ASPECT_RATIO_VALIDITY
    if( this->isPanelAspectRatioValid() == false ){
        cerr << "ERROR: Panel aspect ratio of FASTCAP file is not valid." << endl;
        exit(1);
    }
    #endif

}


void Caplet::loadCapletFile(const std::string filename){
    //* Caplet format
	//* Line 1: nWire
	//* Line 2: nShape for each wire
	//* Line 3: nTotalShape
    //* Line n: shape description lines, 12 numbers per line
	//* shapeType, indexIncrement, XL, XU, YL, YU, ZL, ZU, dir, shapeDir, shapeNormalDistance, shapeShift
	//* shapeType: either F for Flat, A for Arch, and S for Side arch
	//* indexIncrement: 
	//*   a basis function can consists of multiple shapes. The value here can be 0 or 1.
	//*   1: the beginning of a new basis function
	//*   0: not increment, meaning the current shape is combined with the shape before with value 1
	//* XL, XU, YL, YU, ZL, ZU: the range limit of each direction
	//* dir: normal direction of the rectangle that supports the shape
	//*   0 for x-dir, 1 for y-dir, 2 for z-dir
	//* shapeDir: shape varying direction
	//*   0 for x-dir, 1 for y-dir, 2 for z-dir
	//* shapeNormalDistrance: signed distance between the support rectangle and the neighborhoold rectangle
	//*   positive: decaying in the positive direction
	//*   negative: decaying in the negative direction
	//* shapeShift: shift parameter that indicates the shift of starting point of 



	this->clear();
    this->isInstantiable = true;

    ifstream ifile(filename.c_str());
    //* Check if open successfully
	if (!ifile){
        cerr << "ERROR: cannot open the file: " << filename << endl;
		return;
	}

	ifile >> this->nWires;
	this->nWirePanels 		= new int[this->nWires];
	this->Cmat				= new float[this->nWires*this->nWires];
	for ( int i=0; i<this->nWires; i++){
		ifile >> this->nWirePanels[i];
	}
	ifile >> this->nPanels;
	this->panels			= new float	[this->nPanels][3][4];
	this->dirs				= new int	[this->nPanels];
	this->areas	 			= new float	[this->nPanels];
	this->indexIncrements 	= new int	[this->nPanels];
	this->basisTypes		= new char	[this->nPanels];
	this->basisDirs			= new int	[this->nPanels];
	this->basisZs			= new float	[this->nPanels];
	this->basisShifts		= new float	[this->nPanels];

	this->nCoefs = 0;

	for (int i=0; i<this->nPanels; i++){
		ifile >> this->basisTypes[i];
		ifile >> this->indexIncrements[i];
		this->nCoefs += this->indexIncrements[i];

		ifile >> this->panels[i][X][MIN];
		ifile >> this->panels[i][X][MAX];
		this->panels[i][X][LENGTH]
			= this->panels[i][X][MAX] - this->panels[i][X][MIN];
		this->panels[i][X][CENTER]
			= (this->panels[i][X][MAX]+this->panels[i][X][MIN])/2;

		ifile >> this->panels[i][Y][MIN];
		ifile >> this->panels[i][Y][MAX];
		this->panels[i][Y][LENGTH]
			= this->panels[i][Y][MAX] - this->panels[i][Y][MIN];
		this->panels[i][Y][CENTER]
			= (this->panels[i][Y][MAX]+this->panels[i][Y][MIN])/2;

		ifile >> this->panels[i][Z][MIN];
		ifile >> this->panels[i][Z][MAX];
		this->panels[i][Z][LENGTH]
			= this->panels[i][Z][MAX] - this->panels[i][Z][MIN];
		this->panels[i][Z][CENTER]
			= (this->panels[i][Z][MAX]+this->panels[i][Z][MIN])/2;

		ifile >> this->dirs[i];
		switch (this->dirs[i]){
		case X:
			this->areas[i] = this->panels[i][Y][LENGTH] * this->panels[i][Z][LENGTH];
			break;
		case Y:
			this->areas[i] = this->panels[i][Z][LENGTH] * this->panels[i][X][LENGTH];
			break;
		case Z:
			this->areas[i] = this->panels[i][X][LENGTH] * this->panels[i][Y][LENGTH];
		}

		ifile >> this->basisDirs[i];
		ifile >> this->basisZs[i];
		ifile >> this->basisShifts[i];
	}
	this->P 	= new float[this->nCoefs*this->nCoefs];
	this->rhs 	= new float[this->nCoefs*this->nWires];
	this->coefs = new float[this->nCoefs*this->nWires];

	// construct nWireCoefs
	this->nWireCoefs = new int[this->nWires];
	int ind = -1;
	for (int i=0; i<this->nWires; i++){
		this->nWireCoefs[i] = 0;
		for (int j=0; j<this->nWirePanels[i]; j++){
			ind++;
			this->nWireCoefs[i] += this->indexIncrements[ind];
		}
	}

	this->isSolved = false;
	this->isLoaded = true;

	ifile.close();

    #ifdef DEBUG_ASPECT_RATIO_VALIDITY
    if( this->isPanelAspectRatioValid() == false ){
        cerr << "ERROR: Panel aspect ratio of CAPLET file is not valid." << endl;
        exit(1);
    }
    #endif
}


void Caplet::saveCmat(const std::string filename){
	if (this->isSolved){
		std::ofstream ofile(filename.c_str());
		if (!ofile.is_open()){
            cerr << "ERROR: cannot write Cmat file: " << filename << endl;
			return;
		}
		//* DOUBLE block was commented before. Double check this.
        switch(this->mode){
        case FAST_GALERKIN:
            for (int i=0; i<this->nWires; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->Cmat[i + this->nWires * j] << " ";
                }
                ofile << endl;
            }
            break;
        case DOUBLE_COLLOCATION:
            for (int i=0; i<this->nWires; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->dCmat[i + this->nWires * j] << " ";
                }
                ofile << endl;
            }
            break;
        default:
            break;
        }

		ofile.close();
	}
}


void Caplet::saveCoefs(const std::string filename){
	if (this->isSolved){
        ofstream ofile(filename.c_str());
		if (!ofile){
            cerr << "ERROR: cannot write the file: " << filename << endl;
			return;
		}
        switch(this->mode){
        case FAST_GALERKIN:
            for (int i=0; i<this->nCoefs; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->coefs[i + this->nCoefs * j]  * 4*pi*epsilon0 << " ";
                }
                ofile << endl;
            }
            break;
        case DOUBLE_COLLOCATION:
            for (int i=0; i<this->nCoefs; i++){
                for (int j=0; j<this->nWires; j++){
                    ofile << this->dcoefs[i + this->nCoefs * j] * 4*pi*epsilon0 << " ";
                }
                ofile << endl;
            }
            break;
        default:
            ;
        }

		ofile.close();
	}
}

void Caplet::extractC(MODE mode){
    #ifdef CAPLET_TIMER
	this->fillingTime = 0;
	this->solvingTime = 0;
	this->totalTime = 0;
    #endif

    #ifdef CAPLET_INIT_ATAN_LOG
    caplet::atan(1);
    caplet::log(1);
    #endif

    //* Setup mode
    if (isInstantiable==true){
        this->mode = FAST_GALERKIN;
    }
    else{
        this->mode = mode;
    }

    //* Init MPI
    MPI::Init();
    int rank = MPI::COMM_WORLD.Get_rank();

    //* Subdivide panels if aspect ratio is too large
    this->modifyPanelAspectRatio();

    if ( rank==0 ){
		std::cout << "Number of conductors        : " << this->nWires << std::endl;
		std::cout << "Number of basis functions   : " << this->nCoefs << std::endl;
        std::cout << "Number of basis shapes      : " << this->nPanels << std::endl;
    }

	for (int iter = 0; iter < N_ITER; iter++){
        switch ( this->mode ){
        case FAST_GALERKIN:
			this->extractCGalerkin();
			break;
        case DOUBLE_COLLOCATION:
            this->extractCCollocationDouble();
            break;
		}
	}
    this->isSolved = true;

    MPI::Finalize();
    if( rank!= 0 ){
		return;
	}


    #ifdef CAPLET_TIMER
    std::cout << "Total extraction time (s)   : " << (this->totalTime)/N_ITER << std::endl;
    std::cout << "  Setup time (s)            : " << (this->fillingTime)/N_ITER << std::endl;
    std::cout << "  Solving time (s)          : " << (this->solvingTime)/N_ITER << std::endl;
    #endif

    //* Print Cmat
    this->printCmat();
	std::cout << endl;
}



//***************************
//*
//* DOUBLE COOLLOCATION MODE
//*
//*
void Caplet::extractCCollocationDouble(){
	if ( this->isLoaded == true ){
        //* Allocate system memory for double precision
        if (MPI::COMM_WORLD.Get_rank()==0){
            this->dP 	= new double[this->nCoefs*this->nCoefs];
            this->drhs 	= new double[this->nCoefs*this->nWires];
            this->dcoefs= new double[this->nCoefs*this->nWires];
            this->dCmat = new double[this->nWires*this->nWires];
        }
        else{
            this->dP 	= new double[1];
            this->drhs 	= new double[1];
            this->dcoefs= new double[1];
            this->dCmat = new double[1];
            return;
        }

	}else{
		std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
	}
    #ifdef CAPLET_TIMER
    this->timeStart = MPI::Wtime();
    #endif

    this->generateCollocationPMatrixDouble();
    this->generateRHSDouble();

    #ifdef CAPLET_TIMER
    this->timeAfterFilling = MPI::Wtime();
    #endif


    //* Solve the system
	int* 	ipiv = new int[this->nCoefs];
	int  	info;

    dgesv_(&this->nCoefs, &this->nWires, this->dP, &this->nCoefs, ipiv,
            this->dcoefs, &this->nCoefs, &info);

	delete[] ipiv;

	char 	transA 	= 't';
	char 	transB 	= 'n';
	double 	alpha 	= 4*pi*epsilon0;
	double 	beta 	= 0.0;

    //* Use matrix-matrix product to compute Cmat from coefs
    dgemm_(&transA, &transB,
            &this->nWires, &this->nWires, &this->nCoefs,
            &alpha, this->drhs, &this->nCoefs,
            this->dcoefs, &this->nCoefs,
            &beta, this->dCmat, &this->nWires);

    #ifdef CAPLET_TIMER
    this->timeAfterSolving = MPI::Wtime();
 	this->fillingTime 	+= this->timeAfterFilling - this->timeStart;
 	this->solvingTime 	+= this->timeAfterSolving - this->timeAfterFilling;
 	this->totalTime		+= this->timeAfterSolving - this->timeStart;
    #endif
}


void Caplet::generateCollocationPMatrixDouble(){

    //* Set dP to zero
    double zero = 0.0;
    int    inc  = 1;
    int    nC   = nCoefs*nCoefs;
    dscal_(&nC, &zero, dP, &inc);

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }
    const int nK = nPanels*nPanels;

    #ifdef CAPLET_OPENMP
        #pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
    #endif
    for ( int k=0; k < nK; k++ ){
        int j = k/nPanels;
        int i = k - j*nPanels;
        dP[ ind[i] + nCoefs*ind[j] ] += calGalerkinPEntry(i,j);
    }

    delete[] ind;
}


double Caplet::calCollocationPEntryDouble(int panel1, int panel2){
    //* Test point : panel1.center
    //  Integration: panel2

    //* A nDim-element array of pointers pointing to a nBit-element array
    float *coord_ptr_1[3][4];
    float *coord_ptr_2[3][4];

    //* Rotate, mirror, and call proper integrals
    switch( dirs[panel2] ){
    case X:
        this->rotateX2Z(coord_ptr_1, panel1);
        this->rotateX2Z(coord_ptr_2, panel2);
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    case Y:
        this->rotateY2Z(coord_ptr_1, panel1);
        this->rotateY2Z(coord_ptr_2, panel2);
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    case Z:
        return calColD( coord_ptr_1, coord_ptr_2);
        break;
    }
    return 0.0; // dummy return
}


void Caplet::generateRHSDouble(){
    int inc = 1;
    int nRHS = this->nWires*this->nCoefs;
    double alpha = 0.0f;
    dscal_(&nRHS, &alpha, this->dcoefs, &inc);

    int rowIndex = -1;
    int columnIndex = 0;
    int coefCounter = 0;
    for ( int i=0; i < this->nPanels; i++ ){
        rowIndex 	+= this->indexIncrements[i];
        coefCounter += this->indexIncrements[i];

        if ( coefCounter > this->nWireCoefs[columnIndex] ){
            coefCounter = 1;
            columnIndex++;
        }

        float (*func)(float x, float w) = &arch;
        switch (this->basisTypes[i]){
        case 'S':
            func = &side;
        case 'A':
            float b = this->panels[i][this->basisDirs[i]][LENGTH];
            float xm = (b)/2;
            float xr = (b)/2;
            int dir = (this->basisDirs[i]+1)%3;
            if ( dir==dirs[i] ){
                dir = (dir+1)%3;
            }
            float w  = this->panels[i][dir][LENGTH];

            //* Gauss quad for computing area integral of shapes
            float val = 0;
            int init_i = 0;
            if ( (gauss_n%2)==1 ){ // odd n
                init_i = 1;
                val += (*gauss::w[gauss_n])[0] * func(xm, w);
            }
            int gauss_n2 = (gauss_n+1)/2;
            for ( int gi = init_i; gi < gauss_n2; gi++ ) {
                float dx = xr * (*gauss::p[gauss_n])[gi];
                val += (*gauss::w[gauss_n])[gi] * ( func(xm+dx, w) + func(xm-dx, w) );
            }
            //* End of Gauss quad

            this->areas[i] = val*xr*w;
        }
        this->dcoefs[ columnIndex*this->nCoefs + rowIndex ] += this->areas[i];
    }
    dcopy_(&nRHS, this->dcoefs, &inc, this->drhs, &inc);
}


//*********************
//*
//* FAST GALERKIN MODE
//*
//*
void Caplet::extractCGalerkin(){

    int rank = MPI::COMM_WORLD.Get_rank();

    if ( this->isLoaded == true ){
        if( rank==0 ){
            this->P 	= new float[this->nCoefs*this->nCoefs];
        }else{
            this->P		= new float[1];
        }

        this->rhs 	= new float[this->nCoefs*this->nWires];
        this->coefs = new float[this->nCoefs*this->nWires];
        this->Cmat  = new float[this->nWires*this->nWires];

    }else{
        std::cerr << "ERROR: structure file is not yet loaded" << std::endl;
    }
    #ifdef CAPLET_TIMER
    this->timeStart = MPI::Wtime();;
    #endif


    #ifdef CAPLET_MPI
    this->generateGalerkinPMatrixMPI();
    #endif

    #ifndef CAPLET_MPI
    this->generateGalerkinPMatrix();
    #endif

    if (MPI::COMM_WORLD.Get_rank()!=0){
        return;
    }
    this->generateRHS();

    #ifdef DEBUG_PRINT_P
    ofstream fout("pmatrix");
    print_matrix(this->P, this->nCoefs, this->nCoefs, "P", 'u', fout);
    fout.close();
    #endif

    #ifdef CAPLET_TIMER
    this->timeAfterFilling = MPI::Wtime();;
    #endif


    //* Solve the system
    int  	info;
    char uplo = 'u';

    //* Query optimal workspace size
    float*	work = new float[1];
    int		lwork = -1;
    int*	ipiv = new int[this->nCoefs];
    ssysv_(&uplo, &nCoefs, &nWires, P, &nCoefs, ipiv, coefs, &nCoefs, work, &lwork, &info);
    lwork = work[0];
    delete[] work;

    //* Solve system using optimal work length
    work = new float[lwork];
    ssysv_(&uplo, &nCoefs, &nWires, P, &nCoefs, ipiv, coefs, &nCoefs, work, &lwork, &info);
    delete[] ipiv;
    delete[] work;


    //* Use matrix-matrix product to compute Cmat from coefs
    char 	transA 	= 't';
    char 	transB 	= 'n';
    float 	alpha 	= 4*pi*epsilon0;
    float 	beta 	= 0.0;
    sgemm_(&transA, &transB,
            &this->nWires, &this->nWires, &this->nCoefs,
            &alpha, this->rhs, &this->nCoefs,
            this->coefs, &this->nCoefs,
            &beta, this->Cmat, &this->nWires);

    #ifdef CAPLET_TIMER
    this->timeAfterSolving = MPI::Wtime();

    this->fillingTime 	+= this->timeAfterFilling - this->timeStart;
    this->solvingTime 	+= this->timeAfterSolving - this->timeAfterFilling;
    this->totalTime		+= this->timeAfterSolving - this->timeStart;
    #endif
}


void Caplet::generateGalerkinPMatrix(){

    float zero = 0.0f;
    int   inc  = 1;
    int   nC   = nCoefs*nCoefs;
    sscal_(&nC, &zero, P, &inc);

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }
    const int nK = nPanels*(nPanels+1)/2;

    #ifdef CAPLET_OPENMP
        #pragma omp parallel for num_threads(CAPLET_OPENMP_NUM_THREADS)
    #endif
    for ( int k=0; k < nK; k++ ){
        int j = int((sqrt(double(1+8*k))-1)/2);
        int i = k - j*(j+1)/2;

        float result = calGalerkinPEntry(i,j);

        if ( (i!=j) && (ind[i]==ind[j]) ){
            P[ ind[i] + nCoefs*ind[j] ] += result*2;
        }else{
            P[ ind[i] + nCoefs*ind[j] ] += result;
        }
    }

    delete[] ind;
}


//* Convert lower triangular matrix (column major) index k to subscript i,j
inline void ltind2sub(int k, int& i, int& j){
    j = int((std::sqrt(double(1+8*k))-1)/2);
    i = k - j*(j+1)/2;
}
void Caplet::generateGalerkinPMatrixMPI(){

    int rank = MPI::COMM_WORLD.Get_rank();
    int numproc = MPI::COMM_WORLD.Get_size();

    int totalK = nPanels*(nPanels+1)/2;
    int nK = totalK/numproc;

    int* startK 	= new int[numproc];
    int* lastK 		= new int[numproc];
    int* startC		= new int[numproc];
    int* lastC		= new int[numproc];

    for( int i=0; i<numproc; i++ ){
        startK[i] = i*nK;
        lastK[i] = (i+1)*nK-1;
    }
    lastK[numproc-1] = totalK-1;

    //* Construct ind_vec from indexIncrements
    int* ind = new int[nPanels];
    ind[0] = 0;
    for ( int i=1; i<nPanels; i++){
        ind[i] = ind[i-1] + indexIncrements[i];
    }

    int maxNC = -1;
    for( int r=0; r<numproc; r++ ){
        int i, startj, lastj;
        ltind2sub(startK[r], i, startj);
        ltind2sub(lastK[r], i,  lastj );
        startC[r] = ind[startj];
        lastC[r]  = ind[lastj];

        int nC = lastC[r]-startC[r]+1;
        if ( maxNC < nC ){
            maxNC = nC;
        }
    }

    float* ptrP;
    float* tempP;
    if ( rank==0 ){
        //* init this->P
        float zero = 0.0f;
        int   inc  = 1;
        int   nP   = nCoefs*nCoefs;
        sscal_(&nP, &zero, P, &inc);
        //* init tempP to cover the largest number of columns among all proccesses.
        int   nTempP = maxNC*nCoefs;
        tempP= new float[nTempP];
        ptrP = this->P;
    }else{
        float zero 	= 0.0f;
        int   inc  	= 1;
        int   nTempP= ( lastC[rank]-startC[rank]+1 )*nCoefs;
        tempP= new float[ nTempP ];
        sscal_(&nTempP, &zero, tempP, &inc);
        ptrP = tempP;
    }

    //* For each k index
    for ( int k = startK[rank] ; k <= lastK[rank] ; k++ ){
        //* Convert k index to i,j subscript
        int i,j;
        ltind2sub(k, i, j);
        //* Compute the P entry
        float result = calGalerkinPEntry(i,j);
        //* Combine rows or columns in place if consecutive panels
        //  belong to the same basis function
        if ( (i!=j) && (ind[i]==ind[j]) ){
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result*2;
        }else{
            ptrP[ ind[i] + nCoefs*(ind[j]-startC[rank]) ] += result;
        }
    }

    //* Combine sub-matrices to rank0 node
    if( rank==0 ){
        MPI::Status status;
        for ( int i=1; i<numproc; i++ ){
            int copylen = (lastC[i]-startC[i]+1)*nCoefs;
            MPI::COMM_WORLD.Recv(tempP, copylen, MPI::FLOAT, i, 0, status);
            float alpha = 1.0f;
            int inc = 1;
            ptrP = P+(startC[i]*nCoefs);
            saxpy_(&copylen, &alpha, tempP, &inc, ptrP, &inc);
        }
    }else{
        int copylen = (lastC[rank]-startC[rank]+1)*nCoefs;
        MPI::COMM_WORLD.Send(tempP, copylen, MPI::FLOAT, 0, 0);
    }


    delete[] startK;
    delete[] lastK;
    delete[] startC;
    delete[] lastC;
    delete[] tempP;
    delete[] ind;

}


float Caplet::calGalerkinPEntry(int panel1, int panel2){

    //* A nDim-element array of pointers pointing to a nBit-element array
    float *coord_ptr_1[3][4];
    float *coord_ptr_2[3][4];

    if ( this->basisTypes[panel1] == 'F' ){
        int temp = panel1;
        panel1 = panel2;
        panel2 = temp;
    }

    //* Select shapes
    shape_t shape1 = selectShape(panel1);
    shape_t shape2 = selectShape(panel2);

    //* Rotate, mirror, and call proper integrals
    switch( dirs[panel1] ){
    case X:
        switch( basisDirs[panel1] ){
        case Y: // Xy
            this->rotateX2Z(coord_ptr_1, panel1);
            this->rotateX2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Xy_X
                switch( basisDirs[panel2] ){
                case Y: 	// Xy_Xy
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z: 	// Xy_Xz
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT: 	// Xy_Xf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Xy_Y
                switch( basisDirs[panel2] ){
                case X:		// Xy_Yx
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xy_Yz
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Xy_Yf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z:
                switch( basisDirs[panel2] ){
                case X:		// Xy_Zx
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y: 	// Xy_Zy
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT: ;
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Z: // Xz
            this->mirrorX2Z(coord_ptr_1, panel1);
            this->mirrorX2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X:
                switch( basisDirs[panel2] ){
                case Y:		// Xz_Xy
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xz_Xz
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Xz_Xf
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y:
                switch( basisDirs[panel2] ){
                case X:		// Xz_Yx
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Xz_Yz
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
                break;
            case Z:
                switch( basisDirs[panel2] ){
                case X:		// Xz_Zx
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Xz_Zy
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT: // Xf-?f
            switch( dirs[panel2] ){
            case X: // Xf_Xf
                this->mirrorX2Z(coord_ptr_1, panel1);
                this->mirrorX2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            case Y: // Xf_Yf
                this->rotateX2Z(coord_ptr_1, panel1);
                this->rotateX2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z: // Xf_Zf
                this->mirrorX2Z(coord_ptr_1, panel1);
                this->mirrorX2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            }
        }
    case Y:
        switch( basisDirs[panel1] ){
        case X: 	// Yx
            this->mirrorY2Z(coord_ptr_1, panel1);
            this->mirrorY2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Yx_X
                switch( basisDirs[panel2] ){
                case Y: 	// Yx_Xy
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yx_Xz
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Yx_Y
                switch( basisDirs[panel2] ){
                case X:		// Yx_Yx
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Yx_Z
                switch( basisDirs[panel2] ){
                case X:		// Yx_Zx
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Yx_Zy
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Z: 	// Yz
            this->rotateY2Z(coord_ptr_1, panel1);
            this->rotateY2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Yz_X
                switch( basisDirs[panel2] ){
                case Y:		// Yz_Xy
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yz_Xz
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Yz_Xf
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Yz_Y
                switch( basisDirs[panel2] ){
                case X:		// Yz_Yx
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Yz_Yz
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Yz_Z
                switch( basisDirs[panel2] ){
                case X:		// Yz_Zx
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Yz_Zy
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT: 	// Yf_?f
            switch( dirs[panel2] ){
            case X:
                this->mirrorY2Z(coord_ptr_1, panel1);
                this->mirrorY2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
                this->mirrorY2Z(coord_ptr_1, panel1);
                this->mirrorY2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            case Z:
                this->rotateY2Z(coord_ptr_1, panel1);
                this->rotateY2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            }
        }
    case Z: // Z
        switch( basisDirs[panel1] ){
        case X:		// Zx
            this->rotateZ2Z(coord_ptr_1, panel1);
            this->rotateZ2Z(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Zx_X
                switch( basisDirs[panel2] ){
                case Y:		// Zx_Xy
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zx_Xz
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:	// Zx_Xf
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y:	// Zx_Y
                switch( basisDirs[panel2] ){
                case X:		// Zx_Yx
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zx_Yz
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z: // Zx_Z
                switch( basisDirs[panel2] ){
                case X:		// Zx_Zx
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Zx_Zy
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case Y:		// Zy
            this->mirrorY2X(coord_ptr_1, panel1);
            this->mirrorY2X(coord_ptr_2, panel2);
            switch( dirs[panel2] ){
            case X: // Zy_X
                switch( basisDirs[panel2] ){
                case Y:		// Zy_Xy
                    return intZXYX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zy_Xz
                    return intZXYZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXYF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Y: // Zy_Y
                switch( basisDirs[panel2] ){
                case X:		// Zy_Yx
                    return intZXXY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Z:		// Zy_Yz
                    return intZXXZ(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXXF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            case Z:	// Zy_Z
                switch( basisDirs[panel2] ){
                case X:		// Zy_Zx
                    return intZXZY(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case Y:		// Zy_Zy
                    return intZXZX(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2, basisZs[panel2], basisShifts[panel2], shape2);
                case FLAT:
                    return intZXZF(	coord_ptr_1, basisZs[panel1], basisShifts[panel1], shape1,
                                    coord_ptr_2);
                }
            }
        case FLAT:	// Zf_?f
            switch( dirs[panel2] ){
            case X:
                this->rotateZ2Z(coord_ptr_1, panel1);
                this->rotateZ2Z(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Y:
                this->mirrorY2X(coord_ptr_1, panel1);
                this->mirrorY2X(coord_ptr_2, panel2);
                return intZFXF(coord_ptr_1, coord_ptr_2);
            case Z:
                this->rotateZ2Z(coord_ptr_1, panel1);
                this->rotateZ2Z(coord_ptr_2, panel2);
                return intZFZF(coord_ptr_1, coord_ptr_2);
            }
        }
    }

    return 0.0; // dummy return
}


void Caplet::generateRHS(){
	int inc = 1;
	int nRHS = this->nWires*this->nCoefs;
	float alpha = 0.0f;
	sscal_(&nRHS, &alpha, this->coefs, &inc);

	int rowIndex = -1;
	int columnIndex = 0;
	int coefCounter = 0;
	for ( int i=0; i < this->nPanels; i++ ){
		rowIndex 	+= this->indexIncrements[i];
		coefCounter += this->indexIncrements[i];

		if ( coefCounter > this->nWireCoefs[columnIndex] ){
			coefCounter = 1;
			columnIndex++;
		}

		float (*func)(float x, float w) = &arch;
		switch (this->basisTypes[i]){
		case 'S':
			func = &side;
		case 'A':
			float b = this->panels[i][this->basisDirs[i]][LENGTH];
			float xm = (b)/2;
			float xr = (b)/2;
			int dir = (this->basisDirs[i]+1)%3;
			if ( dir==dirs[i] ){
				dir = (dir+1)%3;
			}
			float w  = this->panels[i][dir][LENGTH];

            //* Gauss quad for computing area integral of shapes
            float val = 0;
            int init_i = 0;
            if ( (gauss_n%2)==1 ){ // odd n
                init_i = 1;
                val += (*gauss::w[gauss_n])[0] * func(xm, w);
            }
            int gauss_n2 = (gauss_n+1)/2;
            for ( int gi = init_i; gi < gauss_n2; gi++ ) {
                float dx = xr * (*gauss::p[gauss_n])[gi];
                val += (*gauss::w[gauss_n])[gi] * ( func(xm+dx, w) + func(xm-dx, w) );
            }
            //* End of Gauss quad

            this->areas[i] = val*xr*w;
		}

		this->coefs[ columnIndex*this->nCoefs + rowIndex ] += this->areas[i];
	}
	scopy_(&nRHS, this->coefs, &inc, this->rhs, &inc);
}


void Caplet::modifyPanelAspectRatio(){
    //* Aspect ratio of all panels cannot exceed MAX_ASPECT_RAIO
	if ( this->isLoaded == false ){
		std::cerr << "ERROR: not yet load any structure " << std::endl;
		std::cerr << "       Unable to check structure validity" << std::endl;

		std::exit(1);
	}
    const float maxRatio = MAX_ASPECT_RATIO;

	std::list<float> 	panelsTemp[3][2];
	std::list<int>		dirsTemp;
	std::list<float>	areasTemp;
	std::list<int>		indexIncrementsTemp;
	std::list<char>		basisTypesTemp;
	std::list<int>		basisDirsTemp;
	std::list<float>	basisZsTemp;
	std::list<float>	basisShiftsTemp;

	int currentWireIndex = 0;
	int currentWirePanelCount = 0;
	int nCurrentWirePanels = this->nWirePanels[0];
	int nFinalPanels = this->nPanels;

	for ( int i=0; i<this->nPanels; i++ ){

		currentWirePanelCount ++;
		if ( currentWirePanelCount > this->nWirePanels[currentWireIndex] ){
			currentWirePanelCount = 1;
			this->nWirePanels[currentWireIndex] = nCurrentWirePanels;
			currentWireIndex ++;
			nCurrentWirePanels = this->nWirePanels[currentWireIndex];
		}
		float ba = this->panels[i][ (this->dirs[i]+1)%3 ][LENGTH] / this->panels[i][ (this->dirs[i]+2)%3 ][LENGTH];
        if (  ( ba > maxRatio )  ||  ( 1/ba > maxRatio )  ){
            //* need to split to reduce aspact ratio of a panel
			int nSplit = 2;
			while ( ( ba/nSplit > maxRatio )  ||  ( 1/ba/nSplit > maxRatio ) ){
				nSplit++;
			}
			int   splitDir = ( ( ba>1 ) ? ( (dirs[i]+1)%3 ) : ( (dirs[i]+2)%3 ) );
			float currentCoord = this->panels[i][splitDir][0];
			float sublength = this->panels[i][splitDir][LENGTH] / nSplit;

			for ( int j=0; j<nSplit; j++){
				panelsTemp[splitDir][0].push_back( currentCoord );
				currentCoord += sublength;
				panelsTemp[splitDir][1].push_back( currentCoord );

				panelsTemp[ (splitDir+1)%3 ][0].push_back( this->panels[i][ (splitDir+1)%3 ][0] );
				panelsTemp[ (splitDir+1)%3 ][1].push_back( this->panels[i][ (splitDir+1)%3 ][1] );
				panelsTemp[ (splitDir+2)%3 ][0].push_back( this->panels[i][ (splitDir+2)%3 ][0] );
				panelsTemp[ (splitDir+2)%3 ][1].push_back( this->panels[i][ (splitDir+2)%3 ][1] );

				dirsTemp.push_back(this->dirs[i]);
				areasTemp.push_back(this->areas[i]/nSplit);
				if ( j==0 ){
                    //* The first sub-panel after split. Needs index increment.
					indexIncrementsTemp.push_back(this->indexIncrements[i]);
				}else{
                    //* The rest of sub-panel after split. No index increment needed.
					indexIncrementsTemp.push_back(0);
				}
				basisTypesTemp.push_back(this->basisTypes[i]);
				basisDirsTemp.push_back(this->basisDirs[i]);
				basisZsTemp.push_back(this->basisZs[i]);

				if ( this->basisDirs[i] == splitDir ){
					if ( this->basisZs[i] > 0 ){ // decaying in the positive direction
						basisShiftsTemp.push_back(this->basisShifts[i] + j*sublength);
					}else{ // decaying in the negative direction
						basisShiftsTemp.push_back(this->basisShifts[i] + (nSplit-1-j)*sublength);
					}
				}else{
					basisShiftsTemp.push_back(this->basisShifts[i]);
				}

				if ( j!=0 ){
					nFinalPanels ++;
					nCurrentWirePanels ++;
				}
			}
        }
        else{
            //* no need to split. simply push back
			panelsTemp[X][0].push_back(this->panels[i][X][0]);
			panelsTemp[X][1].push_back(this->panels[i][X][1]);
			panelsTemp[Y][0].push_back(this->panels[i][Y][0]);
			panelsTemp[Y][1].push_back(this->panels[i][Y][1]);
			panelsTemp[Z][0].push_back(this->panels[i][Z][0]);
			panelsTemp[Z][1].push_back(this->panels[i][Z][1]);

			dirsTemp.push_back(this->dirs[i]);
			areasTemp.push_back(this->areas[i]);
			indexIncrementsTemp.push_back(this->indexIncrements[i]);
			basisTypesTemp.push_back(this->basisTypes[i]);
			basisDirsTemp.push_back(this->basisDirs[i]);
			basisZsTemp.push_back(this->basisZs[i]);
			basisShiftsTemp.push_back(this->basisShifts[i]);
		}
	}
	this->nPanels = nFinalPanels;
	this->nWirePanels[currentWireIndex] = nCurrentWirePanels;

    //* Rebuild panel information based on new panel list with split panels
	delete[] this->panels;
	delete[] this->dirs;
    delete[] this->indexIncrements;
	delete[] this->basisTypes;
	delete[] this->basisDirs;
	delete[] this->basisZs;
	delete[] this->basisShifts;

	this->panels			= new float	[this->nPanels][3][4];
	this->dirs				= new int	[this->nPanels];
	this->areas	 			= new float	[this->nPanels];
	this->indexIncrements 	= new int	[this->nPanels];
	this->basisTypes		= new char	[this->nPanels];
	this->basisDirs			= new int	[this->nPanels];
	this->basisZs			= new float	[this->nPanels];
	this->basisShifts		= new float	[this->nPanels];

	for ( int i=0; i<this->nPanels; i++ ){
		this->panels[i][X][0] = panelsTemp[X][0].front(); panelsTemp[X][0].pop_front();
		this->panels[i][X][1] = panelsTemp[X][1].front(); panelsTemp[X][1].pop_front();
		this->panels[i][Y][0] = panelsTemp[Y][0].front(); panelsTemp[Y][0].pop_front();
		this->panels[i][Y][1] = panelsTemp[Y][1].front(); panelsTemp[Y][1].pop_front();
		this->panels[i][Z][0] = panelsTemp[Z][0].front(); panelsTemp[Z][0].pop_front();
		this->panels[i][Z][1] = panelsTemp[Z][1].front(); panelsTemp[Z][1].pop_front();

		this->panels[i][X][LENGTH] = this->panels[i][X][1] - this->panels[i][X][0];
		this->panels[i][Y][LENGTH] = this->panels[i][Y][1] - this->panels[i][Y][0];
		this->panels[i][Z][LENGTH] = this->panels[i][Z][1] - this->panels[i][Z][0];

		this->panels[i][X][CENTER] = (this->panels[i][X][1] + this->panels[i][X][0])/2;
		this->panels[i][Y][CENTER] = (this->panels[i][Y][1] + this->panels[i][Y][0])/2;
		this->panels[i][Z][CENTER] = (this->panels[i][Z][1] + this->panels[i][Z][0])/2;

		this->dirs[i] = dirsTemp.front(); dirsTemp.pop_front();
		this->areas[i] = areasTemp.front(); areasTemp.pop_front();

		this->indexIncrements[i] = indexIncrementsTemp.front(); indexIncrementsTemp.pop_front();

		this->basisTypes[i] = basisTypesTemp.front(); basisTypesTemp.pop_front();
		this->basisDirs[i] = basisDirsTemp.front(); basisDirsTemp.pop_front();
		this->basisZs[i] = basisZsTemp.front(); basisZsTemp.pop_front();
		this->basisShifts[i] = basisShiftsTemp.front(); basisShiftsTemp.pop_front();
	}
}


shape_t Caplet::selectShape(int panel){
	switch ( this->basisTypes[panel] ){
	case 'A':
		return &arch;
		break;
	case 'S':
		return &side;
		break;
	default:
		return 0;
	}
}


int  Caplet::getNPanels() const{
	return this->nPanels;
}


int  Caplet::getNCoefs() const{
	return this->nCoefs;
}


int  Caplet::getSizeP() const{
	return this->nCoefs * this->nCoefs;
}


int  Caplet::getSizeCoefs() const{
	return this->nCoefs * this->nWires;
}


int  Caplet::getSizeCmat() const{
	return this->nWires * this->nWires;
}


const float* const Caplet::getCmat() const{
	return this->Cmat;
}


float Caplet::compareCmatError(const float* const cmatRef, ERROR_REF option) const{
	using std::abs;

	float* const cmatError = new float[this->getSizeCmat()];
	int n = this->nWires;
	float maxError = 0.0f;
	int maxErrorI = 0;
	int maxErrorJ = 0;

	for ( int i=0; i<n; i++){
		for ( int j=0; j<n; j++ ){

			float errRef;
			switch (option){
			case DIAGONAL:
				errRef = (cmatRef[ i + n*i ]);
				break;
			case SELF:
				errRef = (cmatRef[ i + n*j ]);
			}
			float thisError = abs( (this->Cmat[ i + n*j] - cmatRef[ i + n*j ])/errRef * 100 );
			cmatError[ i + n*j ] = thisError;

			if ( maxError < thisError ){
				maxError  = thisError;
				maxErrorI = i;
				maxErrorJ = j;
			}
		}
	}
	print_matrix(cmatError, n, n, "Cmat Error (%)", 'a');

	std::cout << "Max error = " << std::fixed << maxError << "% at (" << maxErrorI << "," << maxErrorJ << ")" << std::endl;

	switch( option ){
	case DIAGONAL:
		std::cout << "    (w.r.t. the row diagonal)" << std::endl;
		break;
	case SELF:
		std::cout << "    (w.r.t. itself)" << std::endl;
	}
	std::cout << std::endl;

	delete [] cmatError;
	return maxError;
}


float Caplet::compareCmatError(const Caplet *const caplet, ERROR_REF option) const{
	const float* const cmatRef = caplet->getCmat();
	return this->compareCmatError(cmatRef, option);
}


float Caplet::compareCmatError(const std::string filename, ERROR_REF option) const{
	const int n = this->nWires;

	float* cmatRef = new float[n*n];

	std::ifstream ifile(filename.c_str());
	if (!ifile){
		std::cerr << "ERROR: cannot open the file: " << filename << std::endl;
		std::exit(1);
	}

	std::string		lineTemp;
	for ( int i=0; i<n; i++ ){
		getline(ifile, lineTemp);
		std::stringstream	stringTokenizer(lineTemp);
		for (int j=0; j<n; j++){
			stringTokenizer >> cmatRef[ i + n*j ];
		}
	}

	float maxError = this->compareCmatError(cmatRef, option);
	delete [] cmatRef;

	return maxError;
}



//******************
//*
//* PRINT UTILITIES
//*
//*
void Caplet::printP(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->P, this->nCoefs, this->nCoefs, "P(float)");
            break;
        case DOUBLE_COLLOCATION:
            print_matrix(this->dP, this->nCoefs, this->nCoefs, "P(double)");
            break;
        default:
            ;
        }
    }
}


void Caplet::printCoefs(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->coefs, this->nCoefs, this->nWires, "coefs(float)");
            break;
        case DOUBLE_COLLOCATION:
            print_matrix(this->dcoefs, this->nCoefs, this->nWires, "coefs(double)");
            break;
        default:
            ;
        }
    }
}


void Caplet::printRHS(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->rhs, this->nCoefs, this->nWires, "rhs(float)");
            break;
        case DOUBLE_COLLOCATION:
            print_matrix(this->drhs, this->nCoefs, this->nWires, "rhs(double)");
            break;
        default:
            ;
        }
    }
}


void Caplet::printCmat(){
    if ( this->isLoaded && this->isSolved ){
        switch(this->mode){
        case FAST_GALERKIN:
            print_matrix(this->Cmat, this->nWires, this->nWires, "Cmat(float)");
            break;
        case DOUBLE_COLLOCATION:
            print_matrix(this->dCmat, this->nWires, this->nWires, "Cmat(double)");
            break;
        default:
            ;
        }
    }
}


void Caplet::printPanel(int index){
    if (isLoaded == false){
        return;
    }
    if (index >= nPanels){
        return;
    }
    for (int i=0; i<3; ++i){
        for (int j=0; j<4; ++j){
            cout << setw(12) << panels[index][i][j];
        }
        cout << endl;
    }
}



//******************
//*
//* DEBUG UTILITIES
//*
//*
bool Caplet::isPanelAspectRatioValid(){
    if ( this->isLoaded == false ){
        std::cerr << "ERROR: not yet load any structure " << std::endl;
        std::cerr << "       Unable to check structure validity" << std::endl;

        return false;
    }
    const float maxRatio = 51;

    for ( int i=0; i< this->nPanels; i++ ){
        float ratio = panels[i][(dirs[i]+1)%3][LENGTH] / panels[i][(dirs[i]+2)%3][LENGTH];

        if ( ratio > maxRatio || 1/ratio > maxRatio ){
            std::cerr << "ERROR: panel ratio check fails" << std::endl;
            std::cerr << "    Panel: " << i << std::endl;
            std::cerr << "    x:" << panels[i][X][0] << ", " << panels[i][X][1] << ", " << panels[i][X][2] << ", " << panels[i][X][3] << std::endl;
            std::cerr << "    y:" << panels[i][Y][0] << ", " << panels[i][Y][1] << ", " << panels[i][Y][2] << ", " << panels[i][Y][3] << std::endl;
            std::cerr << "    z:" << panels[i][Z][0] << ", " << panels[i][Z][1] << ", " << panels[i][Z][2] << ", " << panels[i][Z][3] << std::endl;

            return false;
        }
    }
    return true;
}

} // end of namespace caplet
