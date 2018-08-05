#include "mex.h"
/*
 * RG.cpp
 * Sets each voxel of b to the corresponding voxel of either bA or bB
 * using a multi-seeded safest-first region growing scheme
 *
 * This is a MEX-file for MATLAB.
 */

#include <queue>
#include <cmath>

using namespace std;

unsigned char* status;	// status image; 0=undetermined, 1=bA, 2=bB
float* bA_arg;			// image with phase of bA
float* bB_arg;			// image with phase of bB
float* mw;				// magnitude weight image
double dx, dy, dz;		// voxel size
long nx, ny, nz;		// image size
int numQueues = 100;	// number of FIFO queues
queue<long>* q;			// vector of FIFO queues
float maxPrio;			// maximum Q value, as defined by eqs 9 & 10

// picks an unchecked voxel from the safest (highest Q) nonempty queue
long popQueue() {
	long next;
	for (int k=numQueues-1; k>=0; k--)
		while (!q[k].empty()) {
			next = q[k].front();
			q[k].pop();
			if (status[next]==0) // assert undetermined
				return next;
		}
	return -1; // if all queues are empty
}

// get indexes of the 6-neighbourhood of voxel 'index' and put into Ngb and number into numNgb
void getNgb(const long index, int& numNgb, long* Ngb) {
	numNgb = 0;
	long x = index%nx;
	if (x>0) {Ngb[numNgb] = index-1; numNgb++;}
	if (x<nx-1) {Ngb[numNgb] = index+1; numNgb++;}
	long index2=(index-x)/nx;
	long y = index2%ny;
	if (y>0) {Ngb[numNgb] = index-nx; numNgb++;}
	if (y<ny-1) {Ngb[numNgb] = index+nx; numNgb++;}
	long z = (index2-y)/ny;
	if (z>0) {Ngb[numNgb] = index-nx*ny; numNgb++;}
	if (z<nz-1) {Ngb[numNgb] = index+nx*ny; numNgb++;}
}

// calculate current value Q of voxel according to eqs 9 & 10, and put on corresponding queue
void examine(const long voxel) {
	int numNgb;
	long* Ngb = new long[6];
	getNgb(voxel,numNgb,Ngb);
	
	float CA=0;
	float CB=0;
	
	float bA=bA_arg[voxel];
	float bB=bB_arg[voxel];
	float bk;
	double dk;
	for (int n=0; n<numNgb; n++)
		if (status[Ngb[n]]>0) {
			long indexDiff = abs(Ngb[n]-voxel);
			if (indexDiff==1) dk=dx;		//neighbours along x direction
			else if(indexDiff==nx) dk=dy;	//neighbours along y direction
			else dk=dz;						//neighbours along z direction
			
			float mwk=mw[Ngb[n]];
			status[Ngb[n]]==1 ? bk = bA_arg[Ngb[n]] : bk = bB_arg[Ngb[n]];
			CA+=mwk/dk*cos(bk-bA);
			CB+=mwk/dk*cos(bk-bB);
		}
	
	int Q = int(abs(CA-CB)/maxPrio*float(numQueues-1));
	q[Q].push(voxel);
	
	delete[] Ngb;
}

// identify the b value of voxel, according to eq 9, then examine its neighbours
bool determine(const long voxel) {
	if (voxel<0) return false;
	int numNgb;
	long* Ngb = new long[6];
	getNgb(voxel,numNgb,Ngb);
	
	float C1 = 0;
	float C2 = 0;
	
	float bA=bA_arg[voxel];
	float bB=bB_arg[voxel];
	float bk;
	double dk;
	
    int n;
	for (n=0; n<numNgb; n++) 
		if (status[Ngb[n]]>0) {
			long indexDiff = abs(Ngb[n]-voxel);
			if (indexDiff==1) dk=dx;
			else if(indexDiff==nx) dk=dy;
			else dk=dz;
			
			float mwk=mw[Ngb[n]];
			status[Ngb[n]]==1 ? bk = bA_arg[Ngb[n]] : bk = bB_arg[Ngb[n]];
			C1+=mwk/dk*cos(bk-bA);
			C2+=mwk/dk*cos(bk-bB);
		}
	C1>C2 ? status[voxel]=1 : status[voxel]=2;
	
	for (n=0; n<numNgb; n++)
		if (status[Ngb[n]]==0) examine(Ngb[n]);
	
 	delete[] Ngb;
	return true;
}

void RegionGrowing(float* b_arg ) {
	q = new queue<long>[numQueues];
	nz==1 ? maxPrio = 4*(dx+dy) : maxPrio = 4*(dx+dy+dz);
	
	int numNgb;
	long* Ngb = new long[6];
	
	// examine neighbours of the seeds
    long i;
	for (i=0; i<nx*ny*nz; i++)
		if (status[i]>0) {
			getNgb(i,numNgb,Ngb);
			for (int n=0; n<numNgb; n++)
				if (status[Ngb[n]]==0) examine(Ngb[n]);
		}
	delete[] Ngb;
	
	// identify voxels and examine their neighbours until all queues are empty
	while (determine(popQueue())) {}
  	delete[] q;
	
	// set output data according to status image
	for (i=0; i<nx*ny*nz; i++)
		status[i]==1 ? b_arg[i] = bA_arg[i] : b_arg[i] = bB_arg[i];
}

// MATLAB Gateway function
void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[] ) {
    status = (unsigned char*) mxGetData(prhs[0]);						// 1st input: status image
    bA_arg = (float*) mxGetData(prhs[1]);								// 2nd input: bA candidate phase image
    bB_arg = (float*) mxGetData(prhs[2]);								// 3rd input: bB candidate phase image
    mw = (float*) mxGetData(prhs[3]);									// 4th input: magnitude weight image
    double* voxelsize = mxGetPr(prhs[4]);								// 5th input: voxel size
	dx = voxelsize[0]; dy = voxelsize[1]; dz = voxelsize[2];
    
    const int* size = mxGetDimensions(prhs[0]);							// get image size
    int n = mxGetNumberOfDimensions(prhs[0]);
    nx = size[0]; ny = size[1]; (n>2) ? nz = size[2] : nz = 1;
    
    plhs[0] = mxCreateNumericArray(n, size, mxSINGLE_CLASS, mxREAL);	// create output array
    float* b_arg = (float*) mxGetData(plhs[0]);							// get data pointer to output
    
    RegionGrowing(b_arg);												// do the region growing
}
