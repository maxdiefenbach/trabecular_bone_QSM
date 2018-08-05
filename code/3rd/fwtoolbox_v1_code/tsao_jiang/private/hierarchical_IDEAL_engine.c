/*
 * Calculation engine for Hierarchical IDEAL
 * This version sequentially determines B0, then T2*, which is faster.
 *
 * Description: Fat-water separation by hierarchical decomposition and
 *              direct estimation of phase offset to locate signal null.
 *
 * Tsao J, Jiang Y. Hierarchical IDEAL: robust water?fat separation at high
 * field by multiresolution field map estimation. In: Proceedings of the 18th
 * Annual Meeting of ISMRM, Toronto, ON, Canada, 2008. p 653
 *
 * function [phasemap,t2correctmap] = hierarchical_IDEAL_engine(img,minsize,DegreeRange,DegreeOffset,MaxDepth);
 *  img           : Array of x * y * z * 3
 *  minsize       : Minimum size to subdivide until. default = 0.1
 *  DegreeRange,
 *  DegreeOffset  : Check from angleoffset-anglerange to angleoffset+anglerange
 *                  (in degrees)
 *  MaxDepth      : maximum decomposition depth
 *
 * AUTHOR: Jeffrey Tsao
 *
 * VERSION HISTORY
 *  2011.09.19 Added CompareString_IgnoreCase
 *  2008.06.03 Added MaxDepth argument
 *  2008.06.05 Bug fix with normalization of phasemap, reorg of order between
 *             decomposition & sending back results from current level.
 *  2008.06.06 Added T2* correction.
 *
 */

#include "mex.h"
#include <math.h>
#include <stdlib.h>     /* rand, malloc */
#ifndef SQR
#define SQR(x)          ((x)*(x))
#endif
#ifndef M_PI
#define M_PI    3.1415926536
#endif
#define DegreeRangeContraction  0.8 /* 0.5 */
#define DegreeRangeMax          180  /* 55  less than 1/3 of 180deg */
#define DegreeRangeMax_depth    2
#define DegreeRangeMin          10  /* less than 1/3 of 180deg */

void IDEAL_engine_core(void *img_r, void *img_i, int DataType,
                       double *phasemap_r,double *phasemap_i,
                       double *t2correctmap_r,double *t2correctmap_i,
                       int NumDims, int *ArraySize,
                       int depth, int MaxDepth, /* hierarchical depth */
                       int MinSizeX,    /* Min size along each dimension */
                       int MinSizeY,    /* Min size along each dimension */
                       int MinSizeZ,    /* Min size along each dimension */
                       int FromX, int ToX,
                       int FromY, int ToY,
                       int FromZ, int ToZ,
                       double DegreeRange, double DegreeOffset, double T2Correct, double *NoiseVar);
int CompareString_IgnoreCase(char *s1, char *s2);

/*----------------------------------------------------------------*/
/* Simplex optimization */
#define SIMPLEX_UNBIASED_DIRECTION        0
typedef struct {
	double Vertex[4][3];
	double FuncValue[4];
	double Centroid[4];
	double MinParam[3];
	double MaxParam[3];
	int NumEvals;
} SimplexDataType;
void SimplexSetup(double          *ClientData,               /*::: Input :::*/
                  SimplexDataType *Simplex, /*::: Input & Output :::*/
                  int NumDim,                            /*::: Input WARNING: NumDim <= 3 :::*/
                  double           *StartParam, /*::: Input ::: */
                  double           *ParamStepSize, /*::: Input :::*/
                  double           *MinParam, /*::: Input :::*/
                  double           *MaxParam); /*::: Input :::*/
int SimplexMinimize(double *ClientData,                  /*::: Input :::*/
                    SimplexDataType *Simplex, /*::: Input & Output :::*/
                    int NumDim,                          /*::: Input :::*/
                    double           *MinParamStepSize, /*::: Input :::*/
                    int MaxNumEvals);       /*::: Input :::*/
/*----------------------------------------------------------------*/

double EstimateNoise(void *img_r, void *img_i, int DataType, int NumDims, int *ArraySize);
void showhelp() {
	printf("Expected input arguments:\n");
	printf(" 1. Image - array of (x,y,z,3) of floats or doubles.\n");
	printf(" 2. (optional) MinSize - minimum fractional size to stop sub-divisions.\n");
	printf(" 3. (optional) DegreeRange - angle range (degrees) to search.\n");
	printf(" 4. (optional) DegreeOffset - starting angle (degrees) to search.\n");
	printf(" 5. (optional) MaxDepth - maximum number of hierarchical subdivisions.\n");
}
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
	/*----------------------------------------------------------------*/
	mxArray *tmpMxArray;
	mxArray *Output_PhaseMap;
	mxArray *Output_T2CorrectMap;
	double *pr;
	double DegreeRange, DegreeOffset, T2Correct, MinSize;
	void *img_r, *img_i;
	int DataType;
	double *phasemap_r, *phasemap_i;
	double *t2correctmap_r, *t2correctmap_i;
	int *ArraySize;
	int MaxDepth;
	int MinSizeX, MinSizeY, MinSizeZ;
	double NoiseVar=0.0;

	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	if(nrhs < 1) {
		showhelp();
		mexErrMsgTxt("Incorrect no. input arguments.");
		return;
	}
	if (nrhs == 1 && mxIsChar((mxArray *)prhs[0])) {
		char str[101];
		if (mxGetString(prhs[0], str, 100)) {
			showhelp();
			mexErrMsgTxt("Incorrect input arguments.");
			return;
		}
		if (CompareString_IgnoreCase(str,"name")==0 || CompareString_IgnoreCase(str,"title")==0 || CompareString_IgnoreCase(str,"algorithm")==0) {
			const char *name = "Hierarchical IDEAL with T2* correction";
			if (nlhs>=1) plhs[0] = mxCreateString(name); else printf(name);
			return;
		}
		if (CompareString_IgnoreCase(str,"authors")==0 || CompareString_IgnoreCase(str,"author")==0 || CompareString_IgnoreCase(str,"programmers")==0 || CompareString_IgnoreCase(str,"programmer")==0) {
			const char *name = "Jeffrey Tsao, Yun Jiang";
			if (nlhs>=1) plhs[0] = mxCreateString(name); else printf(name);
			return;
		}
		if (CompareString_IgnoreCase(str,"version")==0) {
			const char *name = "2011.10.17";
			if (nlhs>=1) plhs[0] = mxCreateString(name); else printf(name);
			return;
		}
		printf("Unknown command: %s",str);
		return;
	}

	/*----------------------------------------------------------------*/
	/* Input parameter 1: Image */
	{
		tmpMxArray = (mxArray *)prhs[0];
		img_r = (void *)mxGetPr(tmpMxArray);
		img_i = (void *)(mxIsComplex(tmpMxArray) ? mxGetPi(tmpMxArray) : NULL);
		if (mxGetNumberOfDimensions(tmpMxArray)!=4) {
			mexErrMsgTxt("img (1st parameter) should have 4 dimensions only (x,y,z,3).");
			return;
		}
		ArraySize = (int *)mxGetDimensions(tmpMxArray);
		if (ArraySize[3]!=3) {
			mexErrMsgTxt("img (1st parameter) should have a 4th dimensions of size 3 only.");
			return;
		}
		switch (mxGetClassID(tmpMxArray)) {
		case mxSINGLE_CLASS: DataType = 1; break;
		case mxDOUBLE_CLASS: DataType = 2; break;
		default:             DataType = -1; break;
		}
		if (DataType<0) {
			mexErrMsgTxt("img (1st parameter) should be either floats or doubles.");
			return;
		}
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Input parameter 2 (optional): MinSize */
	while(1) {
		if (nrhs>=2) {
			tmpMxArray = (mxArray *)prhs[1];
			if (!mxIsEmpty(tmpMxArray)) {
				if (mxGetNumberOfElements(tmpMxArray)!=1
				    ||           !mxIsDouble(tmpMxArray)
				    ||           mxIsComplex(tmpMxArray)) {
					mexErrMsgTxt("MinSize (2nd parameter) should have a real-valued scalar.");
					return;
				}
				pr = mxGetPr(tmpMxArray);
				MinSize = pr[0];
				break;
			}
		}
		MinSize = 0.1;
		break;
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Input parameter 3 (optional): DegreeRange */
	while(1) {
		if (nrhs>=3) {
			tmpMxArray = (mxArray *)prhs[2];
			if (!mxIsEmpty(tmpMxArray)) {
				if (mxGetNumberOfElements(tmpMxArray)!=1
				    ||           !mxIsDouble(tmpMxArray)
				    ||           mxIsComplex(tmpMxArray)) {
					mexErrMsgTxt("DegreeRange (3rd parameter) should have a real-valued scalar.");
					return;
				}
				pr = mxGetPr(tmpMxArray);
				DegreeRange = pr[0];
				break;
			}
		}
		DegreeRange = 180.0;
		break;
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Input parameter 4 (optional): DegreeOffset */
	while(1) {
		if (nrhs>=4) {
			tmpMxArray = (mxArray *)prhs[3];
			if (!mxIsEmpty(tmpMxArray)) {
				if (mxGetNumberOfElements(tmpMxArray)!=1
				    ||           !mxIsDouble(tmpMxArray)
				    ||           mxIsComplex(tmpMxArray)) {
					mexErrMsgTxt("DegreeOffset (4th parameter) should have a real-valued scalar.");
					return;
				}
				pr = mxGetPr(tmpMxArray);
				DegreeOffset = pr[0];
				break;
			}
		}
		DegreeOffset = 0.0;
		break;
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Input parameter 5 (optional): MaxDepth */
	while(1) {
		if (nrhs>=5) {
			tmpMxArray = (mxArray *)prhs[4];
			if (!mxIsEmpty(tmpMxArray)) {
				if (mxGetNumberOfElements(tmpMxArray)!=1
				    ||           !mxIsDouble(tmpMxArray)
				    ||           mxIsComplex(tmpMxArray)) {
					mexErrMsgTxt("MaxDepth (5th parameter) should have a real-valued scalar.");
					return;
				}
				pr = mxGetPr(tmpMxArray);
				MaxDepth = pr[0];
				break;
			}
		}
		MaxDepth = -1; /* Default: no depth limit */
		break;
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Create output buffer */
	Output_PhaseMap = mxCreateDoubleMatrix(ArraySize[0]*ArraySize[1]*ArraySize[2], 1, mxCOMPLEX);
	phasemap_r = mxGetPr(Output_PhaseMap);
	phasemap_i = mxGetPi(Output_PhaseMap);
	if(mxSetDimensions(Output_PhaseMap, ArraySize, 3)) {
		mxDestroyArray(Output_PhaseMap);
		mexErrMsgTxt("Failed to set the size of output");
		return;
	}
	{ int i;
	  for(i=ArraySize[0]*ArraySize[1]*ArraySize[2]-1; i>=0; i--) phasemap_r[i] = phasemap_i[i] = 0.0; }

	/* Create output buffer */
	Output_T2CorrectMap = mxCreateDoubleMatrix(ArraySize[0]*ArraySize[1]*ArraySize[2], 1, mxCOMPLEX);
	t2correctmap_r = mxGetPr(Output_T2CorrectMap);
	t2correctmap_i = mxGetPi(Output_T2CorrectMap);
	if(mxSetDimensions(Output_T2CorrectMap, ArraySize, 3)) {
		mxDestroyArray(Output_T2CorrectMap);
		mexErrMsgTxt("Failed to set the size of output");
		return;
	}
	{ int i;
	  for(i=ArraySize[0]*ArraySize[1]*ArraySize[2]-1; i>=0; i--) t2correctmap_r[i] = t2correctmap_i[i] = 0.0; }

	MinSizeX = (MinSize<=1 ? ceil(ArraySize[0]*MinSize) : MinSize);
	MinSizeY = (MinSize<=1 ? ceil(ArraySize[1]*MinSize) : MinSize);
	MinSizeZ = (MinSize<=1 ? ceil(ArraySize[2]*MinSize) : MinSize);
	if (MinSizeX<1) MinSizeX=1;
	if (MinSizeY<1) MinSizeY=1;
	if (MinSizeZ<1) MinSizeZ=1;

	T2Correct = 1.0;
	/* NoiseVar = EstimateNoise(img_r,img_i,DataType,4,ArraySize); */ /* Take off noise estimate */

/*  printf("%d %d %d %d\n",ArraySize[0],ArraySize[1],ArraySize[2],ArraySize[3]);
 */
	/*----------------------------------------------------------------*/
	IDEAL_engine_core(img_r,img_i,DataType,
	                  phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
	                  4,ArraySize, /* NumDims, ArraySize */
	                  0,MaxDepth,    /* depth, MaxDepth */
	                  MinSizeX,MinSizeY,MinSizeZ, /* MinSize */
	                  0, ArraySize[0]-1,
	                  0, ArraySize[1]-1,
	                  0, ArraySize[2]-1,
	                  DegreeRange, DegreeOffset, T2Correct, &NoiseVar);
	if (nlhs>=1) plhs[0] = Output_PhaseMap;
	if (nlhs>=2) plhs[1] = Output_T2CorrectMap;

	return;
}

#define COMPLEXCONJ(rr,ri,ar,ai)                {rr =(ar);        ri =-(ai); }
#define COMPLEXCONJ_ACCUM(rr,ri,ar,ai)          {rr+=(ar);        ri-= (ai); }
#define COMPLEXADD(rr,ri,ar,ai,br,bi)           {rr =((ar)+(br)); ri =((ai)+(bi)); }
#define COMPLEXADD_ACCUM(rr,ri,ar,ai,br,bi)     {rr+=((ar)+(br)); ri+=((ai)+(bi)); }
#define COMPLEXMULT(rr,ri,ar,ai,br,bi)          {rr =((ar)*(br)-(ai)*(bi)); ri =((ar)*(bi)+(ai)*(br)); }
#define COMPLEXMULT_r(ar,ai,br,bi)              ((ar)*(br)-(ai)*(bi))
#define COMPLEXMULT_i(ar,ai,br,bi)              ((ar)*(bi)+(ai)*(br))
#define COMPLEXMULT_ACCUM(rr,ri,ar,ai,br,bi)    {rr+=((ar)*(br)-(ai)*(bi)); ri+=((ar)*(bi)+(ai)*(br)); }
#define COMPLEXREALMULT(rr,ri,ar,ai,br)         {rr =((ar)*(br)); ri =((ai)*(br)); }
#define COMPLEXREALMULT_ACCUM(rr,ri,ar,ai,br)   {rr+=((ar)*(br)); ri+=((ai)*(br)); }
#define COMPLEXCONJ_INLINE(ar,ai)               {ai= -(ai); }
#define COMPLEXADD_INLINE(ar,ai,br,bi)          {ar+=(br); ai+=(bi); }
#define COMPLEXMULT_INLINE(ar,ai,br,bi)         {double origar,origai; origar=(ar); origai=(ai); ar=((origar)*(br)-(origai)*(bi)); ai=((origar)*(bi)+(origai)*(br)); }
#define EXP_I(rr,ri,radians)                    {rr=cos(radians); ri=sin(radians); }
#define ABS2(ar,ai)                             ((ar)*(ar)+(ai)*(ai))
void IDEAL_engine_core(void *img_r, void *img_i, int DataType,
                       double *phasemap_r,double *phasemap_i,
                       double *t2correctmap_r, double *t2correctmap_i,
                       int NumDims, int *ArraySize,
                       int depth, int MaxDepth, /* hierarchical depth */
                       int MinSizeX,int MinSizeY,int MinSizeZ,  /* Min size along each dimension */
                       int FromX, int ToX,
                       int FromY, int ToY,
                       int FromZ, int ToZ,
                       double DegreeRange, double DegreeOffset,
                       double t2correct, double *NoiseVar)
{
	double ClientData[9], Param[2];
	double *WeightZ, *WeightY, *WeightX;
	double NoiseMargin = 2.0;
	int Optimize = 1;

	/*----------------------------------------------------------------*/
	/* Find best angle */
	{
		int NumPixels;
		double aA, aB_r, aB_i, aC_r, aC_i, bB, bC_r, bC_i, cC;

		NumPixels = ArraySize[0]*ArraySize[1]*ArraySize[2];

		/*----------------------------------------------------*/
		/* Pre-calculation */
		/* Weights */
		{ int x,y,z,i;
			 WeightZ = (double *)malloc( sizeof(double)* (int)((ToZ-FromZ+1)+(ToY-FromY+1)+(ToX-FromX+1)) );
			 if (WeightZ==0) { printf("Memory allocation error\n"); return; } /* Error */
			 WeightY = &WeightZ[ToZ-FromZ+1];
			 WeightX = &WeightY[ToY-FromY+1];
			 for(x=FromX,i=0; x<=ToX; x++,i++) { WeightX[i] = cos(2.0*M_PI*(((double)x-(ToX+FromX)/2.0)/(double)(ToX-FromX+1)))*0.5+0.5; }
			 for(y=FromY,i=0; y<=ToY; y++,i++) { WeightY[i] = cos(2.0*M_PI*(((double)y-(ToY+FromY)/2.0)/(double)(ToY-FromY+1)))*0.5+0.5; }
			 for(z=FromZ,i=0; z<=ToZ; z++,i++) { WeightZ[i] = cos(2.0*M_PI*(((double)z-(ToZ+FromZ)/2.0)/(double)(ToZ-FromZ+1)))*0.5+0.5; }
    }
		aA = bB = cC = aB_r = aB_i = aC_r = aC_i = bC_r = bC_i = 0.0;
		if (img_i==NULL) { /* Real image only */
			if (DataType == 1) { /* Floats */
				int x,y,z,idxZ,idxY,idxX;
				float *float_r;
				float_r = (float *)img_r;
				idxZ = FromZ*ArraySize[0]*ArraySize[1];
				for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
					for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
						for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
							aA += ABS2((double)float_r[idxX               ],0.0);
							bB += ABS2((double)float_r[idxX+ NumPixels    ],0.0);
							cC += ABS2((double)float_r[idxX+(NumPixels<<1)],0.0);
							aB_r += float_r[idxX               ]*float_r[idxX+ NumPixels    ];
							aC_r += float_r[idxX               ]*float_r[idxX+(NumPixels<<1)];
							bC_r += float_r[idxX+ NumPixels    ]*float_r[idxX+(NumPixels<<1)];
						}
					}
				}
			} else if (DataType == 2) { /* Doubles */
				int x,y,z,idxZ,idxY,idxX;
				double *double_r;
				double_r = (double *)img_r;
				idxZ = FromZ*ArraySize[0]*ArraySize[1];
				for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
					for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
						for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
							aA += ABS2(double_r[idxX               ],0.0);
							bB += ABS2(double_r[idxX+ NumPixels    ],0.0);
							cC += ABS2(double_r[idxX+(NumPixels<<1)],0.0);
							aB_r += double_r[idxX               ]*double_r[idxX+ NumPixels    ];
							aC_r += double_r[idxX               ]*double_r[idxX+(NumPixels<<1)];
							bC_r += double_r[idxX+ NumPixels    ]*double_r[idxX+(NumPixels<<1)];
						}
					}
				}
			}
		} else {
			if (DataType == 1) { /* Floats */
				int x,y,z,idxZ,idxY,idxX;
				float *float_r, *float_i;
				float_r = (float *)img_r;
				float_i = (float *)img_i;
				idxZ = FromZ*ArraySize[0]*ArraySize[1];
				for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
					for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
						for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
							aA += ABS2((double)float_r[idxX               ],(double)float_i[idxX               ]);
							bB += ABS2((double)float_r[idxX+ NumPixels    ],(double)float_i[idxX+ NumPixels    ]);
							cC += ABS2((double)float_r[idxX+(NumPixels<<1)],(double)float_i[idxX+(NumPixels<<1)]);
							aB_r += (((double)float_r[idxX               ]*(double)float_r[idxX+ NumPixels    ]) - ((double)float_i[idxX               ]*(-(double)float_i[idxX+ NumPixels    ])));
							aC_r += (((double)float_r[idxX               ]*(double)float_r[idxX+(NumPixels<<1)]) - ((double)float_i[idxX               ]*(-(double)float_i[idxX+(NumPixels<<1)])));
							bC_r += (((double)float_r[idxX+ NumPixels    ]*(double)float_r[idxX+(NumPixels<<1)]) - ((double)float_i[idxX+ NumPixels    ]*(-(double)float_i[idxX+(NumPixels<<1)])));
							aB_i += (((double)float_r[idxX               ]*(-(double)float_i[idxX+ NumPixels    ])) + ((double)float_i[idxX               ]*(double)float_r[idxX+ NumPixels    ]));
							aC_i += (((double)float_r[idxX               ]*(-(double)float_i[idxX+(NumPixels<<1)])) + ((double)float_i[idxX               ]*(double)float_r[idxX+(NumPixels<<1)]));
							bC_i += (((double)float_r[idxX+ NumPixels    ]*(-(double)float_i[idxX+(NumPixels<<1)])) + ((double)float_i[idxX+ NumPixels    ]*(double)float_r[idxX+(NumPixels<<1)]));
						}
					}
				}
			} else if (DataType == 2) { /* Doubles */
				int x,y,z,idxZ,idxY,idxX;
				double *double_r, *double_i;
				double_r = (double *)img_r;
				double_i = (double *)img_i;
				idxZ = FromZ*ArraySize[0]*ArraySize[1];
				for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
					for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
						for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
							aA += ABS2(double_r[idxX               ],double_i[idxX               ]);
							bB += ABS2(double_r[idxX+ NumPixels    ],double_i[idxX+ NumPixels    ]);
							cC += ABS2(double_r[idxX+(NumPixels<<1)],double_i[idxX+(NumPixels<<1)]);
							aB_r += (double_r[idxX               ]*double_r[idxX+ NumPixels    ] - double_i[idxX               ]*(-double_i[idxX+ NumPixels    ]));
							aC_r += (double_r[idxX               ]*double_r[idxX+(NumPixels<<1)] - double_i[idxX               ]*(-double_i[idxX+(NumPixels<<1)]));
							bC_r += (double_r[idxX+ NumPixels    ]*double_r[idxX+(NumPixels<<1)] - double_i[idxX+ NumPixels    ]*(-double_i[idxX+(NumPixels<<1)]));
							aB_i += (double_r[idxX               ]*(-double_i[idxX+ NumPixels    ]) + double_i[idxX               ]*double_r[idxX+ NumPixels    ]);
							aC_i += (double_r[idxX               ]*(-double_i[idxX+(NumPixels<<1)]) + double_i[idxX               ]*double_r[idxX+(NumPixels<<1)]);
							bC_i += (double_r[idxX+ NumPixels    ]*(-double_i[idxX+(NumPixels<<1)]) + double_i[idxX+ NumPixels    ]*double_r[idxX+(NumPixels<<1)]);
						}
					}
				}
			}
			ClientData[0] = aA;    ClientData[1] = bB; ClientData[2] = cC;
			ClientData[3] = aB_r;  ClientData[4] = aB_i;
			ClientData[5] = aC_r;  ClientData[6] = aC_i;
			ClientData[7] = bC_r;  ClientData[8] = bC_i;
		}
		/* End of pre-calculation */
		/*----------------------------------------------------*/
		if (0) { /* Check noise level is see if there is useful info or not */
			double NumPixelsInVOI;
			NumPixelsInVOI = (FromX-ToX+1)*(FromY-ToY+1)*(FromZ-ToZ+1);
			if (aA < (*NoiseVar)*NumPixelsInVOI) *NoiseVar = aA/NumPixelsInVOI;
			if (bB < (*NoiseVar)*NumPixelsInVOI) *NoiseVar = bB/NumPixelsInVOI;
			if (cC < (*NoiseVar)*NumPixelsInVOI) *NoiseVar = cC/NumPixelsInVOI;
			if ((aA + bB + cC) < NumPixelsInVOI*3*(*NoiseVar)*NoiseMargin) {
				Optimize = 0;
			} else {
				Optimize = 1;
			}
		}
		/*----------------------------------------------------*/
		/* Find best angle (start) */
		if (Optimize) {
			int i;
			double FuncToMinimize(double *Param, double *ClientData);
			double CurrentResult, BestResult;

			if (DegreeRange==180.0) { /* Global search */
		#define NUMANGLEDIV     45
				double BestDegree;
				Param[1] = t2correct;
				BestResult = -1;
				for(i=-NUMANGLEDIV; i<NUMANGLEDIV; i++) {
					double CurrentDegree;
					CurrentDegree = DegreeOffset + (DegreeRange/(double)NUMANGLEDIV)*((double)i);
					Param[0] = CurrentDegree;
					CurrentResult = FuncToMinimize(Param, ClientData);
					if (BestResult==-1 || CurrentResult<BestResult) {
						BestResult = CurrentResult;
						BestDegree = CurrentDegree;
					}
				}
				DegreeOffset = BestDegree;
		#undef NUMANGLEDIV

			}

			if (1) { /* Local search (Not as robust) */
				SimplexDataType Simplex;
				double ParamStepSize[] = {5, 0.05}; /* degrees, t2correct */
				double MinParamStepSize[] = {0.01, 0.0001}; /* degrees, t2correct */
				double MinParam[2], MaxParam[2];
				double ParamMaxChange[] = {40, 5}; /* degrees, t2correct */
				double StartParam[2];
				int MaxNumEvals = 1000;  /* maximum number of iterations. <=0 if unlimited */
				StartParam[0] = DegreeOffset;
				StartParam[1] = t2correct;
				MinParam[0] = DegreeOffset - DegreeRange;
				MaxParam[0] = DegreeOffset + DegreeRange;
				MinParam[1] = 1.0;
				MaxParam[1] = 50.0;
				SimplexSetup(ClientData, &Simplex, 2, StartParam, ParamStepSize, MinParam, MaxParam);
				SimplexMinimize(ClientData, &Simplex, 2, MinParamStepSize, MaxNumEvals);
				if (abs(Simplex.Vertex[0][0]-DegreeOffset)>ParamMaxChange[0] || abs(t2correct-Simplex.Vertex[0][1])>ParamMaxChange[1]) {
					i=1; /* optimization failed */
				} else {
					DegreeOffset = Simplex.Vertex[0][0];
					t2correct    = Simplex.Vertex[0][1];
					i=0;
				}
   		} else i=1;

			if (i) { /* Find best T2* correction & best angle (start) */
				int t2step;
				double BestParam;
				int NumT2CorrectSteps = 20, NumAngleSteps = 9;
				double T2CorrectStep = 0.01;       
				BestParam = t2correct;
				Param[0] = DegreeOffset;
				BestResult=-1;
				for(t2step=0; t2step<=NumT2CorrectSteps; t2step++) {
					Param[1] = t2correct * pow((1+T2CorrectStep),t2step);
					if (Param[1]<1.0) continue;
					CurrentResult = FuncToMinimize(Param, ClientData);
					if (BestResult==-1 || CurrentResult < BestResult) {
						BestResult = CurrentResult;
						BestParam = Param[1];
					}
				}
				for(t2step=1; t2step<=NumT2CorrectSteps; t2step++) {
					Param[1] = t2correct * pow((1+T2CorrectStep),-t2step);
					if (Param[1]<1.0) continue;
					CurrentResult = FuncToMinimize(Param, ClientData);
					if (CurrentResult < BestResult) {
						BestResult = CurrentResult;
						BestParam = Param[1];
					}
				}
				t2correct = BestParam;

				BestParam = DegreeOffset;
				BestResult=-1;
				Param[1] = t2correct;
				for(i=0; i<=NumAngleSteps; i++) {
					Param[0] = DegreeOffset + DegreeRange*((double)i)/((double)NumAngleSteps);
					CurrentResult = FuncToMinimize(Param, ClientData);
					if (BestResult==-1 ||CurrentResult < BestResult) {
						BestResult = CurrentResult;
						BestParam = Param[0];
					}
				}
				for(i=1; i<=NumAngleSteps; i++) {
					Param[0] = DegreeOffset - DegreeRange*((double)i)/((double)NumAngleSteps);
					CurrentResult = FuncToMinimize(Param, ClientData);
					if (CurrentResult < BestResult) {
						BestResult = CurrentResult;
						BestParam = Param[0];
					}
				}
				DegreeOffset = BestParam;
			} /* Find best T2* correction & best angle (end) */

		}
		/*----------------------------------------------------*/

	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Subdivide */
	{
		int tmpX, tmpY, tmpZ;
		int LowerFromX, LowerToX,
		    LowerFromY, LowerToY,
		    LowerFromZ, LowerToZ;
		int UpperFromX, UpperToX,
		    UpperFromY, UpperToY,
		    UpperFromZ, UpperToZ;
		double overlapfraction = 0.33333333333333333;
		tmpX = ToX-FromX+1;
		tmpY = ToY-FromY+1;
		tmpZ = ToZ-FromZ+1;
		LowerFromX = FromX;
		LowerFromY = FromY;
		LowerFromZ = FromZ;
		UpperToX = ToX;
		UpperToY = ToY;
		UpperToZ = ToZ;
		{
			int OverlapPixels, lowersize, uppersize;
			OverlapPixels = ceil(tmpX*overlapfraction);
			lowersize = ceil((tmpX-OverlapPixels)/2); uppersize = (tmpX-OverlapPixels)-lowersize;
			lowersize += OverlapPixels; uppersize += OverlapPixels;
			LowerToX = LowerFromX + lowersize-1;  UpperFromX = LowerToX - OverlapPixels+1;

			OverlapPixels = ceil(tmpY*overlapfraction);
			lowersize = ceil((tmpY-OverlapPixels)/2); uppersize = (tmpY-OverlapPixels)-lowersize;
			lowersize += OverlapPixels; uppersize += OverlapPixels;
			LowerToY = LowerFromY + lowersize-1;  UpperFromY = LowerToY - OverlapPixels+1;

			OverlapPixels = ceil(tmpZ*overlapfraction);
			lowersize = ceil((tmpZ-OverlapPixels)/2); uppersize = (tmpZ-OverlapPixels)-lowersize;
			lowersize += OverlapPixels; uppersize += OverlapPixels;
			LowerToZ = LowerFromZ + lowersize-1;  UpperFromZ = LowerToZ - OverlapPixels+1;
		}
/*
    printf("DEPTH %d subdiv -> %d %d, %d %d\n",depth,UpperFromX-FromX+1,LowerToX-FromX+1,UpperFromY-FromY+1,LowerToY-FromY+1);
 */
		if (0) printf("DEPTH %d subdiv -> x (%d,%d,%d,%d)  y (%d,%d,%d,%d)\n",depth, LowerFromX,LowerToX,UpperFromX,UpperToX, LowerFromY,LowerToY,UpperFromY,UpperToY);

		if (LowerFromX==UpperFromX || LowerToX==UpperToX || tmpX<MinSizeX) {
			LowerFromX = FromX;
			LowerToX   = ToX;
			tmpX = 0; /* No need to divide */
		} else tmpX = 1;
		if (LowerFromY==UpperFromY || LowerToY==UpperToY || tmpY<MinSizeY) {
			LowerFromY = FromY;
			LowerToY   = ToY;
			tmpY = 0; /* No need to divide */
		} else tmpY = 1;
		if (LowerFromZ==UpperFromZ || LowerToZ==UpperToZ || tmpZ<MinSizeZ) {
			LowerFromZ = FromZ;
			LowerToZ   = ToZ;
			tmpZ = 0; /* No need to divide */
		} else tmpZ = 1;
		/*----------------------------------------------------------------*/

		/*----------------------------------------------------------------*/
		/* Next level */
		depth++;
		DegreeRange *= DegreeRangeContraction;
		if (depth>=DegreeRangeMax_depth) { /* 2009.04.14 Avoid restraining range too soon */
			if (DegreeRange>DegreeRangeMax) DegreeRange=DegreeRangeMax;  /* less than 1/3 of 180deg */
		}
		if (DegreeRange<DegreeRangeMin) DegreeRange=DegreeRangeMin;

		/*----------------------------------------------------------------*/

		/*----------------------------------------------------------------*/
		if ((Optimize==0) || (tmpX==0 && tmpY==0 && tmpZ==0) || (depth>MaxDepth && (MaxDepth>=0))) { /* if MaxDepth<0, it means that MaxDepth is not given. */
			/* Small enough, finished recursion */
			int x,y,z,idxZ,idxY,idxX;
			double PixelPhase_r, PixelPhase_i;
			double WeightYZ, WeightXYZ;
/*			EXP_I(PixelPhase_r, PixelPhase_i, (DegreeOffset-120)/180.0*M_PI); */
			EXP_I(PixelPhase_r, PixelPhase_i, DegreeOffset/180.0*M_PI);
			for(z=FromZ, idxZ=z*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
				for(y=FromY, idxY=idxZ+y*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
					WeightYZ = WeightY[y-FromY] * WeightZ[z-FromZ];
					for(x=FromX, idxX=idxY+x; x<=ToX; x++, idxX++) {
						WeightXYZ = WeightYZ * WeightX[x-FromX];
						phasemap_r[idxX] += PixelPhase_r*WeightXYZ;
						phasemap_i[idxX] += PixelPhase_i*WeightXYZ;
						t2correctmap_r[idxX] += t2correct*WeightXYZ;
						t2correctmap_i[idxX] += WeightXYZ;
					}
				}
			}
			free(WeightZ);
		} else {
			free(WeightZ);
			/* Subdivide into small patches */
			IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
			                  NumDims, ArraySize, depth, MaxDepth,
			                  MinSizeX,MinSizeY,MinSizeZ,
			                  LowerFromX, LowerToX,
			                  LowerFromY, LowerToY,
			                  LowerFromZ, LowerToZ,
			                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			if (tmpX) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  UpperFromX, UpperToX,
				                  LowerFromY, LowerToY,
				                  LowerFromZ, LowerToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpY) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  LowerFromX, LowerToX,
				                  UpperFromY, UpperToY,
				                  LowerFromZ, LowerToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpZ) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  LowerFromX, LowerToX,
				                  LowerFromY, LowerToY,
				                  UpperFromZ, UpperToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpX && tmpY) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i,t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  UpperFromX, UpperToX,
				                  UpperFromY, UpperToY,
				                  LowerFromZ, LowerToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpX && tmpZ) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i, t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  UpperFromX, UpperToX,
				                  LowerFromY, LowerToY,
				                  UpperFromZ, UpperToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpY && tmpZ) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i, t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  LowerFromX, LowerToX,
				                  UpperFromY, UpperToY,
				                  UpperFromZ, UpperToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
			if (tmpX && tmpY && tmpZ) {
				IDEAL_engine_core(img_r,img_i,DataType,phasemap_r,phasemap_i, t2correctmap_r,t2correctmap_i,
				                  NumDims, ArraySize, depth, MaxDepth,
				                  MinSizeX,MinSizeY,MinSizeZ,
				                  UpperFromX, UpperToX,
				                  UpperFromY, UpperToY,
				                  UpperFromZ, UpperToZ,
				                  DegreeRange, DegreeOffset, t2correct, NoiseVar);
			}
		}
	}
	/*----------------------------------------------------------------*/

	/*----------------------------------------------------------------*/
	/* Finished. Normalize by phasemap = exp(i*angle(phasemap)); */
	if (depth==1) {
		int x,y,z,idxZ,idxY,idxX;
		double Magnitude;

		idxZ = FromZ*ArraySize[0]*ArraySize[1];
		for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
			for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
				for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
					Magnitude = ABS2(phasemap_r[idxX],phasemap_i[idxX]);
					if (Magnitude>0.0) {
						Magnitude = 1.0/sqrt(Magnitude);
						phasemap_r[idxX] *= Magnitude;
						phasemap_i[idxX] *= Magnitude;
					}
				}
			}
		}

		idxZ = FromZ*ArraySize[0]*ArraySize[1];
		for(z=FromZ, idxZ=FromZ*ArraySize[0]*ArraySize[1]; z<=ToZ; z++,idxZ+=ArraySize[0]*ArraySize[1]) {
			for(y=FromY, idxY=idxZ+FromY*ArraySize[0]; y<=ToY; y++, idxY+=ArraySize[0]) {
				for(x=FromX, idxX=idxY+FromX; x<=ToX; x++, idxX++) {
					if (t2correctmap_i[idxX]!=0.0) {
						t2correctmap_r[idxX]/=t2correctmap_i[idxX];
						t2correctmap_i[idxX] = 0;
					}
				}
			}
		}

	}
	/*----------------------------------------------------------------*/
}

void SimplexSetup(double          *ClientData,             /*::: Input :::*/
                  SimplexDataType *Simplex, /*::: Input & Output :::*/
                  int NumDim,                          /*::: Input WARNING: NumDim <= 3 :::*/
                  double          *StartParam, /*::: Input ::: */
                  double          *ParamStepSize, /*::: Input :::*/
                  double          *MinParam, /*::: Input :::*/
                  double          *MaxParam) /*::: Input :::*/
{
	/*::: Initialize values ::::::::::::::::::::::::::::::::::
	   Vertex[(0 to NumDim-1)][(0 to NumDim)] = array of vectices
	    For NumDim=2,3: Platonic shapes
	    Otherwise:
	      (0,0,0,0,....)
	      (1,0,0,0,....)     NOTE: 0 = 0
	      (0,1,0,0,....)           1 = StepSize
	        :
	      (0,0,0,0,..,1)

	   FuncValue[0 to NumDim] = array of function values evaluated at
	         the above vertices
	   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	double FuncToMinimize(double *Param, double *ClientData);

	switch(NumDim) {
	case 2:
		Simplex->Vertex[0][0]=0.0; Simplex->Vertex[0][1]=0.0;
		Simplex->Vertex[1][0]=1.0; Simplex->Vertex[1][1]=0.0;
		Simplex->Vertex[2][0]=0.5; Simplex->Vertex[2][1]=0.8660254038;
		break;
	case 3:
		Simplex->Vertex[0][0]=0.0;          Simplex->Vertex[0][1]=0.0;          Simplex->Vertex[0][2]=0.0;
		Simplex->Vertex[1][0]=0.0;          Simplex->Vertex[1][1]=0.7071067812; Simplex->Vertex[1][2]=0.7071067812;
		Simplex->Vertex[2][0]=0.7071067812; Simplex->Vertex[2][1]=0.7071067812; Simplex->Vertex[2][2]=0.0;
		Simplex->Vertex[3][0]=0.7071067812; Simplex->Vertex[3][1]=0.0;          Simplex->Vertex[3][2]=0.7071067812;
		break;
	default:
	{ register int i,j;
	  for(j=0; j<=NumDim; j++)
		  for(i=0; i< NumDim; i++) Simplex->Vertex[j][i]=0;
	  for(i=0; i< NumDim; i++) Simplex->Vertex[i+1][i]=1; }
	}

	/*::: Set up vertices as ParamStepSize * Directions ::::::::::::*/
	{ register int i, j;
	  for(i=0; i< NumDim; i++) {
		  if (ParamStepSize[i]==0) continue;
		  for(j=0; j<=NumDim; j++) Simplex->Vertex[j][i]*=ParamStepSize[i];
	  } }

	/*::: Vertices relative to center of mass ::::::::::::::::::::::::::::*/
	if (SIMPLEX_UNBIASED_DIRECTION) { register int i, j;
		                          register double tmp;
		                          for(i=0; i<NumDim; i++) {
						  for(j=0, tmp=0.0; j<=NumDim; j++) tmp+=Simplex->Vertex[j][i];
						  tmp /= (double)(NumDim+1);
						  for(j=0, tmp=0.0; j<=NumDim; j++) Simplex->Vertex[j][i] = (Simplex->Vertex[j][i]-tmp) *2.0;
					  } }

	/*::: Add starting param to vertices ::::::::::::*/
	if (StartParam) {
		register int i, j;
		for(i=0; i< NumDim; i++) for(j=0; j<=NumDim; j++) Simplex->Vertex[j][i]+=StartParam[i];
	}

	/*::: Add min max for params ::::::::::::*/
	{
		register int i;
		for(i=0; i< NumDim; i++) {
			Simplex->MinParam[i]=MinParam[i];
			Simplex->MaxParam[i]=MaxParam[i];
		}
	}

	/*::::::::::::::::::::::::::::::::::::::::::::::::::::::::
	   FuncValue[0 to NumDim] = array of function values evaluated at
	    the above vertices
	   ::::::::::::::::::::::::::::::::::::::::::::::::::::::::*/
	{ register int j;
	  for(j=0; j<=NumDim; j++)
		  Simplex->FuncValue[j] = FuncToMinimize(Simplex->Vertex[j], ClientData); }
	Simplex->NumEvals = NumDim+1;      /*::: Number of function evaluations :::*/

	/*::: Evaluate function at center of mass :::::::::::::::::::::::::::::*/
	if (SIMPLEX_UNBIASED_DIRECTION) { register int i, j;
		                          register double tmpFuncValue;
		                          double tmpVertex[] = {0,0,0,0};
		                          tmpFuncValue = FuncToMinimize(tmpVertex, ClientData);
		                          (Simplex->NumEvals)++;

		                          /*::: Replace worst vertex with this one :::*/
		                          tmpFuncValue = Simplex->FuncValue[0]; j=0;
		                          for(i=1; i<=NumDim; i++)
						  if (tmpFuncValue < Simplex->FuncValue[i]) { tmpFuncValue = Simplex->FuncValue[i]; j=i; }
		                          Simplex->FuncValue[j] = tmpFuncValue;
		                          for(i=0; i<NumDim; i++) Simplex->Vertex[j][i] = 0.0; }

	{ register int i,j; register double temp;
	  for(i=0; i<NumDim; i++) {     /*::: Calculate centroid :::*/
		  for(j=0, temp=0.0; j<=NumDim; j++) temp += Simplex->Vertex[j][i];
		  Simplex->Centroid[i] = temp;
	  } }
}
/***********************************************************/
/***********************************************************/
/***********************************************************/
int SimplexMinimize(double          *ClientData,                       /*::: Input :::*/
                    SimplexDataType *Simplex, /*::: Input & Output :::*/
                    int NumDim,                        /*::: Input :::*/
                    double          *MinParamStepSize, /*::: Input :::*/
                    int MaxNumEvals)       /*::: Input :::*/
{
	double FuncToMinimize(double *Param, double *ClientData);
	double SimplexMinimizeExtrapolateVertices( double *ClientData,
	                                           SimplexDataType *Simplex,
	                                           int NumDim,
	                                           int idxWorst, double fac);
	int idxWorst, idxBest, idx2ndWorst;
	double tempValue;
	int ReturnValue = 0;

	do {
		/*::: Find best, worst and 2nd worst vertices ::::::::*/
		idxBest=0;
		if (Simplex->FuncValue[0]>Simplex->FuncValue[1]) { idxWorst = 0; idx2ndWorst = 1; }
		else { idxWorst = 1; idx2ndWorst = 0; }
		{ register int j;
		  for(j=0; j<=NumDim; j++) {
			  if      (Simplex->FuncValue[j] <= Simplex->FuncValue[idxBest    ]) idxBest=j;
			  if      (Simplex->FuncValue[j] >  Simplex->FuncValue[idxWorst   ]) { idx2ndWorst=idxWorst; idxWorst=j; }
			  else if (Simplex->FuncValue[j] >  Simplex->FuncValue[idx2ndWorst] && j!=idxWorst) idx2ndWorst=j;
		  } }

		if (Simplex->NumEvals >= MaxNumEvals && MaxNumEvals>0) { ReturnValue=1; break; } /*::: Enough -> Quit :::*/

		/*::: Reflection :::::::::::::::::::::::::::::::*/
		tempValue = SimplexMinimizeExtrapolateVertices(ClientData, Simplex, NumDim, idxWorst, -1.0);
		Simplex->NumEvals += 1;

		if        (tempValue <= Simplex->FuncValue[idxBest]) {
			/*::: Reflection & Expansion :::::::::::::::*/
			tempValue = SimplexMinimizeExtrapolateVertices(ClientData, Simplex, NumDim, idxWorst, 2.0);
			Simplex->NumEvals += 1;
		} else if (tempValue >= Simplex->FuncValue[idx2ndWorst]) { /* Either contraction or multiple contraction */
			/*::: See if simplex bounding box is small enough ::::*/
			{ register int i,j;
			  register double maxPos, minPos; /* Calculate bounding box of simplex */
			  for(i=0; i<NumDim; i++) { /* Each coord */
				  maxPos = minPos = Simplex->Vertex[0][i];
				  for(j=1; j<=NumDim; j++) { /* For this coord, find min and max of NumDim+1 vectors */
					  if      (maxPos < Simplex->Vertex[j][i]) maxPos = Simplex->Vertex[j][i];
					  else if (minPos > Simplex->Vertex[j][i]) minPos = Simplex->Vertex[j][i];
				  }
				  if (maxPos-minPos > MinParamStepSize[i]) break;  /*::: Not converge completely :::*/
			  }
			  if (i>=NumDim) break;  /*::: Converged -> Quit :::*/
			}

			/*::: Contraction ::::::::::::::::::::::::::*/
			{ double oldValue;
			  oldValue = Simplex->FuncValue[idxWorst];
			  tempValue = SimplexMinimizeExtrapolateVertices(ClientData, Simplex, NumDim, idxWorst, 0.5);
			  Simplex->NumEvals += 1;
			  if (tempValue >= oldValue) {
				  /*::: Multiple contraction :::::::::::::::::::::*/
				  { register int i,j;
		    for(i=0; i<=NumDim; i++) {
			    if (i==idxBest) continue;
			    for(j=0; j<NumDim; j++)
				    Simplex->Vertex[i][j] =
				            Simplex->Centroid[j] = 0.5*(Simplex->Vertex[i][j] + Simplex->Vertex[idxBest][j]);
			    Simplex->FuncValue[i] = FuncToMinimize(Simplex->Centroid, ClientData);
		    } }
				  Simplex->NumEvals += NumDim;

				  { register int i,j;  register double temp=0.0;
		    for(j=0; j<NumDim; j++) { /*::: Calculate centroid :::*/
			    for(i=0; i<=NumDim; i++) temp += Simplex->Vertex[i][j];
			    Simplex->Centroid[j] = temp;
		    } }
			  } }
		}
	} while (1);

	if (idxBest!=0) { /*::: Make Vertex[0] the best vertex :::*/
		register int i; register double temp;
		for(i=0; i<NumDim; i++) {
			temp = Simplex->Vertex[0][i];
			Simplex->Vertex[0][i] = Simplex->Vertex[idxBest][i];
			Simplex->Vertex[idxBest][i] = temp;
		}
		temp = Simplex->FuncValue[0];
		Simplex->FuncValue[0] = Simplex->FuncValue[idxBest];
		Simplex->FuncValue[idxBest] = temp;
	}

	return ReturnValue;
}

double SimplexMinimizeExtrapolateVertices(double *ClientData,
                                          SimplexDataType *Simplex,
                                          int NumDim,
                                          int idxWorst, double fac)
{
	double FuncToMinimize(double *Param, double *ClientData);
	double NewVertex[3]; /* Assume NumDim <= 3 dimensions */
	double NewFuncValue;
	register int j;

	{ double fac1, fac2;

	  fac1=(1.0-fac)/(double)NumDim;
	  fac2=fac1-fac;
	  for(j=0; j<NumDim; j++) {
		  NewVertex[j] = Simplex->Centroid[j]        *fac1
		                 - Simplex->Vertex[idxWorst][j]*fac2;
		  if      (NewVertex[j]<Simplex->MinParam[j]) NewVertex[j] = Simplex->MinParam[j];  /* New position */
		  else if (NewVertex[j]>Simplex->MaxParam[j]) NewVertex[j] = Simplex->MaxParam[j];  /* New position */
	  } }

	NewFuncValue = FuncToMinimize(NewVertex, ClientData);

	if (Simplex->FuncValue[idxWorst] > NewFuncValue) {
		Simplex->FuncValue[idxWorst] = NewFuncValue;
		for(j=0; j<NumDim; j++) {
			Simplex->Centroid[j] += NewVertex[j] - Simplex->Vertex[idxWorst][j];
			Simplex->Vertex[idxWorst][j]  = NewVertex[j];
		}
	}
	return NewFuncValue;
}

/*
 * Signals at 3 echo times are a, b, c
 * Want to find | a*p' + b*t2 + c*t2^2*phase |
 *            = | ap' + b(t2) + cp(t2^2) |  where p = t2cor*phase
 *            = aa'(pp') + bb'(t2^2) + cc'(pp')(t2^4) + 2*re(ap'b'(t2)) + 2*re(ac'p'p'(t2^2)) + 2*re(bc'p'(t2^3))
 *            = aa' + bb'(t2^2) + cc'(t2^4) + 2*re(ab'(p*t2)') + 2*re(ac'((p*t2)^2)') + 2*re(bc'(p*t2)'(t2^2))
 */
double FuncToMinimize(double *Param, double *ClientData) {
	double aA, aB_r, aB_i, aC_r, aC_i, bB, bC_r, bC_i, cC;

	double DegreeOffset, t2correct, t2correct2;
	double returnvalue;
	DegreeOffset = Param[0];
	t2correct = Param[1];
	aA = ClientData[0];
	bB = ClientData[1];
	cC = ClientData[2];
	aB_r = ClientData[3];  aB_i = ClientData[4];
	aC_r = ClientData[5];  aC_i = ClientData[6];
	bC_r = ClientData[7];  bC_i = ClientData[8];
	t2correct2 = t2correct*t2correct;

  {
    double p_t2_r, p_t2_i, p_t2_sq_r, p_t2_sq_i;
    EXP_I(p_t2_r, p_t2_i, (DegreeOffset-120.0)/180.0*M_PI); p_t2_r *= t2correct; p_t2_i *= t2correct;
	  COMPLEXMULT(p_t2_sq_r, p_t2_sq_i, p_t2_r, p_t2_i, p_t2_r, p_t2_i);  
	  returnvalue = aA + t2correct2*(bB + cC*t2correct2) 
                + 2 * (              COMPLEXMULT_r(aB_r,aB_i, p_t2_r,    -p_t2_i   )    /* 2*re(ab'(p*t2)')       */
	                    +              COMPLEXMULT_r(aC_r,aC_i, p_t2_sq_r, -p_t2_sq_i)    /* 2*re(ac'((p*t2)^2)')   */
	                    +   t2correct2*COMPLEXMULT_r(bC_r,bC_i, p_t2_r,    -p_t2_i   ) ); /* 2*re(bc'(p*t2)'(t2^2)) */
  }

	return returnvalue;
}

double EstimateNoise(void *img_r, void *img_i, int DataType, int NumDims, int *ArraySize)
{
	double NoiseVar = -1;
	int NumPixels;
	NumPixels = ArraySize[0]*ArraySize[1]*ArraySize[2];

	if (img_i==NULL) { /* Real image only */
		if (DataType == 1) { /* Floats */
			int x,y,z,i;
			float *float_r;
			float_r = (float *)img_r;
			for(z=0,i=0; z<ArraySize[2]; z++) {
				for(y=0; y<ArraySize[1]; y++) {
					double TempVar=0.0;
					for(x=0; x<ArraySize[0]; x++,i++) TempVar += ABS2((double)float_r[i],0.0) + ABS2((double)float_r[i+NumPixels],0.0) + ABS2((double)float_r[i+(NumPixels<<1)],0.0);
					if (NoiseVar<0 || TempVar<NoiseVar) NoiseVar = TempVar;
				}
			}
		} else if (DataType == 2) { /* Doubles */
			int x,y,z,i;
			double *double_r;
			double_r = (double *)img_r;
			for(z=0,i=0; z<ArraySize[2]; z++) {
				for(y=0; y<ArraySize[1]; y++) {
					double TempVar=0.0;
					for(x=0; x<ArraySize[0]; x++,i++) TempVar += ABS2(double_r[i],0.0) + ABS2(double_r[i+NumPixels],0.0) + ABS2(double_r[i+(NumPixels<<1)],0);
					if (NoiseVar<0 || TempVar<NoiseVar) NoiseVar = TempVar;
				}
			}
		}
	} else {
		if (DataType == 1) { /* Floats */
			int x,y,z,i;
			float *float_r, *float_i;
			float_r = (float *)img_r;
			float_i = (float *)img_i;
			for(z=0,i=0; z<ArraySize[2]; z++) {
				for(y=0; y<ArraySize[1]; y++) {
					double TempVar=0.0;
					for(x=0; x<ArraySize[0]; x++,i++) TempVar += ABS2((double)float_r[i],(double)float_i[i]) + ABS2((double)float_r[i+NumPixels],(double)float_i[i+NumPixels]) + ABS2((double)float_r[i+(NumPixels<<1)],(double)float_i[i+(NumPixels<<1)]);
					if (NoiseVar<0 || TempVar<NoiseVar) NoiseVar = TempVar;
				}
			}
		} else if (DataType == 2) { /* Doubles */
			int x,y,z,i;
			double *double_r, *double_i;
			double_r = (double *)img_r;
			double_i = (double *)img_i;
			for(z=0,i=0; z<ArraySize[2]; z++) {
				for(y=0; y<ArraySize[1]; y++) {
					double TempVar=0.0;
					for(x=0; x<ArraySize[0]; x++,i++) TempVar += ABS2(double_r[i],double_i[i]) + ABS2(double_r[i+NumPixels],double_i[i+NumPixels]) + ABS2(double_r[i+(NumPixels<<1)],double_i[i+(NumPixels<<1)]);
					if (NoiseVar<0 || TempVar<NoiseVar) NoiseVar = TempVar;
				}
			}
		}
	}

	NoiseVar = NoiseVar/(double)(ArraySize[0]*3); /* Per pixel */
	return NoiseVar;
}

#undef COMPLEXCONJ
#undef COMPLEXCONJ_ACCUM
#undef COMPLEXADD
#undef COMPLEXADD_ACCUM
#undef COMPLEXMULT
#undef COMPLEXMULT_ACCUM
#undef COMPLEXREALMULT
#undef COMPLEXREALMULT_ACCUM
#undef COMPLEXCONJ_INLINE
#undef COMPLEXADD_INLINE
#undef COMPLEXMULT_INLINE
#undef EXP_I
#undef ABS2

char ToLowercase(char c) {
	return ( (c)>='a' && (c)<='z' ? ((c)-'a'+'A') : (c) );
}
int CompareString_IgnoreCase(char *s1, char *s2)
{
	int i;
	for(i=0;; i++) {
		if (ToLowercase(s1[i])==ToLowercase(s2[i])) {
			if (s1[i]=='\0') return 0;
		} else {
			return (ToLowercase(s1[i])<ToLowercase(s2[i]) ? -1 : +1);
		}
	}
}



