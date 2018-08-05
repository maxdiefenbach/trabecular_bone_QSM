/*

 FWT2_PO.C	.MEX file corresponding to fwt2_po.m
		Periodized Wavelet Transform

  The calling syntax is:

			wc = fwt2_po(image,L,qmf)


  David Donoho
  Copyright (c) 1993  David Donoho
  All Rights Reserved

*/

#include <math.h>
#include "mex.h"
#include "wavelab.h"

void dpwt2(double *sig,int nr,int nc,int ell,int J,
   double *hpf,double *lpf,int lenfil,double *wc,double *temp);


#define DOUBLE		double
#define INT			int

/* Input Arguments */

#define	Sig_IN	prhs[0]
#define	Lx_IN	prhs[1]
#define	Ly_IN	prhs[2]
#define LPF_IN	prhs[3]



/* Output Arguments */

#define	WC_OUT	plhs[0]

INT nlhs, nrhs;
mxArray *plhs[], *prhs[];

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	DOUBLE	*hpf,*lpf;
	DOUBLE	*sig,*wcp,*tmp;
	unsigned int	m,n;
	int nr,nc,nn,lenfil,ellx,elly;
	mxArray *temp, *hpfmat;



	/* Check for proper number of arguments */

	if (nrhs != 4) {
		mexErrMsgTxt("FWT2_PO_TT requires 4 input arguments.");
	}
	else if (nlhs != 1) {
		mexErrMsgTxt("FWT2_PO_TT requires one output argument.");
	}


	/* Check the dimensions of signal */

	m  = mxGetM(Sig_IN);
	n = mxGetN(Sig_IN);
 
	nr = (int) m;
	nc = (int) n;


	/* Create a matrix for the return argument */

	WC_OUT = mxCreateDoubleMatrix(nr, nc, mxREAL);
	if(nr<nc)
		temp   = mxCreateDoubleMatrix(nc, 3, mxREAL);
	else 
		temp   = mxCreateDoubleMatrix(nr, 3, mxREAL);	/*temp[0]: takes one row for processing*/
													/*temp[1]: takes one half of a column (lo) for processing*/
													/*temp[2]: takes the other half of the column (hi)*/

	
	/* Assign pointers to the various parameters */

	wcp = mxGetPr(WC_OUT);
	tmp = mxGetPr(temp);

	sig = mxGetPr(Sig_IN);
    ellx = floor ((mxGetPr(Lx_IN))[0] + .5);
	elly = floor ((mxGetPr(Ly_IN))[0] + .5);
	lpf = mxGetPr(LPF_IN);
    lenfil = (int) (mxGetM(LPF_IN) * mxGetN(LPF_IN));
	hpfmat = mxCreateDoubleMatrix((unsigned int) lenfil,  1, mxREAL);
	hpf    = mxGetPr(hpfmat);
	mirrorfilt(lpf,hpf,lenfil);

	
	/*feasibility check*/
	if(nc%(int)pow(2,ellx)!=0){
		mexPrintf("Lx = %d:\n", ellx);
		mexErrMsgTxt("scale is not feasible");
	}
	if(nr%(int)pow(2,elly)!=0){
		mexPrintf("Ly = %d:\n", elly);
		mexErrMsgTxt("scale is not feasible");
	}

	/* Do the actual computations in a subroutine */

	dpwt2(sig,nr,nc,ellx,elly,hpf,lpf,lenfil,wcp,tmp);
	mxDestroyArray(temp);
	mxDestroyArray(hpfmat);
}


void dpwt2(sig,nr,nc,ellx,elly,hpf,lpf,lenfil,wc,temp)
DOUBLE sig[],hpf[],lpf[],wc[],temp[];
int  nr,nc,ellx,elly,lenfil;
{
	DOUBLE *wcplo,*wcphi,*templo,*temphi;
	int k,j,njc,njr,ell_min;
	copydouble(sig,wc,nr*nc);
	if(nr<nc) {
		templo = &temp[nc];
		temphi = &temp[2*nc];
	}
	else {
		templo = &temp[nr];
		temphi = &temp[2*nr];
	}
	
	if(ellx<elly)
		ell_min = ellx;
	else
		ell_min = elly;
		
	njr = nr; 
	njc = nc;
	for( j=0; j < ell_min; j++) {
		for( k=0; k < njc; k++) {	/*transform columns*/
		   wcplo = &wc[k*nr];
		   wcphi = &wc[k*nr + njr/2];
		   copydouble(wcplo,temp,njr); /*src: wcplo, dst: temp, copylength: njr*/
		   downlo(temp, njr, lpf,lenfil,wcplo);
		   downhi(temp, njr, hpf,lenfil,wcphi);
		}
		for( k=0; k < njr; k++) {	/*transform rows*/
		   unpackdouble(wc,njc,nr,k,temp);	/*copy one row into temp, wc: already transformed in column direction*/
		   downlo(temp, njc, lpf,lenfil,templo);
		   downhi(temp, njc, hpf,lenfil,temphi);
		   packdouble(templo,njc/2,nr,k,wc);
		   packdouble(temphi,njc/2,nr,k,&wc[njc/2*nr]);
		}
		njc = njc/2;
		njr = njr/2;
	}
	
	/*perform 1-d transform in 1 direction*/
	
	if(ell_min == ellx)	{	/*continue in y-direction*/
		for(j=0; j<elly-ell_min; j++) {
			for( k=0; k < njc; k++) {	/*transform columns*/
				wcplo = &wc[k*nr];
				wcphi = &wc[k*nr + njr/2];
				copydouble(wcplo,temp,njr); /*src: wcplo, dst: temp, copylength: njr*/
				downlo(temp, njr, lpf,lenfil,wcplo);
				downhi(temp, njr, hpf,lenfil,wcphi);
			}
			njr = njr/2;
		}
	}
	else {	/*continue in x-direction*/
		for(j=0; j<ellx-ell_min; j++) {
			for( k=0; k < njr; k++) {	/*transform rows*/
				unpackdouble(wc,njc,nr,k,temp);	/*copy one row into temp, wc: already transformed in column direction*/
				downlo(temp, njc, lpf,lenfil,templo);
				downhi(temp, njc, hpf,lenfil,temphi);
				packdouble(templo,njc/2,nr,k,wc);
				packdouble(temphi,njc/2,nr,k,&wc[njc/2*nr]);
			}
			njc = njc/2;
		}
	}
}

void unpackdouble(x,n,nc,k,y)
DOUBLE x[],*y;
int n,nc,k;
{  int i;
   for( i=0; i < n; i++)
   		*y++ = x[k+nc*i]; /*reads wc(i,k) i: x-position(col), k: y-position(row)*/
}
 
void packdouble(x,n,nc,k,y)
DOUBLE *x,y[];
int n,nc,k;
{  int i;
   for( i=0; i < n; i++)
		 y[k+nc*i] = *x++; /*write wc(i,k)*/
 }


void copydouble(x,y,n)
DOUBLE *x,*y;
int n;
{
   while(n--) *y++ = *x++;
 }
 
#include "downhi.c"
#include "downlo.c"
#include "mirrorfilt.c"
         


			
          
