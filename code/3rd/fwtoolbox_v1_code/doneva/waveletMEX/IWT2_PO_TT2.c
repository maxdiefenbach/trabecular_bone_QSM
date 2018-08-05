/*

 IWT2_PO.C	.MEX file corresponding to iwt2_po.m
		Periodized Wavelet Transform

  The calling syntax is:

			wc = iwt2_po(image,L,qmf)


  David Donoho
  Copyright (c) 1993  David Donoho
  All Rights Reserved

*/

#include <math.h>
#include "mex.h"
#include "wavelab.h"

void idpwt2(double *wc,int nr,int nc,int ellx,int elly,
   double *hpf,double *lpf,int lenfil,double *img,double *temp);

#define DOUBLE		double
#define INT			int

/* Input Arguments */

#define	WC_IN	prhs[0]
#define	Lx_IN	prhs[1]
#define	Ly_IN	prhs[2]
#define LPF_IN	prhs[3]


/* Output Arguments */

#define	Img_OUT	plhs[0]

INT nlhs, nrhs;
mxArray *plhs[], *prhs[];

void mexFunction(int nlhs,mxArray *plhs[],int nrhs,const mxArray *prhs[])
{
	DOUBLE	*hpf,*lpf;
	DOUBLE	*imp,*wcp,*tmp;
	unsigned int	m,n;
	int nr,nc,nn,lenfil,ellx,elly;
	mxArray *temp, *hpfmat;



	/* Check for proper number of arguments */

	if (nrhs != 4) {
		mexErrMsgTxt("IWT2_PO requires 3 input arguments.");
	} else if (nlhs != 1) {
		mexErrMsgTxt("IWT2_PO requires one output argument.");
	}


	/* Check the dimensions of signal.  signal can be n X 1 or 1 X n. */

	m  = mxGetM(WC_IN);
	n = mxGetN(WC_IN);
 
	nr = (int) m;
	nc = (int) n;


	/* Create a matrix for the return argument */

	Img_OUT = mxCreateDoubleMatrix(nr, nc, mxREAL);
	if(nr<nc)
		temp   = mxCreateDoubleMatrix(nc, 4, mxREAL);
	else 
		temp   = mxCreateDoubleMatrix(nr, 4, mxREAL);

	/* Assign pointers to the various parameters */

	imp = mxGetPr(Img_OUT);
	tmp = mxGetPr(temp);

	wcp = mxGetPr(WC_IN);
    ellx = floor ((mxGetPr(Lx_IN))[0] + .5);
	elly = floor ((mxGetPr(Ly_IN))[0] + .5);
    lpf = mxGetPr(LPF_IN);
    lenfil = (int) (mxGetM(LPF_IN) * mxGetN(LPF_IN));   /* should check this */
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

	idpwt2(wcp,nr,nc,ellx,elly,hpf,lpf,lenfil,imp,tmp);
	mxDestroyArray(temp);
	mxDestroyArray(hpfmat);
}


void idpwt2(wc,nr,nc,ellx,elly,hpf,lpf,lenfil,img,temp)
DOUBLE img[],hpf[],lpf[],wc[],temp[];
int  nr,nc,ellx,elly,lenfil;
{
        DOUBLE *wcplo,*wcphi,*templo,*temphi,*temptop;
        int k,j,njr,njc,ell_min;
        copydouble(wc,img,nr*nc);
		
		if(nr<nc) {
			templo = &temp[nc];
			temphi = &temp[2*nc];
			temptop = &temp[3*nc];
		}
		else {
			templo = &temp[nr];
			temphi = &temp[2*nr];
			temptop = &temp[3*nr];
		}
		
		if(ellx<elly)
			ell_min = ellx;
		else
			ell_min = elly;
		
		njr = nr;
		njc = nc;
		for( k=0; k < ellx; k++)
		   njc /= 2;
		for( k=0; k < elly; k++)
		   njr /= 2;

		/*perform 1-d transform in 1 direction*/
		if(ell_min == ellx)	{	/*continue in y-direction*/
			for(j=0; j<elly-ellx; j++) {
				for( k=0; k < njc; k++) {	/*transform columns*/
				   wcplo = &img[k*nr];
				   wcphi = &img[k*nr + njr];
				   copydouble(wcplo,temp,njr);
				   uplo(wcplo, njr, lpf,lenfil,templo);
				   uphi(wcphi, njr, hpf,lenfil,temphi);
				   adddouble(templo,temphi,njr*2,wcplo);
			   }
			   njr *= 2;
			}
		}
		else {	/*continue in x-direction*/
			for(j=0; j<ellx-elly; j++) {
				for( k=0; k < njr; k++){	/*transform rows*/
				   unpackdouble(img,njc,nr,k,templo);
				   unpackdouble(&img[njc*nr],njc,nr,k,temphi);
				   uplo(templo, njc, lpf,lenfil,temp);
				   uphi(temphi, njc, hpf,lenfil,temptop);
				   adddouble(temp,temptop,njc*2,temp);
				   packdouble(temp,njc*2,nr,k,img);
			   }
				njc *= 2;
			}
		}
		
		
		/* main 2D-transform*/
		for( j=0; j < ell_min; j++){
				for( k=0; k < 2*njr; k++){	/*transform rows*/
				   unpackdouble(img,njc,nr,k,templo);
				   unpackdouble(&img[njc*nr],njc,nr,k,temphi);
				   uplo(templo, njc, lpf,lenfil,temp);
				   uphi(temphi, njc, hpf,lenfil,temptop);
				   adddouble(temp,temptop,njc*2,temp);
				   packdouble(temp,njc*2,nr,k,img);
			   }

			   for( k=0; k < 2*njc; k++){	/*transform columns*/
				   wcplo = &img[k*nr];
				   wcphi = &img[k*nr + njr];
				   copydouble(wcplo,temp,njr);
				   uplo(wcplo, njr, lpf,lenfil,templo);
				   uphi(wcphi, njr, hpf,lenfil,temphi);
				   adddouble(templo,temphi,njr*2,wcplo);
			   }
			   njr *= 2;
			   njc *= 2;
		}
}

void unpackdouble(x,n,nr,k,y)
DOUBLE x[],*y;
int n,nr,k;
{  int i;
   for( i=0; i < n; i++)
   		*y++ = x[k+nr*i];	/*read im(i,k)*/
 }
 
void packdouble(x,n,nr,k,y)
DOUBLE *x,y[];
int n,nr,k;
{  int i;
   for( i=0; i < n; i++)
		 y[k+nr*i] = *x++;
 }


void copydouble(x,y,n)
DOUBLE *x,*y;
int n;
{
   while(n--) *y++ = *x++;
 }
 
void adddouble(x,y,n,z)
DOUBLE *x,*y, *z;
int n;
{
   while(n--) *z++ = *x++ + *y++;
}

 
#include "uphi.c"
#include "uplo.c"
#include "mirrorfilt.c"
