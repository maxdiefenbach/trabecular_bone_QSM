/* MEX function gr2dKB_mex() 
   Usage:
   m = gr2dKB_mex(k, d, w, grsize, kwidth, kernLUT);   
     k -- k-trajectory, Npts complex, scaled -0.5 to 0.5
     d -- k-space data, Npts
     w -- k-space weighting, Npts
     grsize -- grid size (m will be grsize X grsize)
     kwidth -- kernel width (one-sided)
     kernLUT -- sampled kernel lookup table

     m -- gridded k-space data
*/

#include "mex.h"
#include <math.h>

/* Computational subroutine */
/* copied from Brian's gridroutines.c and modified */
void gr2dKB_mex(
/* Gridding function that uses a lookup table for a separable
   convolution kernel, with linear interpolation. */
/* INPUT */
    double *kx,      /* Array of kx locations of samples. */
    double *ky,      /* Array of ky locations of samples. */
    double *s_real,  /* Sampled data, real part. */
    double *s_imag,  /* Sampled data, imag part. */
    double *dcf,     /* Density compensation factors. */
    int nsamples,    /* Number of k-space samples, total. */
    int gridsize,    /* Number of points in kx and ky in grid. */
    double convwidth,/* Kernel width, in grid points. */
    double *kerneltable, /* 1D array of convolution kernel values, starting
                            at 0, and going to convwidth. */
    int nkernelpts,  /* Number of points in kernel lookup-table */
/* OUTPUT */
    double *sg_real, /* OUTPUT array, real parts of data. */
    double *sg_imag  /* OUTPUT array, imag parts of data. */
    )
{    
    int kcount;             /* Counts through k-space sample locations */
    int gcount1, gcount2;	/* Counters for loops */
    double kwidth;			/* Conv kernel width, in k-space units */
    double dkx,dky;         /* Delta in x, y for kernel calc.*/
    int ixmin,ixmax,iymin,iymax;	/* min/max indices that current k may affect*/
    int kernelind;			/* Index of kernel value, for lookup. */
    double fracdk;			/* Fractional part of lookup index. */
    double fracind;			/* Fractional lookup index. */    
    double kernx, kerny;    /* Interpolated Kernel value from LUT */
    int gridsizesq;			/* Square of gridsize */
    int gridcenter;         /* Index of center of grid. */
    int gpt;                /* Target point on k-sp grid */  

    gridcenter = gridsize/2;
    kwidth = convwidth/(double)(gridsize);	/* Width of kernel, in k-space units. */
    
    /* ========= Init Output Points ========== */
    gridsizesq = gridsize*gridsize;
    for(gcount1 = 0; gcount1 < gridsizesq; gcount1++)
    {
        sg_real[gcount1] = 0;
        sg_imag[gcount1] = 0;
    }

    /* ========= Loop Through k-space Samples ========= */
    for(kcount = 0; kcount < nsamples; kcount++)
    {
        /* ----- Find limit indices of grid points that current
             sample may affect (and check they are within grid) ----- */
        ixmin = (int) ((kx[kcount] - kwidth)*gridsize +gridcenter);
        if (ixmin < 0) ixmin=0;
        ixmax = (int) ((kx[kcount] + kwidth)*gridsize +gridcenter)+1;
        if (ixmax >= gridsize) ixmax=gridsize-1;
        iymin = (int) ((ky[kcount] - kwidth)*gridsize +gridcenter);
        if (iymin < 0) iymin=0;
        iymax = (int) ((ky[kcount] + kwidth)*gridsize +gridcenter)+1;
        if (iymax >= gridsize) iymax=gridsize-1;

        for(gcount1 = ixmin; gcount1 <= ixmax; gcount1++)
        {
            dkx = (double)(gcount1-gridcenter)/(double)gridsize  - kx[kcount];
            dkx = fabs(dkx); /* in k-space units */
            
            for(gcount2 = iymin; gcount2 <= iymax; gcount2++)
            {
                dky = (double)(gcount2-gridcenter)/(double)gridsize - ky[kcount];
                dky = fabs(dky); /* in k-space units */
                
                if(dkx<=kwidth && dky<=kwidth) /* sample affects this grid point */
                {
                    /* Separable Kernel */
                    /* KX: Find index in kernel lookup table */
                    fracind = dkx/kwidth*(double)(nkernelpts-1);
                    kernelind = (int)fracind;
                    fracdk = fracind-(double)kernelind;
                    /* KX: Linearly interpolate in kernel lut */
                    kernx = kerneltable[kernelind]*(1-fracdk)+
                            kerneltable[kernelind+1]*fracdk;
                    
                    /* KY: Find index in kernel lookup table */
                    fracind = dky/kwidth*(double)(nkernelpts-1);
                    kernelind = (int)fracind;
                    fracdk = fracind-(double)kernelind;
                    /* KY: Linearly interpolate in kernel lut */
                    kerny = kerneltable[kernelind]*(1-fracdk)+
                            kerneltable[kernelind+1]*fracdk;
                    
                    gpt = gcount1*gridsize + gcount2;
                    sg_real[gpt] += kernx*kerny * s_real[kcount] * dcf[kcount];
                    sg_imag[gpt] += kernx*kerny * s_imag[kcount] * dcf[kcount];
                }
            }
        }
    }    

    return;
}    


/* The gateway routine */
/* copied from Brian's gridlut_mex.c and modified */
void mexFunction(int nlhs, mxArray *plhs[],
                 int nrhs, const mxArray *prhs[])
{
    /* Input */
    double *kx;          /* k-space trajectory x locations.	*/
    double *ky;          /* k-space trajectory y locations.	*/
    double *s_real;		 /* real-part of data samples.	*/
    double *s_imag;		 /* imaginary-part of data samples. */
    double *dcf;		 /* (pre) density correction factor at each traj loc.*/    
    int nsamples;		 /* number of data samples (or trajectory points) */
    int gridsize;		 /* size of grid, in samples. ie gridsize x gridsize */
    double convwidth;	 /* width of convolution kernel for gridding. */
    double *kerneltable; /* lookup-table for linearly-interpolated kernel.  */
    int nkernelpts;		 /* number of points in kernel lookup table. */
    /* Output */
    double *sg_real;	 /* real-part of gridded data sample. */
    double *sg_imag;	 /* imaginary-part of gridded data sample. */

    /*mxArray *dens_mat;*/

    /* Parse Input */
    /* ---------------------------------------------------------------- */
    /* Samples may be passed as 1xN, Nx1 or 2D array */
    nsamples = mxGetM(prhs[0]) * mxGetN(prhs[0]);    
    
    kx = mxGetPr(prhs[0]);		/* Get kx locations */
    ky = mxGetPi(prhs[0]);		/* Get ky locations */
    s_real = mxGetPr(prhs[1]);	/* Get real parts of data samples. */
    s_imag = mxGetPi(prhs[1]);	/* Get imaginary parts of data samples. */
    dcf = mxGetPr(prhs[2]);	    /* Get density correction factors. */
    gridsize = (int)(*mxGetPr(prhs[3])); /* Get grid size. */

    if (nrhs > 4)
        convwidth = *mxGetPr(prhs[4]);	 /* Get conv width */
    else
        convwidth = 1.0;				 /* Assign default */

    kerneltable = mxGetPr(prhs[5]);		            /* Get kernel table.*/
    nkernelpts = mxGetM(prhs[5]) * mxGetN(prhs[5]);	/* and # points. */
    /* ---------------------------------------------------------------- */    

    /* Init Output */
    /* ---------------------------------------------------------------- */    
    plhs[0] = mxCreateDoubleMatrix(gridsize,gridsize,mxCOMPLEX);
    sg_real = mxGetPr(plhs[0]);
    sg_imag = mxGetPi(plhs[0]);    
    /* ---------------------------------------------------------------- */    
    
    gr2dKB_mex(kx, ky, s_real, s_imag, dcf, nsamples,
	           gridsize, convwidth, kerneltable, nkernelpts,
               sg_real, sg_imag);
    
    return;
}
