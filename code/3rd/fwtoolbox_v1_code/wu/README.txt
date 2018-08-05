=====================================================================

CCRingsFW is a collection of MATLAB scripts for
fat/water Separation using concentric rings k-space sampling.

Distributed as part of the ISMRM Fat-Water Imaging Toolbox.

Written by Holden H. Wu, 2011
(c) Board of Trustees, Leland Stanford Junior University

Please send comments and suggestions to
holdenhwu@stanford.edu

=====================================================================

Modified:
  2012/02/15 now combining coils inside fw_k2cs1i_ccr_gcut_wu()
  2012/01/16 now using ccr_k2im()
             included templates for general image-space algorithms
  2012/01/12 updated fdLS and gcut code
  2011/12/28 updated README.txt	       

=====================================================================

***** Setup instructions *****

Make sure the following are in your MATLAB path:
CCRingsFW
CCRingsFW/ccrlib/
CCRingsFW/grid2D/

Compile gr2dKB_mex.c (see info below).

To use the graph cut approach, please make sure Diego Hernando's
directories are also in your MATLAB path.


***** Notes *****

The concentric rings enable a 2-fold reduction in the number of 
readouts (excitations) compared to 2D Cartesian encoding [1]. 

Chemical shift information is encoded by a time-efficient retracing
design that traverses each ring over multiple revolutions [1], 
similar to a multi-echo 2DFT acquisition. As the rings always rotate 
in the same direction, there is no need to align "even/odd" echoes.

The (k,t)-space spectroscopic algorithm [1] is not included in 
this release, but is planned for the next version.


***** Package contents *****

CCRingsFW/

  ccrings_FW_fdLS.m
    Concentric rings fat/water separation using a 
    frequency-demodulated iterative least-squares algorithm [1,2].

  ccrings_FW_gcut.m
    Concentric rings fat/water separation using Hernando's
    graph cut algorithm [3] for water recon and frequency-demodulated 
    least squares [1,2] for fat recon.

  ccrings_FW_template.m
    Template for concentric rings fat/water separation using any
    general image-space algorithm. See ccrings_FW_gcut.m for an example.

  README.txt
    What you are reading right now.

CCRingsFW/ccrlib/

  fw_k2cs0i_ccr_fdLS_wu.m
    Main function for freq-demodulated iterative LS [1].
    Type 'help fw_k2cs0i_ccr_fdLS_wu' for more information.

  fdIDEAL.m
    Modified version of IDEAL [2] that accounts for frequency
    demodulation of the k-space data. This simple implementation
    uses only low-pass filtering for field map processing, and can
    be improved with more sophisticated processing algorithms.
    Type 'help fdIDEAL' for more information.

  fw_k2cs1i_ccr_gcut_wu.m
    Main function for the graph cut approach [3].
    Type 'help fw_k2cs1i_ccr_gcut_wu' and consult Hernando's files
    for more information.

  fw_tttttt_ccr_imsp.wu.m
    Template for using any image-space algorithm for the concentric rings.
    Type 'help fw_tttttt_ccr_imsp_wu' for more information.    
    fw_k2cs1i_ccr_gcut_wu.m is an example.

  ccr_k2im.m
    Split multi-rev rings k-space data into individual revs,
    demodulate at the desired frequency, and reconstruct using
    2D gridding.

  getsubset.m
    Utility function for extracting a subset of the k-space data 
    for reconstruction.

CCRingsFW/grid2D/

  grecon2d.m
    Main function for 2D gridding reconstruction.

  gr2dKB_mex.c
    Core MEX function for 2D Kaiser-Bessel gridding.
    Precompiled for Linux 32/64 bit, Windows 32/64 bit, and Mac OSX.
    Please recompile for your specific platform if necessary.

CCRingsFW/ccrdata/

  20090520_3DHeadB_3sl.mat
    3D stack-of-rings IR-SPGR head scan [4].
    3 revolutions for each ring.
    Only including 3 of 180 1-mm slices.
    Already Fourier-transformed in the slice direction.

  20110819_SAX.mat
    Concentric rings 2D GRE cardiac short-axis cine scan [5].
    3 revolutions for each ring.
    This a fully sampled dataset, no k-t BLAST acceleration.

  20110819_LAX.mat
    Concentric rings 2D GRE cardiac long-axis cine scan [5].
    3 revolutions for each ring.
    This a fully sampled dataset, no k-t BLAST acceleration.


***** References *****

[1] Wu HH, Lee J, Nishimura DG. Fat/Water Separation Using 
a Concentric Rings Trajectory. Magn Reson Med 2009;61:639-649.

[2] Reeder SB, Wen Z, Yu H, Pineda AR, Gold GE, Markl M, Pelc NJ.
Multicoil Dixon chemical species separation with an iterative 
least squares estimation method. Magn Reson Med 2004;51:35¡V45.

[3] Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat 
separation in the presence of large field inhomogeneities using a 
graph cut algorithm. Magn Reson Med 2010 63:79-90.

[4] Wu HH, Nishimura DG. 3D Magnetization-Prepared Imaging Using 
a Stack-of-Rings Trajectory Magn Reson Med 2010;63:1210-1218.

[5] Wu HH, Shin T, Nishimura DG, McConnell MV. Rapid Fat-Water-
Separated Cardiac Cine Imaging Using Concentric Rings and k-t BLAST, 
In: Proceedings of the ISMRM 19th Annual Meeting, Montreal, Canada, 
2011, p4366.

=====================================================================
