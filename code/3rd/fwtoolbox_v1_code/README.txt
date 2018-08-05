%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%  Fat-Water Toolbox, v1  %%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%%%%%%%%%%%%%%%%      Important        %%%%%%%%%%%%%%%%%%%%%%%%

THIS SOFTWARE IS FOR INVESTIGATIONAL PURPOSES. NOT INTENDED FOR CLINICAL USE.


%%%%%%%%%%%%%%%%%%%%       Warning         %%%%%%%%%%%%%%%%%%%%%%%%

Please note that some of the functions included in this toolbox are
implemented in MEX. These have been compiled to work on a variety of
platforms, but may require re-compilation for your particular
operating system/Matlab version. 


%%%%%%%%%%%%%%%%%%%%     Introduction      %%%%%%%%%%%%%%%%%%%%%%%%

This is a toolbox of Matlab functions containing a set of recently proposed 
methods and algorithms related to fat-water separation in MRI. The goal of 
the toolbox is to provide easy access to these advanced algorithms. Also, we 
aim to define a unified interface for these various algorithms (same structure 
of inputs/outputs for different Matlab functions), to simplify their use. 
Additionally, the toolbox includes a large number of test datasets acquired using 
various acquisition parameters in different anatomical regions (cardiac, abdomen, 
liver, thigh,...)

This toolbox has been developed as an initiative of the ISMRM Fat-Water MRI Workshop 
(February 19-22, Los Angeles, CA, USA). The initiative has been coordinated by D. Hernando 
(dhernando@wisc.edu) and Houchun H. Hu (houchunh@usc.edu). The toolbox includes 
code/datasets by the following contributors:

Johan Berglund
Emily Bice
Mark Bydder
Mariya Doneva
Holger Eggers 
Diego Hernando
Houchun H. Hu
Yun Jiang
Peter Kellman
Wenmiao Lu
Angel Pineda
Samir Sharma
Jeff Tsao
Holden Wu

Note that some of the functions require compilation of the Mex source code using 
either C or C++. Although we have attempted to include compiled versions of these 
functions for a variety of platforms (including Linux, Mac and Windows), your 
platform/Matlab version may require compilation. Please see each individual 
contributor's folder for directions. 



%%%%%%%%%%%%%%%%%%%%         Usage         %%%%%%%%%%%%%%%%%%%%%%%%

To test the methods contained in this toolbox, either run the test Graphical 
User Interface in this same folder ("testGUI"), or directly access the desired 
code (see list below -- most folders include individual test and README files). 




%%%%%%%%%%%%%%%%%%    Basic code structure     %%%%%%%%%%%%%%%%%%%%

Most functions have two input structures (dataParams and algoParams) and one
output structure (outParams). For instance, for image-based fat-water 
separation algorithms, dataParams contains the following fields:

  dataParams.FieldStrength (in Tesla)
  dataParams.TE (echo times in seconds; vector of length nte)
  dataParams.PrecessionIsClockwise (1/-1)
  dataParams.images (the actual data, array of dimensions nx X ny X nz X ncoils X nte)

The algoParams structure contains information to solve the separation problem, including:

  algoParams.species structures (typically two of them, i=1 for water and i=2 for fat)
  algoParams.species(i).name  
  algoParams.species(i).frequency (in ppm)
  algoParams.species(i).relAmps (typically normalized to add up to 1)

plus any number of algorithm-specific parameters (regularization parameters, use of R2*, etc)

The outParams structure contains the resulting separated water/fat images, plus possibly 
fieldmap estimates, R2* maps,...

  outParams.species(i).amps (map of amplitudes of species "i", where typically i=1 for water, i=2 for fat)
  outParams.fieldmap (B0 field map in Hz; not all algorithms return one)
  outParams.r2starmap (R2* map in 1/s; not all algorithms return one)


* File names for reconstruction algorithms:

For fat-water separation functions, the algorithm type is "encoded" in the function name: 

for instance
   fw_i2cm1i_3pluspoint_hernando_graphcut
has the following properties:

 - i: image-space ('k' for k-space)
 - 2: 2D (processes a single slice; '3' for 3D)
 - c: complex-valued data processing ('m' for magnitude processing, 'x' for mixed)
 - m: multi-peak fat model ('s' for single-peak)
 - 1: for single-R2* model ('0' for no-R2*)
 - i: for independent initial phase for water and fat ('c' for common phase)

 - 3pluspoint: works with any number of echoes >= 3 (some choices much better than others!)
               ('2point' and '3point' would require exactly 2 or 3 points, respectively)

 - hernando: authorname

 - graphcut: algorithm keyword

%%%%%%%%%%%%%%%%%%%%  Available functions  %%%%%%%%%%%%%%%%%%%%%%%%



   %%%%%%%%%%% Image-space fat-water separation %%%%%%%%%%%


* Two-point fat-water separation with flexible echo times (J. Berglund)

berglund/FlexibleTwoPoint/fw_i3cm0c_2flexiblepoint_berglund.m (May require Mex C++ compilation)


* Three-point fat-water separation with multi-seed region growing (J. Berglund)

berglund/MultiSeedRegionGrowing/fw_i3cm0i_3point_berglund.m (May require Mex C++ compilation)


* Multi-point fat-water separation with R2* using a graphcut field map estimation algorithm (D. Hernando)

hernando/graphcut/fw_i2cm1i_3pluspoint_hernando_graphcut.m


* Multi-point fat-water separation with R2* and phase error correction using mixed magnitude/complex fitting (D. Hernando) 

hernando/mixed_fitting/fw_i2xm1c_3pluspoint_hernando_mixedfit.m (May require C compilation/instalation of Matlab-BGL)


* Multi-point fat-water separation with regularized field map using an optimization transfer algorithm (D. Hernando based on algorithm by Huh et al)

hernando/descent/fw_i2cm0i_3plusploint_hernando_optimtransfer.m


* Fat-water separation with multi-peak fat using region-growing (H. Hu based on an algorithm by H. Yu et al)

hu/fw_i2cm0c_3pluspoint_RGsinglecoil_hu.m


* Three-point fat-water separation with a multi-resolution golden-section search algorithm (W. Lu)

lu/fw_3point_wm_goldSect.m


* Multi-point ater-fat separation using a restricted subspace field map estimation (S. Sharma)

sharma/fw_i2cm0i_3pluspoint_sharma.m


* Multi-point fat-water separation using a hierarchical field map estimation (J. Tsao, Y. Jiang)

tsao_jiang/fw_i2cm0c_3pluspoint_tsaojiang.m


* Three-point fat-water separation using a hierarchical field map estimation (J. Tsao, Y. Jiang)

tsao_jiang/fw_i2cs0c_3point_tsaojiang.m


 
      %%%%%%%%%%% k-space fat-water separation %%%%%%%%%%%


* Compressed-sensing based fat-water separation (M. Doneva)

doneva/fw_k2cs0i_3point_doneva_cswf.m (May require Mex C compilation)


* Concentric circles acquisition fat-water separation  (H. Wu)

wu/ccrlib/fw_k2cs0i_ccr_fdLS_wu.m (No R2*; may require Mex C compilation)
wu/ccrlib/fw_k2cs1i_ccr_gcut_wu.m (With R2*; may require Mex C compilation)


* Spiral (deblurring) fat-water imaging (H. Eggers)

eggers/cpr.m (Actually used as deblurring postprocessing for spiral acquisitions processed using image-space algorithm; see example in ./eggers/)



       %%%%%%%%%%% Modeling/analysis/misc %%%%%%%%%%%


* Calculation of fat spectral model based on the number of double bonds (M. Bydder)

bydder/getFatWaterSpectrum_massScale_fatByDoubleBonds.m


* Generation of synthetic image-space datasets with arbitrary parameters (D. Hernando)

hernando/create_synthetic/createSynthetic_imageSpace.m


* Generation of synthetic k-space datasets (M. Doneva)

doneva/create_synthetic/createSynthetic_kSpace.m


* Graphical user interface for many of the functions in this toolbox (most functions except k-space reconstructions; D. Hernando)

gui/fw_gui.m (Accessible from this folder using testGUI.m)


* Number of Signal Averages (NSA) calculation and mapping (A. Pineda, E. Bice)

pineda_bice/fw_nsa_2pluspoint.m (Calculate theoretical NSA over a range of fat/water fractions)
pineda_bice/fw_mcsim_2pluspoint.m (Measure simulated -- Monte-Carlo -- NSA)
pineda_bice/fw_nsa_map.m (Calculate NSA maps based on water/fat images and echo times)


* Fat-fraction calculation (D. Hernando based on method by C.Y. Liu et al)

hernando/common/computeFF.m




%%%%%%%%%%%%%%%%%%%%    Acknowledgements     %%%%%%%%%%%%%%%%%%%%%%
 
This package utilizes routines from MatlabBGL, a Matlab wrapper
around the Boost graph library written by David F. Gleich.  If you
publish anything utilizing the MaxFlow/MinCut routines (D. Hernando's 
graphcut algorithm in this package), be sure to cite both this 
package and MatlabBGL.

