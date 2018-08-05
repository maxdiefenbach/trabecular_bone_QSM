Author: Samir Sharma
Date: December 20, 2011


This folder contains functions for non-accelerated image-space water-fat separation.  The main functions contained are:

* fw_i2cm0i_3pluspoint_sharma.m
Water-fat separation using a restricted subspace field map estimation
- Reference: Sharma SD, Hu HH, Nayak KS. Accelerated water-fat imaging using restricted subspace field map estimation and compressed sensing. Magn Reson Med (in early view)


* utils/coilCombine.m
Coil combination for image sequences
- Reference: Walsh DO, Gmitro AF, Marcellin MW. Adaptive reconstruction of phased array MR imagery. Magn Reson Med 2000;43:682-690


* utils/calculate_chemical_shift_encoding_matrix
Calculate water-fat modeling matrix using a six-peak fat spectrum. 


* utils/fnlCg.m
Conjugate gradient implementation using for updating the field map estimate.  A modification of the function provided by Lustig.
- Reference: Lustig M, Donoho D, Pauly J. Sparse MRI: the application of compressed sensing for rapid MR imaging. Magn Reson Med 2007;58:1182-1195. 



The top-level folder contains a test script that shows how to call fw_i2cm0i_3pluspoint_sharma.m 

* test_sharma_20111220.m (test restricted subspace field map estimation on randomly-chosen dataset from Peter Kellman)



Some issues: 

* So far, code has been tested succesfully on:
1) Windows - 64bit - R2010b
