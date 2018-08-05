Author: Holger Eggers
Date: January 11, 2012

This folder contains one sample dataset for image-space fat-water separation 
and field inhomogeneity correction. 

The mat-file consists of three single echo images for a single, coronal slice 
through the abdomen of a volunteer from a 3D spiral acquisition. The data are 
stored as data.images(row,col,slice,coil,echo), where #row = #col = 320, 
#slice = 1, #coil = 1 (already combined), #echo = 3. For relevant acquisition 
parameters, see the provided test program.

The field strength is data.FieldStrength = 1.5, and the echo times are data.TE 
= [0.001, 0.0025, 0.004] s.
