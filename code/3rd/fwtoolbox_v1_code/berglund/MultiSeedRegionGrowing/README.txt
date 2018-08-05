
This is a MATLAB implementation of the three-point Dixon method described in:

Johan Berglund, Lars Johansson, HŒkan Ahlstršm, and Joel Kullberg
'Three-point Dixon Method Enables Whole-Body Water and Fat Imaging of Obese Subjects'
Magnetic Resonance in Medicine 63 (6): 1659-1668 (2010)

Using some generalizaed equations described in:

Johan Berglund, HŒkan Ahlstršm, Lars Johansson, and Joel Kullberg
'Closed-form solution for the three-point Dixon method with advanced spectrum modeling'
Submitted to Proc. ISMRM 2011

Please cite these papers if you use the code.
The code is written by Johan Berglund, Uppsala University, Sweden.

A subroutine of the algorithm is implemented in c++ for efficiency. The c++ file (RG.cpp) must be compiled before MATLAB can use it. This can be done in MATLAB by navigating to the source code folder, then typing in the MATLAB command window:

	mex RG.cpp

If this doesn't work, you may have to specify a compiler by typing: mex -setup