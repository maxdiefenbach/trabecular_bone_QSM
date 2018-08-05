% Function: getFatWaterSpectrum_massScale_fatByDoubleBonds
%
% Description: Function that produces the fat-water spectrum (including multipeak
% fat frequencies and relative amplitudes), ready to be fed to a water/fat 
% decomposition algorithm (output of this function can be directly assigned to the 
% "species" field of the algoParams parameter structure.)
%
% Inputs:
%   NDB (number of double bonds - NB. per trigliceride molecule)
%   H2O (water freq in ppm)
% Ouputs:
%   species: species structure array, defining frequencies and relative amplitudes of water and fat
%      This array will have content similar to:
%      - species(1).name = 'water' % Water
%      - species(1).frequency = [0] 
%      - species(1).relAmps = [1]   
%      - species(2).name = 'fat' % Fat
%      - species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%      - species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
% 
% References:
%
% In vivo characterization of the liver fat 1H MR spectrum.
% G Hamilton, T Yokoo, M Bydder, I Cruite, ME Schroeder, CB Sirlin, MS Middleton.
% NMR in Biomedicine 2011, Volume 24, Issue 7, pages 784–790.
%
% Mapping the double bonds in triglycerides.
% M Bydder, O Girard, G Hamilton.
% Magnetic Resonance Imaging 2011, Volume 29, Issue 8, Pages 1041-1046 
%
% Author: Mark Bydder
% Date created: November 14, 2011
% Date last modified: November 14, 2011


function species = getFatWaterSpectrum_massScale_fatByDoubleBonds(NDB,H2O)

% do checks here
if nargin<1
    NDB = 3;
    disp(['Warning: fat_basis.m assuming NDB=' num2str(NDB)])
end
if nargin<2
    H2O = 4.7;
    disp(['Warning: fat_basis.m assuming H2O=' num2str(H2O)])
end

% fat chemical shifts in ppm
d = [5.29 5.19 4.2 2.75 2.2 2.02 1.6 1.3 0.9];

% heuristic formulas relating CL and NDDB to NDB
CL = 16.8+0.25*NDB;
NDDB = 0.093*NDB^2;

% no. protons per molecule
awater = 2;
a(1) = NDB*2;
a(2) = 1;
a(3) = 4;
a(4) = NDDB*2;
a(5) = 6;
a(6) = (NDB-NDDB)*4;
a(7) = 6;
a(8) = (CL-4)*6-NDB*8+NDDB*2;
a(9) = 9;

% the above code counts no. molecules. i.e. 1 unit of water = 2
% protons and 1 unit of fat = (2+6*CL-2*NDB) protons. now scale
% so 18 "units" of water = 2 protons and (134+42*CL-2*NDB) "units"
% of fat = (2+6*CL-2*NDB) protons, i.e. so that w and f count mass.
awater = awater/18;
afat = a/(134+42*CL-2*NDB);


% DH*: Let us normalize so awater = 1
afat = afat./awater;
awater = 1;

% Create species structure array
species(1).name = 'water';
species(1).frequency = 0;
species(1).relAmps = awater;
species(2).name = 'fat';
species(2).frequency = d - H2O;
species(2).relAmps = afat;
