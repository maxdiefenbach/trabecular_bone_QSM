function A = fat_basis(te,NDB,H2O,Tesla)

% Function that produces the fat-water matrix.
%
% Inputs:
%   te (echo times in ms)
%   NDB (number of double bonds - NB. per trigliceride molecule)
%   H2O (water freq in ppm)
%   Tesla (field strength)
% Ouputs:
%   A (matrix [water fat] basis vectors)
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

% do checks here
ne = numel(te);
te = reshape(te,ne,1);
if nargin<2
    NDB = 3;
    disp(['Warning: fat_basis.m assuming NDB=' num2str(NDB)])
end
if nargin<3
    H2O = 4.7;
    disp(['Warning: fat_basis.m assuming H2O=' num2str(H2O)])
end
if nargin<4
    Tesla = 3;
    disp(['Warning: fat_basis.m assuming Tesla=' num2str(Tesla)])
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
a = a/(134+42*CL-2*NDB);

% time evolution matrix (te in ms)
water = repmat(awater,ne,1);
fat = zeros(ne,1);
larmor = 42.576e6*Tesla*1e-6;
for j = 1:numel(d)
    freq = larmor*(d(j)-H2O); % Hz relative to water
    fat = fat + a(j)*exp(-2*pi*i*freq*te*1e-3); % negative sign for conj?
end
A = [water fat];