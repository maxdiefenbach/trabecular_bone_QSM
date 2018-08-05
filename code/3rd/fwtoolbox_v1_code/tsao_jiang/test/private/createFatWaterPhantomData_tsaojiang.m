function [imDataParams,trueParams] = createFatWaterPhantomData_tsaojiang(MatrixSize,FieldStrength,InitTEInMs,singlepeak,params);

% Modify from Diego Hernando's code for generating phantom data
%
% Only 3 TEs is allowed,and TE is calculated based on the fieldstrength to 
% meet the requirement of Hierarchical IDEAL fat-water separation 
% algorithm.  The InitTE is for reference only, usually can't be achieved. 
%
% Inputs:
%   MatrixSize - image matrix
%   FieldStrength - in unit of TESLA
%   InitTE - minimum echo time in unit of MILLISECOND; for reference only, 
%           usually can't be achieved.
if nargin < 1, MatrixSize = []; end
if nargin < 2, FieldStrength = []; end
if nargin < 3, InitTEInMs = []; end %Unit is MILLISECOND
if nargin < 4, singlepeak = []; end
if nargin < 5, params = []; end
if isempty(MatrixSize),    MatrixSize = 128; end
if isempty(FieldStrength), FieldStrength = 1.5; end
if isempty(InitTEInMs),    InitTEInMs = 1.0; end %Unit is MILLISECOND
if isempty(singlepeak),    singlepeak = 1; end

sx = MatrixSize;
sy=sx;
imtest = phantom(sx); % create a Shepp-Logan phantom based on provided 
                      % matrix size
threshold = 0.8;

% [echotimes] = idealgradechotimes(MinEchoTimeInMs,fieldsstrength),
% IDEAL water fat suppression - gradient echo
%  - corrected based on Pineda et al. MRM. 54: 625-635. 2005
%    Mid-echo should be 90-phase offset

ppmdif = 3.4;
TEInSec = idealgradechotimes(InitTEInMs, FieldStrength, ppmdif);

trueParams.species(1).amps = imtest.*(imtest<threshold); % Water
trueParams.species(2).amps = imtest.*(imtest>=threshold); % Fat

if isempty(params), params = randn(8,1); end  % Generate new params if not given
fprintf('Phantom parameters:\nparams=[');
fprintf(' %f',params); fprintf(' ];\n');
[X,Y] = meshgrid(linspace(-1,1,sx),linspace(-1,1,sy));
trueParams.fieldmap = 20*params(1)*ones(sx,sy) + 40*params(2)*X + 40*params(3)*Y + 100*params(4)*X.^2 + 100*params(5)*Y.^2  + 100*params(6)*X.*Y.^2 + 100*params(7)*X.^3  + 100*params(8)*Y.^3;
trueParams.r2starmap = 80 - 40*imtest;
trueParams.ParamsToGeneratePhantom = params;

imDataParams0.TE = TEInSec/1000.0'; % ms -> second
imDataParams0.FieldStrength = FieldStrength;

clear algoParams
%% Set recon parameters
% General parameters
algoParams.species(1).name = 'water';
algoParams.species(1).frequency = 0;
algoParams.species(1).relAmps = 1;


algoParams.species(2).name = 'fat';
if singlepeak,
  % Set fat as single peak model
  algoParams.species(2).frequency = [-3.40];
  algoParams.species(2).relAmps = [1];
else
  % Set fat as multi-peak model
  algoParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
  algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
end

% Simulate data
imDataParams = createSynthetic_imageSpace( imDataParams0, algoParams, trueParams );
