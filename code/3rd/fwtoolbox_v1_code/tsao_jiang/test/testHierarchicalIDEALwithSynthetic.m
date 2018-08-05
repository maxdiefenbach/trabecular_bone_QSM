% testHierarchicalIDEALwithSynthetic.m is a function for creating synthetic
% data to test Hierarchical IDEAL algorithm
%
% As Hierarchical IDEAL algorithm requires TEs to be optimized echo times
% based on Pineda et al. MRM. 54: 625-635. 2005, the synthetic data will be
% generated with 3 optimized TEs 

% Input:
%       
%       FieldStrength - unit in Tesla
%       CoilMode   - 'single' - single coil
%                    'multiple' - generate using 8-channel sensitivity maps



% Jeff Tsao & Yun Jiang --Sept 16, 2011

function [outParams] = testHierarchicalIDEALwithSynthetic(FieldStrength,CoilMode,GenerateSingleFatPeak,PhantomParams);
if nargin<1, FieldStrength = []; end;
if nargin<2, CoilMode = []; end;
if nargin<3, GenerateSingleFatPeak = []; end
if nargin<4, PhantomParams=[]; end
if isempty(FieldStrength), FieldStrength = 1.5;end;
if isempty(CoilMode),
  str1 = 'Single coil';
  str2 = '8-coil array';
  button=questdlg({'Testing hierarchical IDEAL with synthetic data...','Which coil to generate test?'},'Hierarchical IDEAL',str1,str2,str2);
  if isequal(button,str1),
    CoilMode = 'single';
  elseif isequal(button,str2),
    CoilMode = 'multiple';
  else  % cancelled
    outParams = [];
    return;
  end
  clear str1 str2 button;
end;
if isempty(GenerateSingleFatPeak), 
  str1 = 'Single peak';
  str2 = 'Multi-peak';
  button=questdlg({'Which signal model for fat?'},'Hierarchical IDEAL',str1,str2,str2);
  if isempty(button), return; end % Cancelled;
  if     isequal(button,str1), GenerateSingleFatPeak = 1;
  elseif isequal(button,str2), GenerateSingleFatPeak = 0;
  end
end
if GenerateSingleFatPeak,
  datadescription = 'single-peak data';
  Version = 'Original'; % Choose original single-peak version to run
else
  datadescription = 'multi-peak data';
  Version = 'Flexible'; % Choose new multi-peak version to run
end

[BASEPATH,tmpfile] = fileparts(mfilename('fullpath'));clear tmpfile;
tmp = BASEPATH; addpath(tmp); fprintf('Adding to path: %s\n',tmp); clear tmp;

MatrixSize = 128;
  
while 1
  switch lower(CoilMode),
    case 'single'
        [imDataParams,trueParams] = createFatWaterPhantomData_tsaojiang(MatrixSize,[],[],GenerateSingleFatPeak,PhantomParams);
        break;
    case 'multiple'
        [imDataParams,trueParams] = createFatWaterPhantomData_tsaojiang(MatrixSize,[],[],GenerateSingleFatPeak,PhantomParams);
        load coilmap.mat; % 8 channel coil arry sensivity map with matrix 128*128
      
        imDataParams.images = repmat(imDataParams.images,[1,1,1,size(profile_est,3),1]); % [nx,ny,nz,coils,TEs]
        profile_est = reshape(profile_est,[size(profile_est,1),size(profile_est,2),1,size(profile_est,3)]);
        for n =1:size(imDataParams.images,5),
            imDataParams.images(:,:,:,:,n) = imDataParams.images(:,:,:,:,n).*profile_est;
        end;clear n;
        break;
    otherwise
        disp('Please choose right coil mode Single or Multiple...and try again.');
        CoilMode=questdlg({sprintf('Incorrect value for CoilMode ( %s )',CoilMode), ...
                        ' ',...
                        'single - Single coil',...
                        'multiple - Phased array'},...
                       'Hierarchical IDEAL version', ...
                       'single','multiple','multiple');
        if isempty(CoilMode), return; end % Pressed cancel
  end
end;

% Set algoParams
algoParams.MinFractSizeToDivide = 0.01;
algoParams.MaxNumDiv = 7;

% Select version
if isempty(Version),
  Version=questdlg({'Which version do you want to run?', ...
                        ' ',...
                        'Original version - only 3 echoes, requires optimal TE settings',...
                        'Flexible version - 3 or more echoes with uniform spacing'},...
                       'Hierarchical IDEAL version', ...
                       'Original','Flexible','Flexible');
  if isempty(Version), return; end % Pressed cancel
end

% run
if isequal(Version,'Original'), % Original
  recondescription = 'single-peak recon';
  outParams = fw_i2cs0c_3point_tsaojiang(imDataParams,algoParams);
elseif isequal(Version,'Flexible'), % New flexible
  recondescription = 'multi-peak recon';
  outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams,algoParams);
  outParams.water = outParams.species(1).amps;
  outParams.fat   = outParams.species(2).amps;
end

if strcmpi(CoilMode,'multiple'),
  % Normalized by RMS coil map
  rmsCoil = sqrt(mean(abs(profile_est).^2,4));
  outParams.water = outParams.water./rmsCoil;
  outParams.fat = outParams.fat./rmsCoil;
  outParams.fiterror = outParams.fiterror./rmsCoil;
  clear rmsCoil;
end

if strcmpi(CoilMode,'single'),
  figure('name',sprintf('Single Coil (%s, %s)',datadescription,recondescription));
else
  figure('name',sprintf('Multiple Coils (%s, %s)',datadescription,recondescription));
end

subplot(3,4,1);imagesc(abs(trueParams.species(1,1).amps));caxis([0 1]);title('True WATER image');colorbar;axis image;
subplot(3,4,2);imagesc(abs(trueParams.species(1,2).amps));caxis([0 1]);title('True FAT image');colorbar;axis image;
subplot(3,4,3);imagesc(abs(trueParams.r2starmap));title('True R2* map');r2star_caxisval = caxis; colorbar;axis image;
subplot(3,4,4);imagesc(abs(trueParams.fieldmap));title('True field map');colorbar;axis image;

% plot output data
subplot(3,4,5);imagesc(abs(outParams.water));caxis([0 1]);title('WATER image');colorbar;axis image;
subplot(3,4,6);imagesc(abs(outParams.fat));caxis([0 1]);title('FAT image');colorbar;axis image;
subplot(3,4,7);imagesc(abs(outParams.r2starmap));title('R2* map');caxis(r2star_caxisval); colorbar;axis image;
subplot(3,4,8);imagesc(angle(outParams.phasemap));title('Phase map'); phase_caxisval = caxis; colorbar;axis image;
subplot(3,4,12);imagesc(angle(exp(1i*(2*pi*trueParams.fieldmap.*(imDataParams.TE(2)-imDataParams.TE(1)))))); caxis(phase_caxisval);title('True phase map');colorbar;axis image;

subplot(3,4,9);imagesc(abs(abs(outParams.water)-abs(trueParams.species(1,1).amps)));caxis([0 0.3]);title('WATER error');colorbar;axis image;
subplot(3,4,10);imagesc(abs(abs(outParams.fat)-abs(trueParams.species(1,2).amps)));caxis([0 0.3]);title('FAT error');colorbar;axis image;
subplot(3,4,11);imagesc(abs(outParams.fiterror));caxis([0 0.3]);title('Fitting Error');colorbar;axis image;
clear r2star_caxisval phase_caxisval;