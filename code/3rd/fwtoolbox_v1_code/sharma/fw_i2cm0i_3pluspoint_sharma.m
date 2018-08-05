function [outParams] = fw_i2cm0i_3pluspoint_sharma(data,algoParams)

% Author:  Samir Sharma
% Created: December 2011
% This function controls the iterative decomposition

% INPUTS
%   data        the .mat file contents
%   algoParams  algorithm parameters

% OUTPUTS
%   outParams   contains the water, fat, and field map estimates    


%* BEGIN: Add path *%
addpath(strcat(pwd,'/utils'));
%* END: Add path *%

algoParams.species(1).ppm = algoParams.species(1).frequency; 
algoParams.species(2).ppm = algoParams.species(2).frequency; 


%* BEGIN: Check for errors, initialize some values *
data = error_check(data,algoParams);
if data.PrecessionIsClockwise <= 0 % DH*: made condition "<=0" because some datasets have  = -1 
    data.images = conj(data.images);
    data.PrecessionIsClockwise = 1;
end
[Nx,Ny,Nz,Ncoils,~] = size(data.images);
TE = data.TE;
%* END: Check for errors, initialize some values *


%* BEGIN: Calculate chemical shift encoding matrix, and coil combine *
A = calculate_chemical_shift_encoding_matrix(algoParams,data);

if (Ncoils > 1)
    image_comb = coilCombine(squeeze(data.images));
else
    image_comb = squeeze(data.images);
end
image_comb = permute(image_comb,[1 2 4 3]); % Nx x Ny x Nz x NTE
%* END: Calculate chemical shift encoding matrix, and coil combine *


%* BEGIN: Initialization *
field_map = zeros(Nx,Ny);
normalization = max(abs(image_comb(:))); 
image_comb = (10^3)*image_comb/normalization;
current_support = [1 1];
XFM = Bspline(current_support,[Nx Ny]);
%* END: Initialization *


%* BEGIN: Reconstruction *
counter = 1;
while (1)
    counter2 = 1;
    while (1)
        fprintf('Current iteration number is %d\n',counter2);
        %** BEGIN: Water-fat estimate **
        fprintf('Estimating water and fat images using current field map estimate... ');
        [wf_est,residual] = estimate_water_fat(image_comb,A,field_map,TE);
        fprintf('done\n');
        %** END: Water-fat estimate **

        %** BEGIN: Field map update **
        fprintf('Updating field map estimate with B-spline support size %dx%d... ',current_support(1),current_support(2));
        error_est = update_field_map(residual,A,field_map,TE,wf_est,XFM);
        fprintf('done\n\n');
        field_map_update = real(error_est(:,:,:,1));
        field_map = field_map + field_map_update;
        field_map_record(:,:,:,counter) = field_map; %#ok<AGROW>
        support_size_record(counter) = current_support(1); %#ok<AGROW>
        counter = counter + 1; counter2 = counter2 + 1;
% $$$         figure(1); subplot(1,3,1); imshow(field_map,[-400 400]); drawnow;
% $$$         subplot(1,3,2); imshow(abs(wf_est(:,:,:,1)),[]); drawnow;
% $$$         subplot(1,3,3); imshow(abs(wf_est(:,:,:,2)),[]); drawnow;
        if (counter2 > algoParams.MaxIter) || max(abs(field_map_update(:))) < 1
            break;
        end
        %** END: Field map update **
    end
    
    %** BEGIN: Update scale **
    if (current_support == algoParams.min_win_size)
        break;
    end
    if current_support == 1
        current_support = [Nx Ny];
    else
        current_support = round(algoParams.stepsize*current_support);
    end
    if sum(current_support<algoParams.min_win_size) >= 1
        current_support = [algoParams.min_win_size algoParams.min_win_size];
    end
    XFM = Bspline(current_support,[Nx Ny]);
    %** END: Update scale **
end
%* END: Reconstruction *



outParams.water = wf_est(:,:,:,1)*normalization/(10^3);
outParams.fat = wf_est(:,:,:,2)*normalization/(10^3);
outParams.fieldmap = field_map;

%% DH* Added this for improved interfacing with other functions and ability to extend to more species
outParams.species = algoParams.species;
outParams.species(1).amps = outParams.water;
outParams.species(2).amps = outParams.fat;
outParams.r2starmap = zeros(size(field_map));
