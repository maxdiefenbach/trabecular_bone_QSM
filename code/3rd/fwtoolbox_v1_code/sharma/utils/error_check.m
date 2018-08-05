function [data] = error_check(data,algoParams)

[Nx,Ny,Nz,Nc,NTE] = size(data.images); % DH* Let's get the number of coils too

% check data parameters
if NTE < 3
    error('ERROR: At least three echoes are required');
end
if Nz > 1
    disp('Only 2D datasets are permitted... reconstructing only the middle slice');
    data.images = data.images(:,:,floor((Nz+1)/2),:,:);
end

% DH* If more than one channel, coil combine
if Nc > 1
  disp('Multi-coil data: coil-combining');
  data.images = coilCombine(data.images);
end
  
if numel(algoParams.species) ~= 2
    error('ERROR: Sorry, this code only support two species');
end

% check algorithm specific parameters
if ~exist('algoParams.stepsize')
    algoParams.stepsize = 0.75;     % support size scaling 
end
if ~exist('algoParams.min_win_size')
    algoParams.min_win_size = 16;   % minimum 1D support size for B-spline
end
if ~exist('algoParams.MaxIter')
    algoParams.MaxIter = 10;        % maximum iterations per scale
end
pause(1); % for user to see messages