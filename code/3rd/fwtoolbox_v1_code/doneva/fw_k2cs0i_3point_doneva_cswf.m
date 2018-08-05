%% Function name: fw_k2cs0i_3point_doneva_cswf
%%
%% Description: Integrated Compressed Sensing and water fat separation
%%
%% Doneva M, Börnert P, Eggers H, Mertins A, Pauly J, Lustig M
%% Compressed sensing for chemical shift-based water-fat separation.
%% Magn Reson Med. 2010 Sep;64(6):1749-1759
%% 
%% Properties:
%%   - k-space
%%   - 2 species (water-fat)
%%   - Complex-fitting
%%   - Single/Multi peak fat  
%%   - Independent water/fat phase
%%   - Requires 3 echoes at arbitrary echo times (Field map initialization
%%   currently implemented for 3 echoes with constant dTE)
%%
%% Input: structures kDataParams and algoParams
%%   - kDataParams.kSpaceValues: acquired k-space data, array of size[nx,ny,1,ncoils,nTE]
%%   - kDataParams.kSpaceLocations:
%%   - kDataParams.TE: echo times (in seconds)
%%   - kDataParams.FieldStrength: (in Tesla)
%%
%%   - algoParams.species(ii).name = name of species ii (string)
%%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
%%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
%%   Example
%%      - algoParams.species(1).name = 'water' % Water
%%      - algoParams.species(1).frequency = [0] 
%%      - algoParams.species(1).relAmps = [1]   
%%      - algoParams.species(2).name = 'fat' % Fat
%%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
%%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
%% 
%%   - algoParams.Thresh = signal threshold     
%%   - algoParams.largeFM = (0,1) flag enabling few 1/dTE periods for the  field map estimation
%%   Optional parameters:
%%   - algoParams.num_outer_iter_wfcs = number of outer iterations for wf compressed sensing
%%   - algoParams.num_inner_iter_wfcs = number of inner iterations for wf compressed sensing
%%   - algoParams.lambda_tv_cs         = TV regularization parameter for cs
%%   - algoParams.lambda_wavelet_cs    = wavelets regularization parameter
%%   - algoParams.lambda_tv_wfcs       = TV regularization parameter for wfcs
%%   - algoParams.lambda_wavelet_wfcs  = wavelets regularization
%%
%% Output: structure outParams
%%   - outParams.species(ii).name: name of the species (taken from algoParams)
%%   - outParams.species(ii).amps: estimated water/fat images, size [nx,ny,ncoils] 
%%   - outParams.fieldmap: field map (in Hz, size [nx,ny])
%%
%%
%% Author: Mariya Doneva
%% Date created: 11 October 2011
%% Date last modified: 12 January 2012

function outParams = fw_k2cs0i_3point_doneva_cswf( kDataParams, algoParams )

%% Initialize the CS water-fat reconstruction 
%% Perform CS on each echo image
%% Estimate initial field map 

DEBUG = 0;

% Check validity of params, and set default algorithm parameters if not provided
[validParams,algoParams] = checkParamsAndSetDefaults_cswf( kDataParams,algoParams );
if validParams==0
  disp('Exiting -- data not processed');
  outParams = [];
  return;
end


%% Reorder k-space data and set sampling mask

sx = kDataParams.FOV(1);
sy = kDataParams.FOV(2);

t     = kDataParams.kSpaceTimes;
TE    = unique(t);
dTE   = TE(2)-TE(1);
N     = length(TE);  
kdata = zeros(kDataParams.FOV(1),kDataParams.FOV(2),N);

for kt = 1:N
    locations_idx = find(kDataParams.kSpaceTimes==TE(kt));
    idx = sub2ind([sx,sy],kDataParams.kSpaceLocations(locations_idx,1)+ sx/2 +1, kDataParams.kSpaceLocations(locations_idx,2)+ sy/2 +1 );
    temp = kdata(:,:,kt);
    temp(idx) = kDataParams.kSpaceValues(locations_idx);
    kdata(:,:,kt) = temp;    
end

mask = abs(kdata)>0;



%% Initialize CS parameters
params = initializecg_cs(sx, sy);

if (isfield(algoParams,'lambda_wavelet_cs'))
  params.WTWeight = algoParams.lambda_wavelet_cs;
end
if (isfield(algoParams,'lambda_tv_cs'))
  params.TVWeight = algoParams.lambda_tv_cs;
end


%% Perform CS reconstruction for each echo
disp ('Initial CS iterations for echo images');
image = zeros(sx,sy,N);

for kt = 1:N
    
    
    ii = ifft2c(kdata(:,:,kt));
    factor = max(abs(ii(:)));
    image(:,:,kt) = ii/factor ;
    params.data = kdata(:,:,kt)/factor;
    params.image = image(:,:,kt);
    params.FT = p2DFT(mask(:,:,kt), size(mask(:,:,kt)), 1, 2);
    
    for ki = 1:params.n_outer
        params.image = image(:,:,kt);
        image(:,:,kt) = nonlin_cg(params);
    end
   
    image(:,:,kt) = image(:,:,kt)*factor;
    
end

if (DEBUG)
figure(4); 
subplot(1,3,1);
imshow(abs(image(:,:,1)),[]); title('TE1');
subplot(1,3,2);
imshow(abs(image(:,:,2)),[]); title('TE2');
subplot(1,3,3);
imshow(abs(image(:,:,1)),[]); title('TE3');
drawnow;
end

%% Estimate initial field map
% This is a heutistic method which does not guarantee the correct field map estimation 
% It may be replaced by other method for field map estimation 
%% Compute frequencies in Hz
gyro = 42.58;
fieldStrength = kDataParams.FieldStrength;
f_wf = [];
rel_amp = [];
for i = 1:length(algoParams.species)
f_wf = [f_wf; gyro*fieldStrength*algoParams.species(i).frequency(:)];
rel_amp = [rel_amp; algoParams.species(i).relAmps(:)];
end

disp ('Compute maps');
% Estimate all possible field map values
maps = locatemins_outer1(image,algoParams,dTE,f_wf);
% Estimate map of large image amplitudes
amp_map = sum(abs(image),3)>algoParams.Thresh*max(abs(image(:)));
%%%%%%%%%%%%%%%%%% Free up some memory %%%%%%%
clear image; 
clear params;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%% Select initial field map 
disp ('Select map');
% Initialize field map search
m1 = (maps(:,:,1));
m2 = (maps(:,:,2));
ind = abs(m1)<abs(m2);
fmap = m1.*ind + m2.*(1-ind);
fmap = fmap.*amp_map;
fmap = ones(size(fmap)).*amp_map*median(fmap(fmap~=0));
clear m1;
clear m2;


% Augment search region to 3 periods if large field map option is selected
if (algoParams.largeFM ==1)
   maps1(:,:,1:2) =  maps(:,:,1:2) - 1/dTE;
   maps1(:,:,3:4) =  maps(:,:,1:2);
   maps1(:,:,5:6) =  maps(:,:,1:2) + 1/dTE;  
else
   maps1 = maps(:,:,1:2);
end



% Perform search on subsampled region
fmap1 = fmap(1:4:end,1:4:end);
% There is no great idea behind the choice of these points. The goal is to
% start from few different points,there might be a more intelligent way to
% set these
seeds = [size(fmap1,1)*7/8  size(fmap1,1)/8 size(fmap1,1)/7 size(fmap1,1)/8 size(fmap1,1)*7/8 size(fmap1,1)*7/8; size(fmap1,2)/2   size(fmap1,2)/8  size(fmap1,2)/7 size(fmap1,2)*7/8   size(fmap1,2)/8 size(fmap1,2)*7/8 ];

dist = 1;
i = 1;
time1 = tic;
while((dist > 0)&&(i<6))
y_seed = round(seeds(2,i));
x_seed = round(seeds(1,i));
fmap_new1 = rm2d2( maps1(1:4:end,1:4:end,:), fmap1, [x_seed,y_seed], fmap1(x_seed,y_seed),amp_map(1:4:end,1:4:end), 3, 1, DEBUG);
dist = (fmap_new1(:) - fmap1(:))'*(fmap_new1(:) - fmap1(:));
fmap1 = fmap_new1;
i = i+1;
end

fmap1 = imresize(fmap1,2).*amp_map(1:2:end,1:2:end);

% selection on full resolution map
x_seed = round(size(fmap1,1)/2);
y_seed = round(size(fmap1,2)/2);
fmap1 = rm2d2( maps1(1:2:end,1:2:end,:), fmap1, [x_seed,y_seed], fmap(x_seed,y_seed),amp_map(1:2:end,1:2:end), 3, 0, DEBUG);

fmap1 = imresize(fmap1,2).*amp_map;

% selection on full resolution map
x_seed = round(size(fmap1,1)/2);
y_seed = round(size(fmap1,2)/2);
fmap = rm2d2( maps1, fmap1, [x_seed,y_seed], fmap(x_seed,y_seed),amp_map, 3, 0, DEBUG);

time1 = toc(time1);
disp('Time for field map estimation'); disp(time1);

%% Actual CS-WF code starts here
% Initialize CSWF reconstruction

% initialization x-vector 
% x contains the water image, fat image and field map
% initialize water and fat images with 0 

disp ('Perform CS-WF iterations');
x0(:,:,1:2) = zeros(sx,sy,2);  	                   
x0(:,:,3)   = fmap;

cg_param = initializecg_wf([sx,sy], f_wf , rel_amp,TE, mask);
% set user input parameters

if (isfield(algoParams,'num_inner_iter_wfcs'))
  cg_param.num_inner_iter_wfcs = algoParams.num_inner_iter_wfcs;
end
 
if (isfield(algoParams,'num_outer_iter_wfcs'))
  cg_param.num_outer_iter_wfcs = algoParams.num_outer_iter_wfcs;
end

if (isfield(algoParams,'lambda_wavelet_wfcs'))
  cg_param.WTWeight = algoParams.lambda_wavelet_wfcs;
end

if (isfield(algoParams,'lambda_tv_wfcs'))
  cg_param.FD1Weight = algoParams.lambda_tv_wfcs;
end

if(DEBUG)
    cg_param.display = 1;
end



time1 = tic;
[w,f,fm] = wfcs(kdata, x0, cg_param);
time1 = toc(time1);
disp('Time for CS-WF iterations'); disp(time1);



% Put results in outParams structure
try
  outParams.species(1).name = algoParams.species(1).name;
  outParams.species(2).name = algoParams.species(2).name;
catch
  outParams.species(1).name = 'water';
  outParams.species(2).name = 'fat';
end  
  
outParams.species(1).amps = w;
outParams.species(2).amps = f;
outParams.fieldmap = real(fm);
outParams.r2starmap = imag(fm)*2*pi;




