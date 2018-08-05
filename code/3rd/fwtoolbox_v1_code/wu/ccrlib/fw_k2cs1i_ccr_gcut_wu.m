function [outParams] = fw_k2cs1i_ccr_gcut_wu( kDataParams, gcutParams, algoParams )
% Project: Concentric Rings Fat/Water Imaging
% Filename: fw_k2cs1i_ccr_gcut_wu.m
%   Concentric rings fat/water separation using
%   graph-cut and frequency-demodulated least-squares estimation.
%   (cf. Wu et al., MRM 2009;61:639-649)
%   (cf. Hernando et al., MRM 2010;63:79-90)
% Signal model (k2cs1i):
%   k-space, 2D, complex, single peak fat, single R2*, indep initial phase
%   fat has negative frequency -3.5 ppm
%   field map effect is exp(1i*2*pi*fMap*TE)
%   (NOTE1: this is a hybrid algorithm where k-space data is demodulated
%    at the frequencies of interest prior to image-space processing)
%   (NOTE2: 3D stack of rings is processed slice by slice)
% Echo time requirements:
%   3+ revolutions (TEs)
%   only tested for uniform TE separation, but should also work for
%   non-uniform TE separation
%
% function [outParams] = fw_k2cs1i_ccr_gcut_wu( kDataParams, gcutParams, algoParams )
%     =INPUT=
%     kDataParams (struct) --
%       B0      -- field strength, in Tesla
%       FOV     -- image field of view, in cm
%       Ncr     -- no. of concentric rings
%       imgN    -- nominal image matrix size (square), ==2*Ncr
%       Tsamp   -- data sampling rate, in micro sec
%       ginfo (struct) -- details of the gradient waveforms
%       TE0     -- time of first readout sample after RF, in sec
%       ksp     -- Nrd x Ncr complex, k-space sampling locations
%       dcf     -- Nrd x Ncr, k-space density compensation values
%       kdata   -- Nrd x Ncr x Nph x Nsl x Nch complex, k-space data
%       Nrd     -- no. of readout points along each ring
%       Nph     -- no. of temporal (e.g., cardiac) phases
%       Nsl     -- no. of slices
%       Nch     -- no. of channels
%     gcutParams (struct) --
%       see fw_i2cm1i_3pluspoint_hernando_graphcut() for more info
%     algoParams (struct) --
%       padN    -- zero-interpolate images from imgNximgN to padNxpadN
%       TEs     -- N_TE x 1, echo time for each source image, in sec
%       f0      -- off-res freq for target species (0:water, df_FW:fat), in Hz
%       df_FW   -- off-res freq for fat, in Hz
%       cFOV    -- 1:impose circular FOV, 0:off
%       rMask   -- padN x padN, mask used for cFOV
%       imrot   -- rotate images by imrot*90 deg
%
%     =OUTPUT=
%     outParams (struct) --
%       species(ii).name -- name of the species, e.g. 'water','fat'
%       species(ii).amps -- padN x padN x Nph x Nsl, estimation results (coils combined)
%       fMap    -- padN x padN x Nph x Nsl, field map
%       r2star  -- padN x padN x Nph x Nsl, r2star map
%
% Modified:
%   2012/02/15 Hernando's function now combines coils beforehand
%                therefore, adjust code accordingly here
%   2012/01/16 now using ccr_k2im()
%   2012/01/12 disable R2* estimation in default settings
%              outParams now using 'species' convention
%   2011/12/28 using Diego Hernando's graph cut function for water recon
% Last Modified: 2012/01/16
% Created:       2011/12/28
% Holden H. Wu
% holdenhwu@stanford.edu

% init output
outParams = [];
outParams.species(1).name = 'water';
outParams.species(2).name = 'fat';

% Extract data sizes
Ncr  = kDataParams.Ncr;  % # concentric rings
Nph  = kDataParams.Nph;  % # temporal phases
Nsl  = kDataParams.Nsl;  % # slices
Nch  = kDataParams.Nch;  % # channels
imgN = kDataParams.imgN; % nominal size of image (square)

% Set default values if gcutParams doesn't exist
if( nargin<2 )
% ------------------------------------------------------------------------
fprintf(1, 'Warning: gcutParams not specified, using default values.\n');

% General parameters
gcutParams.species(1).name = 'water';
gcutParams.species(1).frequency = 0;
gcutParams.species(1).relAmps = 1;
gcutParams.species(2).name = 'fat';
%gcutParams.species(2).frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];
%gcutParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
% Diego's function supports multi-peak fat, but use single peak for now
gcutParams.species(2).frequency = [-3.50]; % -3.5 ppm
gcutParams.species(2).relAmps = [1.00]; % single peak for now

% Algorithm-specific parameters
gcutParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
%gcutParams.range_r2star = [0 100]; % Range of R2* values
%gcutParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
% Hernando's function supports R2*, but turn off for now
gcutParams.range_r2star = [0 0]; 
gcutParams.NUM_R2STARS = 1;
gcutParams.range_fm = [-400 400]; % Range of field map values
gcutParams.NUM_FMS = 301; % Number of field map values to discretize
gcutParams.NUM_ITERS = 40; % Number of graph cut iterations
gcutParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
gcutParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
gcutParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
gcutParams.lambda = 0.05; % Regularization parameter
gcutParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
gcutParams.TRY_PERIODIC_RESIDUAL = 0;
% ------------------------------------------------------------------------
end

% Set default values if algoParams doesn't exist
if( nargin<3 )
% ------------------------------------------------------------------------
fprintf(1, 'Warning: algoParams not specified, using default values.\n');

zpadX = 1; % zero-interpolate the images zpadX-fold
algoParams.padN = zpadX * imgN;

GAM = 42.576*1e6; % Hz/Tesla, for protons
algoParams.df_FW  = -3.5*1e-6 * GAM * kDataParams.B0; % Hz
algoParams.f0     = 0;
algoParams.cFOV   = 1;

% impose circular FOV
xaxis = -algoParams.padN/2:(algoParams.padN/2-1);
[X, Y] = meshgrid( xaxis, xaxis );
R = sqrt( X.^2 + Y.^2 );
algoParams.rMask = zeros(algoParams.padN, algoParams.padN);
algoParams.rMask( R < algoParams.padN/2 ) = 1;

algoParams.TEs    = [];
algoParams.imrot  = 0;
% ------------------------------------------------------------------------
end    

padN = algoParams.padN;


% Flags
% ========================================================================
% For Fat/Water separation
RECONFAT  = 1; % whether to recon fat on fat-res
% PLOT
PLOT_SRC0 = 0; % plot source images for all CH (ph1, rev1)
PLOT_SRC  = 0; % plot source images
PLOT_DBG  = 0; % plot intermediate results for debugging
%PLOT_FMAP = 0; % plot final fMap from iterative calculations
% ========================================================================


% Recon all source images at water and fat center frequency
% ========================================================================
NTE = kDataParams.ginfo.REVLN(1);
if( NTE<3 )
  fprintf(1, 'Error: NTE<3, please load a compatible dataset.\n');
  return;
end

imsrcW = ccr_k2im(kDataParams, algoParams, 0);
algoParams.TEs = imsrcW.TEs;
fprintf(1,'Reconstructed all source images (water demod)\n');

if( RECONFAT )
  imsrcF = ccr_k2im(kDataParams, algoParams, algoParams.df_FW);
  fprintf(1,'Reconstructed all source images (fat demod)\n');
end

if( PLOT_SRC0 )
  ph=1; sl=1;
  imgplot = abs( reshape(imsrcW.images(:,:,1,ph,sl,:), imgN,imgN,1,Nch) );
  figure(201); 
  montage(imgplot, [min(imgplot(:)) max(imgplot(:))]);
  title(sprintf('imsrcW (Nch=%d), slice %d/%d, phase %d/%d, revln %d/%d', Nch,sl,Nsl,ph,Nph,1,NTE));

  %fprintf(1,'Please specify appropriate CHset.\n'); pause;
end    
% ========================================================================


% Fat/water separation using Hernando's graph cut algorithm for water recon
%   and modified IDEAL for fat recon (w/known field map from water)
% init image-space structures for graph cut
imDataParams.TE = imsrcW.TEs;
imDataParams.images = zeros(padN, padN, 1, Nch, NTE);
imDataParams.PrecessionIsClockwise = -1;
imDataParams.FieldStrength = kDataParams.B0;
% init output structures
outParams.species(1).amps = zeros( padN, padN, Nph, Nsl ); % water
outParams.species(2).amps = zeros( padN, padN, Nph, Nsl ); % fat
outParams.fMap  = zeros( padN, padN, Nph, Nsl );
fatAllCh = zeros( padN, padN, Nph, Nsl, Nch );
% SLICE loop
for sl=1:Nsl,
% ************************************************************************
% PHASE loop
for ph=1:Nph,
% ========================================================================

% WATER RECON
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
% Hernando's function supports multi-peak fat, but use single peak for now
% all channels are processed in the graph cut function

if( kDataParams.ginfo.REVLN(1)<4 && gcutParams.NUM_R2STARS>1 )
  fprintf(1, 'Warning: Estimating R2* with <4 source images.\n');
end  

% imsrcW [padN, padN, NTE, Nph, Nsl, Nch]
% --> images [padN, padN, 1, Nch, NTE]
imDataParams.images = permute( imsrcW.images(:,:,:,ph,sl,:), [1 2 4 6 3 5] );
imDataParams.images = reshape( imDataParams.images, [padN,padN,1,Nch,NTE] );

gcutW = fw_i2cm1i_3pluspoint_hernando_graphcut( imDataParams, gcutParams );

outParams.species(1).amps( :, :, ph,sl ) = gcutW.species(1).amps;
outParams.fMap( :, :, ph,sl )    = -gcutW.fieldmap;
outParams.r2star( :, :, ph,sl )  = gcutW.r2starmap;
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-


% FAT RECON
if( RECONFAT )
% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-

% COIL loop
for ch=1:Nch, 
% ========================================================================
% Source Images 
% ------------------------------------------------------------------------
img_mat = reshape( imsrcF.images(:,:,:,ph,sl,ch), padN^2,NTE ).';

if( PLOT_SRC )
    w_max = max( abs(img_mat(:)) );
    w_min = min( abs(img_mat(:)) );
    
    for im=1:NTE,
    figure; imagesc( abs( reshape(img_mat(im,:), padN,padN) ), [w_min, w_max] ); 
    colormap gray; truesize; axis image; title(sprintf('Src Image %d (Fat)',im));
    end
end
% ------------------------------------------------------------------------

% iterative least-squares estimation (modified IDEAL)
% use known field map from water recon, don't need any iterations
algoParams.f0      = algoParams.df_FW;
algoParams.fMap0   = -gcutW.fieldmap;
algoParams.mlevel  = 0.1;
algoParams.maxiter = 0;
outParamsF = fdIDEAL(img_mat, algoParams);

% Plot
if( PLOT_DBG )
%figure; imagesc(abs(outParamsF.img_W)); colormap gray; truesize; axis image; title('Water Image (F)');
figure; imagesc(abs(outParamsF.img_F), [w_minW, w_maxW]); colormap gray; truesize; axis image; title('Fat Image WL (F)');
%figure; imagesc(abs(outParamsF.img_R)); colormap gray; truesize; axis image; title('Residual (F)');
figure; imagesc(abs(outParamsF.img_F) + abs(outParamsW.img_W)); colormap gray; truesize; axis image; title('|Fat|+|Water|');
end

% keep track of results
%outParams.species(2).amps( :, :, ph,sl,ch ) = outParamsF.img_F;
fatAllCh( :, :, ph,sl,ch ) = outParamsF.img_F;
% ========================================================================
end % COIL loop

% combine coils using root of sum of squares
outParams.species(2).amps( :, :, ph,sl ) = ...
    sqrt( sum( abs(fatAllCh(:, :, ph,sl,:)).^2, 5 ) );
% NOTE: could also use Hernando's coilCombine() for source images beforehand ...

% =-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-=-
end % end of RECONFAT

fprintf(1,'Fat/water separation: sl %d/%d, ph %d/%d, Nch=%d\n', sl,Nsl,ph,Nph,Nch);
% ========================================================================
end % PHASE loop
% ************************************************************************
end % SLICE loop


return;
