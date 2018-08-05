function [imDataParams] = ccr_k2im(kDataParams, algoParams, fd)
% Project: Concentric Rings Fat/Water Imaging
% Filename: ccr_k2im.m
%     Split concentric rings multi-rev k-space data into individual revs,
%     demodulate k-space data at desired 'fd',
%     and reconstruct each rev with 2D gridding.
%
% function [imDataParams] = ccr_k2im(kDataParams, algoParams, fd)
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
%     algoParams (struct) --
%       padN    -- zero-interpolate images from imgNximgN to padNxpadN
%       cFOV    -- 1:impose circular FOV, 0:off
%       rMask   -- padN x padN, mask used for cFOV
%       imrot   -- rotate images by imrot*90 deg
%     fd -- demodulation frequency
%
%     =OUTPUT=
%     imDataParams (struct) --
%       images  -- padN x padN x NTE x Nph x Nsl x Nch complex, recon results
%       TEs     -- N_TE x 1, effective TE for each source image
%
% Last Modified: 2012/01/16
% Created:       2012/01/16
% Holden H. Wu
% holdenhwu@stanford.edu

% Extract data sizes
Ncr  = kDataParams.Ncr;  % # concentric rings
Nph  = kDataParams.Nph;  % # temporal phases
Nsl  = kDataParams.Nsl;  % # slices
Nch  = kDataParams.Nch;  % # channels
imgN = kDataParams.imgN; % nominal size of image (square)
padN = algoParams.padN;  % zero-pad images from imgN^2 to padN^2

if( nargin<3 )
    fd = 0;
end    

% Split multi-rev acq into each individual rev
% ------------------------------------------------------------------------
NTE   = kDataParams.ginfo.REVLN(1); % should match numel(Eidx)
NpRev = (kDataParams.ginfo.lenS.lenRd/NTE); % should be an integer
Eidx  = 1 + NpRev*[0:NTE-1]; % starting sample of each rev
TEshift = kDataParams.TE0 + (Eidx-1)*kDataParams.ginfo.T*1e-6; % corresponding effective TEs
% Reshape
kDataParams.kdata = reshape(kDataParams.kdata, [NpRev, NTE, Ncr, Nph, Nsl, Nch]);
kDataParams.ksp   = reshape(kDataParams.ksp, [NpRev, NTE, Ncr]);
kDataParams.dcf   = reshape(kDataParams.dcf, [NpRev, NTE, Ncr]);

imDataParams.TEs = TEshift;
% ------------------------------------------------------------------------

% demodulate k-space data at fd
if( abs(fd)>0 )
% ------------------------------------------------------------------------
tvec = [0:kDataParams.ginfo.lenS.lenRd-1]*kDataParams.ginfo.Ta*1e-6;
demodM = ones(Ncr,1)*exp(1i*2*pi* fd *tvec);
demodM = reshape(demodM.', [NpRev, NTE, Ncr]);
for sl=1:Nsl,
  for ph=1:Nph,
    for ch=1:Nch,      
      kDataParams.kdata(:,:,:,ph,sl,ch) = kDataParams.kdata(:,:,:,ph,sl,ch) .* demodM;
    end
  end
end
% ------------------------------------------------------------------------
end

% k-space --> image space
% ------------------------------------------------------------------------
% init matrix
imDataParams.images = zeros(padN, padN, NTE, Nph, Nsl, Nch);
% Params for gridding
grdfactor = 1.5; % 2
kwidth    = 4;   % 1.5
for sl=1:Nsl,
  for ph=1:Nph,
    for ch=1:Nch,
      for im=1:NTE,
      % extract rev
      kdata = reshape( kDataParams.kdata( :,im,:,ph,sl,ch ), NpRev, Ncr );    
      % for gridding recon
      kspR = kDataParams.ksp(:,im,:); dcfR = kDataParams.dcf(:,im,:);    

      % gridding recon
      img_TE = grecon2d( kspR(:), kdata(:), dcfR(:), imgN, grdfactor, kwidth, padN );
      img_TE = rot90( img_TE, algoParams.imrot );
      if( algoParams.cFOV ) 
          img_TE = img_TE .* algoParams.rMask; 
      end
      imDataParams.images( :,:,im,ph,sl,ch ) = img_TE;      
      end
    end
  end
end
clear kdata;
% ------------------------------------------------------------------------

return;
