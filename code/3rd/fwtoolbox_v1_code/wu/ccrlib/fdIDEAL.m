function [outParams] = fdIDEAL(img_mat, algoParams)
% Project: Concentric Rings Fat/Water Imaging
% Filename: fdIDEAL.m
%   Frequency-demodulated iterative least-squares fat/water estimation.
%   (cf. Wu et al., MRM 2009;61:639-649)
%   (modified from IDEAL, cf. Reeder et al., MRM 2004;51:35-45)
% Signal Model:
%   fat has negative frequency -3.5 ppm
%   field map effect is exp(1i*2*pi*fMap*TE)
%
% function outParams = fdIDEAL(img_mat, algoParams)
%     =INPUT=
%     img_mat -- N_TE x (imgN^2) complex, source images demodulated at f0
%     algoParams (struct) --
%       TEs     -- N_TE x 1, echo time for each source image, in sec
%       f0      -- off-res freq for target species (0:water, df_FW:fat), in Hz
%       df_FW   -- off-res freq for fat, in Hz
%       fMap0   -- imgN x imgN, initial fMap, in Hz
%       mlevel  -- 0~1, threshold value for magnitude mask
%       maxiter -- max number of iterations for IDEAL
%       cFOV    -- 1:impose circular FOV, 0:off
%       rMask   -- imgN x imgN, mask used for cFOV
%       convk   -- kernel used for low-pass filtering of field map
%
%     =OUTPUT=
%     outParams (struct) --
%       img_W   -- imgN x imgN complex, reconstructed water image
%       img_F   -- imgN x imgN complex, reconstructed fat image
%       img_R   -- imgN x imgN, residual error image
%       fMap    -- imgN x imgN, calculated field map, in Hz
%       mask    -- imgN x imgN, magnitude mask used for fMap processing
%
% Last Modified: 2011/12/26
% Created:       2009, Feb
% Holden H. Wu
% holdenhwu@stanford.edu

% Flags
PLOT_DBG = 0;

% Extract data sizes
N_TE = numel( algoParams.TEs );
imgN = sqrt(size(img_mat, 2));
if( size(img_mat,1) ~= N_TE )
    fprintf(1, 'TEvec and img_mat do not match!\n');
    return;
end

if( isempty(algoParams.fMap0) )
    % initial est. of field map
    % use all zeros 
    algoParams.fMap0 = zeros( imgN, imgN );    
end

% Initialize matrices
% ------------------------------------------------------------------------
img_matD = zeros( N_TE, imgN^2 );
outParams.fMap = algoParams.fMap0;
% demod field map
for im=1:N_TE,
    img_matD(im,:) = img_mat(im,:) .* [exp(-1i*2*pi*outParams.fMap(:)*algoParams.TEs(im)).'];
    %img_matD(im,:) = img_mat(im,:) .* [exp(1i*2*pi*outParams.fMap(:)*algoParams.TEs(im)).'];
end
% Setup input S
S = zeros(N_TE*2, imgN^2);
% Real part
for im=1:N_TE,
    S(im, :) = real(img_matD(im,:));
end    
% Imag part
for im=1:N_TE,
    S(im+N_TE, :) = imag(img_matD(im,:));
end    

% Setup matrix A at f0 (data already demodulated at f0)
c_fww = cos(-2*pi*(0-algoParams.f0)*algoParams.TEs);
d_fww = sin(-2*pi*(0-algoParams.f0)*algoParams.TEs);
c_fwf = cos(-2*pi*(algoParams.df_FW-algoParams.f0)*algoParams.TEs);
d_fwf = sin(-2*pi*(algoParams.df_FW-algoParams.f0)*algoParams.TEs);

A = [c_fww.' -d_fww.' c_fwf.' -d_fwf.';
     d_fww.'  c_fww.' d_fwf.'  c_fwf.'];

% LS solution
FW = A \ S;     
% residual
R = S - A*FW; 

% water img
outParams.img_W = FW(1, :) + 1i*FW(2, :);
% fat img
outParams.img_F = FW(3, :) + 1i*FW(4, :);
% residual img
outParams.img_R = sum(R.*R, 1);

% reshape
outParams.img_W = reshape( outParams.img_W, imgN, imgN );
outParams.img_F = reshape( outParams.img_F, imgN, imgN );
outParams.img_R = reshape( outParams.img_R, imgN, imgN );

if( PLOT_DBG )
figure; imagesc(abs(outParams.img_W)); colormap gray; truesize; axis image; title('Water Img - 0th iter');
figure; imagesc(abs(outParams.img_F)); colormap gray; truesize; axis image; title('Fat Img - 0th iter');
figure; imagesc(abs(outParams.img_R)); colormap gray; truesize; axis image; title('Residual Img - 0th iter');
end

% Water + Fat
img_WF = abs(outParams.img_W) + abs(outParams.img_F);
% Weighting mask for the 6 source channels (N_TE x2)
outParams.mask = abs(img_WF) > algoParams.mlevel*max( abs( img_WF(:) ) );

%figure; imagesc(img_WF); colormap gray; truesize; axis image; title('|W|+|F| (W)');
%figure; imagesc(mask); colormap gray; truesize; axis image; title('mag mask');

clear img_WF;
% ------------------------------------------------------------------------


% IDEAL
B = zeros(6, 5);
B(:, 2:end) = A;
Y = zeros(5, imgN^2);
for iter=1:algoParams.maxiter,
% ------------------------------------------------------------------------
g_R = -2*pi*(algoParams.TEs'*ones(1,3).*[ones(3,1) d_fwf.'  c_fwf.']) *FW([2 3 4], :);
g_I =  2*pi*(algoParams.TEs'*ones(1,3).*[ones(3,1) c_fwf.' -d_fwf.']) *FW([1 3 4], :);
g = [g_R; g_I];

%tt0 = cputime;
% calculate dPSI for each voxel
for kk=1:imgN^2,
    %B = [g(:,kk) A];
    if( ~algoParams.cFOV || algoParams.rMask(kk) )
    B(:,1) = g(:,kk);
    Y(:,kk) = B \ R(:,kk);
    end
end    
%tt = cputime-tt0;

dfMap  = outParams.mask .* reshape( Y(1,:), imgN, imgN );
% remove salt-pepper
dfMap = medfilt2( dfMap );
% smooth
dfMap = filter2(algoParams.convk, dfMap);
if( PLOT_DBG )
    max(abs(dfMap(:)))
end

%img_dW = mask .* reshape( Y(2,:)+i*Y(3,:), imgN, imgN );
%img_dF = mask .* reshape( Y(4,:)+i*Y(5,:), imgN, imgN );

outParams.fMap = outParams.fMap + dfMap;
%outParams.fMap = outParams.fMap - dfMap;
%figure; imagesc(fMap); colormap gray; truesize; axis image; title('freq map');

% demod field map
for im=1:N_TE,
    img_matD(im,:) = img_mat(im,:) .* [exp(-1i*2*pi*outParams.fMap(:)*algoParams.TEs(im)).'];
    %img_matD(im,:) = img_mat(im,:) .* [exp(1i*2*pi*outParams.fMap(:)*algoParams.TEs(im)).'];
end
% Setup input S
S = zeros(N_TE*2, imgN^2);
% Real part
for im=1:N_TE,
    S(im, :) = real(img_matD(im,:));
end    
% Imag part
for im=1:N_TE,
    S(im+N_TE, :) = imag(img_matD(im,:));
end  

% LS solution
FW = A \ S;
% residual
R = S - A*FW;  
% ------------------------------------------------------------------------
end

% Water img
outParams.img_W = FW(1, :) + 1i*FW(2, :);
% Fat img
outParams.img_F = FW(3, :) + 1i*FW(4, :);
% Residual img
outParams.img_R = sqrt( sum(R.*R, 1) );
% reshape
outParams.img_W = reshape( outParams.img_W, imgN, imgN );
outParams.img_F = reshape( outParams.img_F, imgN, imgN );
outParams.img_R = reshape( outParams.img_R, imgN, imgN );

return;
