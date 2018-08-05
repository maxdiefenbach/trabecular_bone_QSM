function [img, gdat] = grecon2d(ktraj, kdata, dcf, imgN, a, kwidth, padN, imrot)

% ========================================================================
% [img, gdat] = grecon2d(ktraj, kdata, dcf, imgN, a, kwidth, padN, imrot)
% Uses mex function gr2dKB_mex()
%     =INPUT= 
%     ktraj  -- Npts x 1 complex, k-trajectory, scaled -0.5 to 0.5
%     kdata  -- Npts x 1 complex, k-space data
%     dcf    -- Npts x 1, k-space weighting, 
%     imgN   -- 1x1, image size (grid will be a*imgN x a*imgN) 
%     a      -- 1x1, oversampling ratio 
%     kwidth -- 1x1, kernel width
%     =optional=
%     padN   -- 1x1, zero pad output image matrix to padN x padN
%     imrot  -- 1x1, number of times to rotate recon img by 90deg CCW
%
%     =OUTPUT= 
%     img    -- (imgN x imgN) or (padN x padN), reconstructed image  
%     gdat   -- a*(imgN x imgN) or a*(padN x padN), gridded k-space data 
% ========================================================================
% Modified:
%           2009/11/08 now takes padN instead of pad
%                      assimilates function of grecon2dPad.m
% 2007, Nov 
% Holden H. Wu
% ========================================================================

if( nargin<5 )
    a = 2; % 2X grid
end
if( nargin<6 )
    kwidth = 5;
end
if( nargin<7 )
    padN = imgN;
end    
if( nargin<8 )
    imrot = -1;
    %imrot = 0;
end    

if( size(ktraj,2)==2 )
    ktraj = ktraj(:,1) + 1i*ktraj(:,2);
end
% NOTE: Should do some more checking.


% Pre-compute KB kernel, cf. Beatty et al., IEEE TMI 2005; 24:799-808
b = pi*sqrt( (kwidth/a)^2*(a-0.5)^2 - 0.8  ); % b=11.441 for a=2&kwidth=5
G = 1;
% sample the kernel and store as a lookup table
max_eps1 = 1e-4;
S = sqrt(0.37/max_eps1)/a; % linear interp
S = 10*ceil( S/10 ); % S samples per grid unit
kx_samp = 0:1/(G*S):kwidth/(2*G); % only need one side
C_samp = (G/kwidth) * besseli( 0, b*sqrt(1-(2*G*kx_samp/kwidth).^2) );

% 2D Gridding
% convert to single column
ktraj = ktraj(:);
kdata = kdata(:) + 1i*eps; % make sure it's complex
dcf   = dcf(:);
% initialize k-space grid 
grsize = ceil(a*imgN); %grcenter = (a*imgN/2+1);
gdat   = zeros(grsize, grsize);
% call mex function
gdat = gr2dKB_mex(ktraj, kdata, dcf, grsize, kwidth/2, C_samp);

% Zero-pad to larger image matrix
if( padN>imgN )
    grsizePad = ceil(a*padN);
    gdat2 = zeros( grsizePad, grsizePad );
    x0 = floor( (grsizePad-grsize)/2 );
    gdat2(x0+[1:grsize], x0+[1:grsize]) = gdat;
    imgN = padN;
    grsize = grsizePad;
    gdat = gdat2;
    clear gdat2;
end

% 2D-IFFT
img = fftshift( ifft2( ifftshift(gdat) ) );
img = rot90( img, imrot );

% Deapodization
% using analytical expression of FT{kernel}, cf. Beatty et al.
argvec = (pi*kwidth*[-grsize/2:grsize/2-1]/grsize).^2 - b^2;
c_1D = sin(sqrt( argvec )) ./ sqrt( argvec );
% create deap matrix
c_2D = c_1D' * c_1D;
c_2D = c_2D + 10000; % baseline offset added to deap function
img = img ./ c_2D;

% Extract FOV
x0  = floor( (a-1)*imgN/2 );
img = img( x0+[1:imgN], x0+[1:imgN] );

return;
