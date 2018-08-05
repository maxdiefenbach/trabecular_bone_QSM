% Project: Concentric Rings, v4
% Recon: ccr_savecine.m
%        save image matrix [N N Nech] as an .avi movie
%
%
% 2011/04/19
% Holden H. Wu

function mov_fps = ccr_savecine(imgt, movfname, WND_MAXALL, PLAY_MOV, cstr)
% imgt is [N N Nech]
% WND_MAXALL is 0~1.0
% cstr is 'None', etc. for movie2avi() codec

if( nargin<2 )
  movfname = 'ccr_cine.avi';
end  
if( nargin<3 )
  WND_MAXALL = 1.0;
end
if( nargin<4 )
  PLAY_MOV = 0;
end  
if( nargin<5 )
  cstr = 'None';
end  

% process input
imgt = abs( imgt );
Nech = size(imgt, 3);
% set saturation level
maxval = WND_MAXALL*max(imgt(:));
imgt( find(imgt>maxval) ) = maxval;

% convert to truecolor RGB
imgt = ( imgt - min(imgt(:)) )/( max(imgt(:)) - min(imgt(:)) );
imgtRGB = zeros(size(imgt,1),size(imgt,2),3,size(imgt,3));
imgtRGB(:,:,1,:) = imgt; 
imgtRGB(:,:,2,:) = imgt;
imgtRGB(:,:,3,:) = imgt;
% make movie from images
mov = immovie( imgtRGB );    
% play movie
if( PLAY_MOV )
  % position: [left, bottom, width, height]
  scrsz = get(0,'ScreenSize');
  %hf = figure('Position', [1, 1+scrsz(4)-size(imgt,1), size(imgt,2), size(imgt,1)]);
  hf = figure(1001);
  set(hf, 'Position', [1, 1+scrsz(4)-size(imgt,1), size(imgt,2), size(imgt,1)]);
  movie(hf, mov, round(Nech/2));
end

% No codecs seem to work on Win7 x64
%movfname = sprintf('%s/%s_cine_FW.avi', SAVE_DIR,PFile);
mov_fps = round(Nech/2);
movie2avi( mov, movfname, 'compression',cstr, 'fps',mov_fps );
    
return;