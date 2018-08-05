%
% [RgbImg]=ComplexToRgb(data, ShowLogMagnitude),
%

% REVISION HISTORY
% 2003.05.09 - Jeffrey Tsao
%     - Converted to callable function based on
%		"script - create k-space movie.m" in
% 2009.05.18 - Jeffrey Tsao
%     - To allow arbitrary colormaps
function [RgbImg]=ComplexToRgb(data, ShowLogMagnitude, cmap),
if nargin<1, help(mfilename); return; end
if nargin<2, ShowLogMagnitude=[]; end
if nargin<3, cmap=[]; end
if isempty(ShowLogMagnitude), ShowLogMagnitude=0; end

if isempty(cmap),
    % Create color table to represent phase
    NumHues = 256;
    KeyHues = [0,  1, 1;...  % H,S,V
        1,  1, 1];...
        HueMap = zeros(NumHues,1);
    tmpx=[0:size(HueMap,1)-1]/(size(HueMap,1)-1)*(size(KeyHues,1)-1)+1;
    HueMap(:,1) = interp1([1:size(KeyHues,1)],KeyHues(:,1),tmpx).';
    HueMap(:,2) = interp1(KeyHues(:,2),tmpx).';
    HueMap(:,3) = 1;                 % Brightness
    clear tmpx KeyHues NumHues;
    cmap = hsv2rgb(HueMap);
    clear HueMap;
end

% Reshape matrix to [no. elements, 1] to simplify handling of arbitrary no. dimensions
OrigSize = size(data);					% Remember the size
data = reshape(data, [prod(OrigSize),1]);

% Calculate phase and magnitude
phase=(angle(data)+pi)/(2*pi);        	  % in [0,1]
phase=round(phase*(size(cmap,1)-1)+1);  % round to nearest entry in hue map
data = abs(data);
tmpval = max(data);
if tmpval~=0,
    data = data./tmpval;
    if ShowLogMagnitude,
        tmpThreshold = 0.0001;
        data = log(data + tmpThreshold);
        TmpHiVal = log(1+tmpThreshold);
        TmpLoVal = min(data);
        if TmpLoVal==TmpHiVal,
            TmpHiVal=TmpHiVal+0.5;
            TmpLoVal=TmpLoVal-0.5;
        end
        data = (data-TmpLoVal)./(TmpHiVal-TmpLoVal);
        clear TmpLoVal TmpHiVal tmpThreshold;
    end
end
clear tmpval;

% Convert first to HSV colors, then to RGB colors
RgbImg = zeros([prod(OrigSize),3],'uint8');   	  % Matrix size = [no. elements, 3]
RgbImg(:,1) = uint8(255*cmap(phase,1).*abs(data));	% Hue
RgbImg(:,2) = uint8(255*cmap(phase,2).*abs(data));	% Saturation
RgbImg(:,3) = uint8(255*cmap(phase,3).*abs(data)); 					% Brightness to range from 0 to 1
clear data phase;
%RgbImg = hsv2rgb(RgbImg);   % Convert to RGB colors
%RgbImg = uint8(RgbImg*255); % Convert to 8-bit to save memory
%clear HueMap;
clear cmap;

% Restore size
RgbImg = reshape(RgbImg, [OrigSize, 3]);
clear OrigSize ShowLogMagnitude;
end
