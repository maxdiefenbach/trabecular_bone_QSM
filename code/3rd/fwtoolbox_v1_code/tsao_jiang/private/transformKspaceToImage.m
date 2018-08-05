%
% [newdata] = transformKspaceToImage(data,dim,DCasFirstElement,NormalizeAsInV)
%
% - To transform k-space data to an image by Inverse Fourier Transform.
% - Input Parameters
%     data: k-space data
%     dim: dimension(s) to transform, if not given, transform all dimensions
%     DCasFirstElement: by default, the DC term of image and k-space are
%                       considered to be in the center (bitshift(datasize,-1)+1).
%                       If DCasFirstElement==1, DC is considered to be at the first
%                       element.
%     NormalizeAsInV: 1 to normalize Fourier transform according to the
%                     convention of V (default). 0 to use native convention, which
%                     may be variable from computer to computer.
%

% VERSION:
%   1.0 - March 27, 2000 Jeffrey Tsao
%   Sept 27, 2002 - Check input parameters for empty (Jeffrey Tsao jtsao2@hotmail.com)
function [newdata] = transformKspaceToImage(data,dim,DCasFirstElement,NormalizeAsInV)
if nargin<1, help mfilename; end
if nargin<2, dim=[]; end
if nargin<3, DCasFirstElement=[]; end
if nargin<4, NormalizeAsInV=[]; end
if isempty(DCasFirstElement), DCasFirstElement = 0; end
if isempty(NormalizeAsInV),   NormalizeAsInV   = 1; end

numDims = ndims(data);   % Number of dimensions

%----------------------------------------------------------------------
% Move DC from center to 1st element, if needed.
%   (fft routine assumes DC to be 1st element)
if DCasFirstElement,
    newdata = data;
else                          % Move DC to 1st element for fft
    idx = cell(1, numDims);
    for k = 1:numDims
        m = size(data, k);
        if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
            p = bitshift(m,-1)+1;   % central pixel (i.e. position of DC)
            idx{k} = [p:m 1:p-1];
            clear p
        else
            idx{k} = [1:m];
        end
    end
    clear k m
    newdata = data(idx{:});      % Perform fft-shift
end
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Perform fft
fftnumelements = 1;
for k = 1:numDims
    m = size(newdata, k);
    if m>1 & (isempty(dim) | ~isempty(find(k==dim))), % either transform all dimensions or just transform this dimension
        newdata = ifft(newdata,[],k);
        fftnumelements = fftnumelements*m;  % Count total number of elements in fft
    end
end
clear k m numDims
%----------------------------------------------------------------------

%----------------------------------------------------------------------
% Move DC from 1st element back to center, if necessary
if DCasFirstElement==0,
    newdata(idx{:}) = newdata;  % Perform fft-shift
    clear idx
end
%----------------------------------------------------------------------
clear numDims

%----------------------------------------------------------------------
% Normalize according to V convention, if necessary
if NormalizeAsInV,
    newdata = newdata*fftnumelements;
end
%----------------------------------------------------------------------

end
