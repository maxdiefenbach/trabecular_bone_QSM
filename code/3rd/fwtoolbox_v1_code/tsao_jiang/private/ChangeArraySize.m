%
% [ResizedData] = ChangeArraySize(Data,NewSize)
%
% - Change the size of Data, by cropping or zero-padding if needed.
% - Center of data is assumed at floor(N/2) where N is the array size.

% April 24, 2003 - Created by Jeffrey Tsao
function [ResizedData] = ChangeArraySize(Data,NewSize)
if nargin<1, help(mfilename); end
if nargin<2, NewSize=[]; end
if isempty(NewSize), ResizedData = Data; clear Data; return; end

NewSize = NewSize(:);
NumDims = ndims(Data);
NumDims = max(NumDims, length(NewSize));

% Set up arrays for size and indices
if length(NewSize)==1, NewSize = ones(NumDims,1)*NewSize; end
DesArraySize = zeros(1,NumDims);
SrcIdx = cell(1, NumDims);
DesIdx = cell(1, NumDims);
for DimNum = [1:NumDims],
  SrcArraySize = size(Data,DimNum);
  if DimNum<=length(NewSize),
    DesArraySize(DimNum) = NewSize(DimNum);
  else
    DesArraySize(DimNum) = SrcArraySize;
  end
  
  MinSize = min(SrcArraySize, DesArraySize(DimNum));
  tmpidx = [1:MinSize]-bitshift(MinSize,-1); clear MinSize;
  SrcIdx{DimNum} = tmpidx + bitshift(SrcArraySize        ,-1);
  DesIdx{DimNum} = tmpidx + bitshift(DesArraySize(DimNum),-1);
  clear tmpidx SrcArraySize;
end; clear DimNum NumDims;

% Resize
ResizedData = zeros(DesArraySize); clear DesArraySize;
ResizedData(DesIdx{:}) = Data(SrcIdx{:}); clear SrcIdx DesIdx;
clear Data NewSize;
