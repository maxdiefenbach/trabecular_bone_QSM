function A = calc_mixMat(TE, df, dfAmp)
%
% TE: a vector containing multiple echo times (for example: 3x1)
%
% Example:
% TE =
% 
%     0.0016    0.0032    0.0048
% 
% >> A = calc_mixMat(TE, -210, 1)
% 
% A =
% 
%    1.0000            -0.5979 - 0.8016i
%    1.0000            -0.2850 + 0.9585i
%    1.0000             0.9387 - 0.3446i

TE = TE(:);
N = size(TE,1); % # of echoes

if nargin<2,
    df = -210;
end
if nargin<3,
    dfAmp = 1;
end

c = zeros(N, 1);
for k=1:N,
  c(k) = sum(exp(i*2*pi*df.*TE(k)).*dfAmp);
end

A = cat(2, ones(N, 1), c(:));
