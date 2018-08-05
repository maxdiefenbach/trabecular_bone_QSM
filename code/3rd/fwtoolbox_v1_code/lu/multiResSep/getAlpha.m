function alpha = getAlpha(dat)
%
% Obtain the factor of linear phase correction for bi-directional
% multi-echo sequences -- alpha
% 
% Output: alpha (linear phase factor)
% Input: dat (ht, wd, nechoes) image data
%
% See: get_linPhaseMat
% 
% Wenmiao Lu Nov. 18, 2006
%

%% get the k-space data
K = cat(3, ft(dat(:,:,1)), ft(dat(:,:,2)), ft(dat(:,:,3)));

%% use golden section search for alpha
alphaRange = -.095; % be careful about how to get k-space data (FT
                   % or IFT?
tol = 1e-5;
if alphaRange>0,
  
alpha = goldenAlpha('alphaFunc', 0, alphaRange/2, alphaRange, tol, ...
                    K, 1);
else
    
alpha = goldenAlpha('alphaFunc', alphaRange, alphaRange/2, 0, tol, ...
                    K, 1);
end