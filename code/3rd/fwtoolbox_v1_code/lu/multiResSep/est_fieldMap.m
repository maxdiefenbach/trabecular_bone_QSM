function psiHat = est_fieldMap(TE, S, df, dfAmp, vlevel)
% Estimate the field map at the coarsest level using multiple golden
% section search
% Inputs:
%        TE (seconds)
%         S (image data)
% Output:
%     psiHat (estimated field map)
%
% Two constraints:
% 1) Minimize the cost function
% 2) Ensure the overall smoothness of the field map
%


%% initialization
[ht, wd, necho] = size(S);
mask = imMagMask(S, .2);

pInitMap = zeros(ht, wd);
pp = zeros(ht, wd, 3);

A = calc_mixMat(TE, df, dfAmp );
deltaTE = abs(TE(2)-TE(1));
psiRange = 1/deltaTE;

%% 1) Minimize the cost function
for k=1:ht,
    for l=1:wd,
        if vlevel>1, disp([num2str(k), ', ', num2str(l)]); end
        
        if mask(k, l),
        s = S(k, l, :); s = s(:); % take one pixel's image data
        tol = 1e-5; % threshold for golden section search
        
        % multiple golden-section search
        [pHat, rMin, tpsi] = multiGolden(psiRange, s, A, TE, tol, vlevel);
        
        % initial field map contains the value minimizing the cost function
        pInitMap(k, l) = pHat; 
        % pp contains three possible field map values
        pp(k, l, :) = tpsi;     
        end

    end
end

%% 2) Ensure the overall smoothness of the field map
psiHat = smooth_fmap(pInitMap, pp, psiRange, mask, vlevel);

return;