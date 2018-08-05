function [w, f, resd] = lsSep_wf(TE, s, psi, df, dfAmp )
% least-square separation given TE and field map

[ht, wd, necho] = size(s);

% demodulate field map
sc = [];
for k=1:necho,
    tmp = s(:,:,k).*exp(-i*2*pi*psi.*TE(k));
    sc = cat(2, sc, tmp(:));
end
% MATLAB transpose is conjugate transpose
sc = sc.'; 

A = calc_mixMat(TE, df, dfAmp );
wf = A\sc;
resd = sc - A*wf;

w = wf(1, :); f = wf(2, :);
w = reshape(w, [ht wd]);
f = reshape(f, [ht wd]);

resd = abs(resd);      
resd = permute(resd, [2 3 1]);
resd = reshape(resd, [ht wd size(resd, 3)]);
resd = sum(resd, 3);
          
return;
