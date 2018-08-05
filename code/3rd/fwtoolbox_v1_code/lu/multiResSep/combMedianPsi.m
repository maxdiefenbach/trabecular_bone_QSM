function psi = combMedianPsi(psiHat, mask)

psi = psiHat(:,:,2);

[ht, wd] = size(psi);

% sliding window for local field map fitting
swSiz = 3; swStep = 3;
qmap = zeros((ht-swSiz)/swStep+1, (wd-swSiz)/swStep+1);
for k=1:(ht-swSiz)/swStep+1,
    for l=1:(wd-swSiz)/swStep+1,
        
        swVrange = (k-1)*swStep+1:(k-1)*swStep+swSiz;
        swHrange = (l-1)*swStep+1:(l-1)*swStep+swSiz;
        
        localPsi = psi(swVrange, swHrange);
        qmap(k, l) = std(localPsi(:));
    end
end
