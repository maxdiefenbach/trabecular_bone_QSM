function psiHat = bestInitMap(psiHatSet, mask)

[ht, wd, N] = size(psiHatSet);

scores = zeros(N, 1);

for ctr = 1:N,
    
    score = 0;
    currPsiHat = psiHatSet(:,:,ctr);
    currPsiHat = filter2(ones(3)/9, currPsiHat, 'same');
    
    for k=1:ht,
        for l=1:wd,
            if mask(k, l),
                % check the first-order 4 neighbors
                dum = 0;
                if mask(max(1, k-1), l)
                    dum = dum + abs(currPsiHat(k, l)-...
                        currPsiHat(max(1, k-1), l));
                end
                if mask(k, max(1, l-1))
                    dum = dum + abs(currPsiHat(k, l)-...
                        currPsiHat(k, max(1, l-1)));
                end
                if mask(min(ht, k+1), l)
                    dum = dum + abs(currPsiHat(k, l)-...
                        currPsiHat(min(ht, k+1), l));
                end
                if mask(k, min(wd, l+1))
                    dum = dum + abs(currPsiHat(k, l)-...
                        currPsiHat(k, min(wd, l+1)));
                end
                score = score + dum;
            end
        end
    end
    
    scores(ctr) = score;
end

[dum, mInd] = min(scores);
psiHat = psiHatSet(:, :, mInd);

return;