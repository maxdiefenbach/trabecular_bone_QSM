function Extrema = find_Extrema(array, iDim, delta)
    
    forwardDiff = finite_difference(array, iDim, 'forward', 'Dirichlet');
    backwardDiff = finite_difference(array, iDim, 'backward', 'Dirichlet');
    
    maxMask = forwardDiff < delta & backwardDiff < delta;
    minMask = forwardDiff > delta & backwardDiff > delta;
    
    indMax = find(maxMask);
    indMin = find(minMask);
    
    subMax = [];
    for iInd = 1:length(indMax)
        sub = myInd2sub(size(array), indMax(iInd));
        subMax = [subMax; sub];
    end
    subMin = [];
    for iInd = 1:length(indMin)
        sub = myInd2sub(size(array), indMin(iInd));
        subMin = [subMin; sub];
    end

    Extrema.Maxima.ind = indMax;
    Extrema.Maxima.sub = subMax;
    Extrema.Maxima.val = array(indMax);
    Extrema.Minima.ind = indMin;
    Extrema.Minima.sub = subMin;
    Extrema.Minima.val = array(indMin);

end