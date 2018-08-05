function matrixSize = get_matrixSize(ImDataParams)
    
    signal = ImDataParams.signal;
    signalDims = ndims(signal);
    signalSize = size(signal);

    if signalDims == 2
        matrixSize = [signalSize(1), 1, 1];
       
    elseif signalDims == 3
        matrixSize = [signalSize(1), signalSize(2), 1];

    else
        matrixSize = signalSize(1:(signalDims-1));
        
    end

end % EOF