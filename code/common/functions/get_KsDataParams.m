function KsDataParams =  get_KsDataParams(ImDataParams)
    
    signalSize = size(ImDataParams.signal);
    matrixSize = signalSize(1:3);
    nTE = signalSize(4);
    
    Signal = complex(zeros(signalSize));
    for iTE = 1:nTE
        s = squeeze(ImDataParams.signal(:, :, :, iTE));
        Signal(:, :, :, iTE) = ifftshift(fftn(s));
    end

    KsDataParams = rmfield(ImDataParams, 'signal');
    KsDataParams.Signal = Signal;
end