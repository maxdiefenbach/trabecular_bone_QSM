function ImDataParamsPlusNoise = addNoise(ImDataParams, SNR)
    
    ImDataParamsPlusNoise = ImDataParams
    matrixSize = get_matrixSize(ImDataParams);
    nDynamics = size(ImDataParams.signal, 4);
    signalSize = [matrixSize, nDynamics];
    noise = 1/SNR .* (randn(signalSize) + 1j * randn(signalSize));
    ImDataParamsPlusNoise.signal = ImDataParamsPlusNoise.signal + noise;

end