function WFIparams = separate_waterFat(ImDataParams, ModelParams, AlgoParams)
% WFIparams = separate_waterFat(ImDataParams, ModelParams, AlgoParams)
% 
    % ImDataParams
    signalSize = size(ImDataParams.signal);
    matrixSize = signalSize(1:3);
    nTE = signalSize(4);
    if isfield(AlgoParams, 'iSlice')
        matrixSize(3) = length(AlgoParams.iSlice);
    end
    TE_s = ImDataParams.TE_s(:)';
    if isfield(AlgoParams, 'iTE')
        TE_s = TE_s(AlgoParams.iTE);
        nTE = length(AlgoParams.iTE);
    end
    centerFreq_Hz = ImDataParams.centerFreq_Hz;
    
    % ModelParams
    fieldmap_Hz = set_option(ModelParams, 'fieldmap_Hz', zeros(matrixSize));
    R2s_Hz = set_option(ModelParams, 'R2s_Hz', zeros(matrixSize));
    waterR2s_Hz = set_option(ModelParams, 'waterR2s_Hz', R2s_Hz);
    fatR2s_Hz = set_option(ModelParams, 'fatR2s_Hz', R2s_Hz);
    phase_rad = set_option(ModelParams, 'phase_rad', zeros(matrixSize));
    ModelParams.FatModel = set_option(ModelParams, 'FatModel', struct());
    % default FatModel like on the scanner
    relAmps = set_option(ModelParams.FatModel, 'relAmps', [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016]);
    freqs_ppm = set_option(ModelParams.FatModel, 'freqs_ppm', -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83]);
    freqs_Hz = freqs_ppm .* 1e-6 .* centerFreq_Hz;
    freqs_Hz = [0, freqs_Hz];
    nSpecies = length(freqs_Hz);
    % constrainMatrices (default model: 7 fat peaks, single T2*)
    C_rho = zeros(nSpecies, nSpecies);
    C_rho(1, 1) = 1;
    C_rho(2:end, 2) = relAmps;
    C_phi = set_option(ModelParams, 'C_phi', double(logical(C_rho)));
    C_omega = zeros(nSpecies, nSpecies);
    C_omega(:, 1) = 1;
    C_omega = set_option(ModelParams, 'C_omega', C_omega);
    C_r = zeros(nSpecies, nSpecies);
    C_r(:, 1) = 1;
    C_r = set_option(ModelParams, 'C_r', C_r);
    constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);
    
    % AlgoParams
    if nargin < 3
        AlgoParams = struct();
    end
    tol = set_option(AlgoParams, 'tol', 1e-3);
    iterMax = set_option(AlgoParams, 'iterMax', 100);
    skipVoxelThreshold_percent = set_option(AlgoParams, 'skipVoxelThreshold_percent', 5);
    iSlice = set_option(AlgoParams, 'iSlice', 1:matrixSize(3));
    iTE = set_option(AlgoParams, 'iTE', 1:nTE);
    signalArray = double(flatten3d(ImDataParams.signal(:, :, iSlice, iTE)));
    nVoxels = size(signalArray, 1);
    
    % paramsMatrixArray
    paramsMatrixArray = zeros(nVoxels, nSpecies, 4);
    % water
    paramsMatrixArray(:, 1, 3) = fieldmap_Hz(:);
    paramsMatrixArray(:, 1, 4) = waterR2s_Hz(:);
    paramsMatrixArray(:, 1, 2) = phase_rad(:);
    % fat
    for iSpecies = 2:nSpecies
        paramsMatrixArray(:, iSpecies, 3) = fieldmap_Hz(:) + freqs_Hz(iSpecies);
        paramsMatrixArray(:, iSpecies, 4) = fatR2s_Hz(:);
        paramsMatrixArray(:, 1, 2) = phase_rad(:);
    end
    
    % run CSS
    tic
    [paramsMatrixArray, residualArray, iterArray, precisionArray, nSkippedVoxels] = ...
        CSS_mex(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax, skipVoxelThreshold_percent);
    elapsedTime_s = toc;
    
    % reshape output
    water = reshape(paramsMatrixArray(:, 1, 1) .* exp(1j * paramsMatrixArray(:, 1, 2)), matrixSize);
    fat = reshape(sum(paramsMatrixArray(:, 2:end, 1), 2) .* exp(1j * paramsMatrixArray(:, 2, 2)), matrixSize);
    fieldmap_Hz = reshape(paramsMatrixArray(:, 1, 3), matrixSize);
    R2s_Hz = reshape(paramsMatrixArray(:, 1, 4), matrixSize);
    residual = reshape(residualArray, matrixSize);
    iterations = reshape(iterArray, matrixSize);
    
    % pack output struct
    WFIparams = struct();
    WFIparams = add_vars2struct(WFIparams, water, fat, fieldmap_Hz);
    if isfield(ModelParams, 'waterR2s_Hz')
        waterR2s_Hz = R2s_Hz;
        WFIparams = add_vars2struct(WFIparams, waterR2s_Hz);
    else
        WFIparams = add_vars2struct(WFIparams, R2s_Hz);
    end
    if isfield(ModelParams, 'fatR2s_Hz')
        fatR2s_Hz = reshape(paramsMatrixArray(:, 2, 4), matrixSize);
        WFIparams = add_vars2struct(WFIparams, fatR2s_Hz);
    end
    WFIparams = add_vars2struct(WFIparams, paramsMatrixArray, residual, iterations);
end