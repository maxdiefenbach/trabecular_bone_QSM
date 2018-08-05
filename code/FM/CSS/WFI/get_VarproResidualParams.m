function ResidualParams = get_VarproResidualParams(signalTrue, ImDataParams, FatModel, Options)
    
    if nargin < 4
        Options = struct();
    end

    TE_s = ImDataParams.TE_s;
    dTE_s = diff(TE_s(1:2));
    centerFreq_Hz = ImDataParams.centerFreq_Hz;
    
    fatFractionRange_percent = set_option(Options, 'fatFractionRange_percent', [0, 1, 2, 10:10:90, 98, 99, 100]);
    freqWrapStep_Hz = 1/dTE_s;
    fieldmapRange_Hz = set_option(Options, 'fieldmapRange_Hz', (-2:0.01:2) * freqWrapStep_Hz);
    % compute in range [-2pi, 2pi] (phase_rad pi <=> fieldmapRange_Hz 1/2/dTE_s)
    waterR2sRange_Hz = set_option(Options, 'waterR2sRange_Hz', 0:10:300);
    waterR2sRange_Hz = set_option(Options, 'R2s_Hz', waterR2sRange_Hz);
    fatR2sRange_Hz = set_option(Options, 'fatR2sRange_Hz', 0:10:300);
    fatR2sRange_Hz = set_option(Options, 'R2s_Hz', fatR2sRange_Hz);
    ModelParamsVary = build_ModelParams(centerFreq_Hz, FatModel, fatFractionRange_percent, fieldmapRange_Hz, waterR2sRange_Hz, fatR2sRange_Hz);
    nPeaks = length(FatModel.freqs_ppm);
    skipStepMat = ones(nPeaks+1, 4);
    skipStepMat(1, :) = 0;
    skipStepMat = set_option(Options, 'skipStepMat', skipStepMat);
    additionalPhase_rad = set_option(Options, 'additionalPhase_rad', zeros(size(TE_s(:)')));
    
    [residualArray, paramsMatrixArray] = compute_varpro_residual_mex(TE_s(:)', ...
                                                      signalTrue, ...
                                                      ModelParamsVary, ...
                                                      skipStepMat, ...
                                                      additionalPhase_rad);
    C_rho = zeros(nPeaks+1, nPeaks+1);  % fat + water
    C_rho(1, 1) = 1;
    C_rho(2, 2:end) = FatModel.relAmps;
    warning('off', 'Coder:MATLAB:rankDeficientMatrix')
    warning('off', 'Coder:MATLAB:singularMatrix')
    paramsMatrixArray = compute_linParamsLandscape_mex(TE_s(:).', signalTrue, paramsMatrixArray, C_rho);
    warning('on', 'Coder:MATLAB:rankDeficientMatrix')
    warning('on', 'Coder:MATLAB:singularMatrix')

    [countMatrix, nParamCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
    residualSize = [length(fatFractionRange_percent), length(fieldmapRange_Hz), length(waterR2sRange_Hz)];
    residualLandscape = reshape(residualArray, residualSize);
    paramsMatrixArray = reshape(paramsMatrixArray, [residualSize, nPeaks + 1, 4]);

    subMin = get_subMinResidual(residualLandscape);
    FF = fatFractionRange_percent(subMin(1));
    fB = fieldmapRange_Hz(subMin(2));
    R2s = waterR2sRange_Hz(subMin(3));

    ResidualParams.TE_s = TE_s;
    ResidualParams.signal = signalTrue;
    ResidualParams.FatModel = FatModel;
    ResidualParams.fatFractionRange_percent = fatFractionRange_percent;
    ResidualParams.fieldmapRange_Hz = fieldmapRange_Hz;
    ResidualParams.waterR2sRange_Hz = waterR2sRange_Hz;
    ResidualParams.fatR2sRange_Hz = fatR2sRange_Hz;
    ResidualParams.residualLandscape = residualLandscape;
    ResidualParams.paramsMatrixArray = paramsMatrixArray;
    ResidualParams.subMin = subMin;
    ResidualParams.fatFraction_percent = FF;
    ResidualParams.fieldmap_Hz = fB;
    ResidualParams.R2s_Hz = R2s;
end