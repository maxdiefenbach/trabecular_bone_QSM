classdef TestCSS < matlab.unittest.TestCase
    
    properties
        TE_s
        centerFreq_Hz
        paramsMatrix
        signal
    end
    
    methods (TestClassSetup)
        function addPath(testCase)
            addpath(genpath(pwd))
            addpath(genpath('~/Programs/dev/BMRR/'))
        end
        function make_mex(testCase)
            codegen CSS
        end
        function set_properties(testCase)
            TE_s = (1:3) * 1e-3;
            centerFreq_Hz = 128e6;
            ModelParams.Species(1).name = 'water';
            ModelParams.Species(1).magnitude = 80;
            ModelParams.Species(1).phase = pi/4;
            ModelParams.Species(1).freq_Hz = 10;
            ModelParams.Species(1).R2s_Hz = 0;
            ModelParams.Species(2).name = 'fat';
            ModelParams.Species(2).magnitude = 20;
            ModelParams.Species(2).phase = pi/4;
            ModelParams.Species(2).freq_Hz = -3.5e-6 * centerFreq_Hz;
            ModelParams.Species(2).R2s_Hz = 0;
            paramsMatrix = ModelParams2paramsMatrix(ModelParams);
            % paramsMatrix = get_paramsMatrix(ModelParams);
            signal = build_signal(TE_s, paramsMatrix);
            
            testCase.TE_s = TE_s;
            testCase.centerFreq_Hz = centerFreq_Hz;
            testCase.paramsMatrix = paramsMatrix;
            testCase.signal = signal;
        end
    end

    methods (Test)
        
        
        function test_build_signal(testCase)
            TE_s = testCase.TE_s;
            paramsMatrix = testCase.paramsMatrix;
            signal = build_signal(TE_s, paramsMatrix);
        end


        function test_get_CSS_Amatrix(testCase)
            TE_s = testCase.TE_s;
            paramsMatrix = testCase.paramsMatrix;
            A = get_CSS_Amatrix(TE_s, paramsMatrix);
            signal = build_signal(TE_s, paramsMatrix);
            
            rhoComplex = A \ signal
            
            testCase.verifyEqual(abs(rhoComplex), paramsMatrix(1:2, 1), 'AbsTol', 1e-6);
            testCase.verifyEqual(angle(rhoComplex), paramsMatrix(1:2, 2), 'AbsTol', 1e-6);
        end
        
        
        function test_compute_multiSpeciesJacobian(testCase)
            TE_s = testCase.TE_s;
            centerFreq_Hz = testCase.centerFreq_Hz;
            paramsMatrix = testCase.paramsMatrix;
            
            C_rho = eye(2);
            C_phi = eye(2);
            C_omega = [1, 0; 1, 0]; %eye(2);
            C_r = eye(2);
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r)
            
            J = compute_multiSpeciesJacobian(TE_s, paramsMatrix, constrainMatrices)
            
            nTE = length(TE_s);
            nSpecies = size(paramsMatrix, 1);
            testCase.verifySize(J, [nTE, 4 * nSpecies]);
        end

        
        function test_CSSvoxel(testCase)
            signal = testCase.signal;
            TE_s = testCase.TE_s;
            centerFreq_Hz = testCase.centerFreq_Hz;
            paramsMatrix = testCase.paramsMatrix;
            
            C_rho = eye(2);
            C_phi = eye(2);
            C_omega = [1, 0; 1, 0]; %eye(2);
            C_r = [1, 0; 1, 0];
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);
            
            tol = 1e-4;
            iterMax = 100;
            [paramsMatrix, residual, iter, precision] = ...
                CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax)
            
            AbsTol = 1e-3;
            testCase.verifyEqual(paramsMatrix, testCase.paramsMatrix, 'AbsTol', AbsTol);
        end

        
        function test_CSS(testCase)
            signalArray = transpose(repmat(testCase.signal, [1, 10]));
            TE_s = testCase.TE_s;
            centerFreq_Hz = testCase.centerFreq_Hz;
            paramsMatrixArray = shiftdim(repmat(testCase.paramsMatrix, [1, 1, 10]), 2);
            
            expected = paramsMatrixArray;

            C_rho = eye(2);
            C_phi = eye(2);
            C_omega = [1, 0; -1, 0]; %eye(2);
            C_r = [1, 0; 1, 0]; %eye(2);
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);

            tol = 1e-3;
            iterMax = 100;
            [paramsMatrixArray, residualArray, iterArray, precisionArray] = CSS(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax);
            
            AbsTol = 1e-3;
            testCase.verifyEqual(paramsMatrixArray, expected, 'AbsTol', AbsTol);
        end

        
        function test_CSS_mex(testCase)
            nVoxels = 1000;
            signalArray = transpose(repmat(testCase.signal, [1, nVoxels]));
            TE_s = testCase.TE_s;
            centerFreq_Hz = testCase.centerFreq_Hz;
            paramsMatrixArray = shiftdim(repmat(testCase.paramsMatrix, [1, 1, nVoxels]), 2);
            
            expected = paramsMatrixArray;

            C_rho = eye(2);
            C_phi = eye(2);
            C_omega = [1, 0; -1, 0]; %eye(2);
            C_r = [1, 0; 1, 0]; %eye(2);
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);
            
            tol = 1e-3;
            iterMax = 100;
            skip_threshold_percent = 5;
            [paramsMatrixArray, residualArray, iterArray, precisionArray, nSkippedVoxels] = CSS_mex(signalArray, TE_s, paramsMatrixArray, constrainMatrices, tol, iterMax, skip_threshold_percent);

            AbsTol = 1e-3;
            testCase.verifyEqual(paramsMatrixArray, expected, 'AbsTol', AbsTol);
        end
        
        
        function test_fatSpectrum(testCase)
            relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016]; % 1;
            freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83]; % -3.4; 
            
            freqs_Hz = freqs_ppm * 1e-6 * testCase.centerFreq_Hz;
            nSpecies = 1 + length(relAmps);

            paramsMatrix = zeros([nSpecies, 4]);
            paramsMatrix(1, 1) = 50;
            paramsMatrix(1, 2) = pi/4;
            paramsMatrix(1, 3) = 10;
            paramsMatrix(1, 4) = 2;
            paramsMatrix(2:end, 1) = (100-paramsMatrix(1, 1)) .* relAmps';
            paramsMatrix(2:end, 2) = pi/4;
            paramsMatrix(2:end, 3) = paramsMatrix(1, 3) + freqs_ppm';
            paramsMatrix(2:end, 4) = 2;
            
            nSpecies = size(paramsMatrix, 1);
            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = relAmps;
            C_phi = double(logical(C_rho));
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = -1;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);

            TE_s = testCase.TE_s;
            centerFreq_Hz = testCase.centerFreq_Hz;
            signal = build_signal(TE_s, paramsMatrix);
            
            tol = 1e-3;
            iterMax = 100;
            
            tic
            [paramsMatrixR, residual, iter, precision] = CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax);
            toc

            paramsMatrixR
            residual
            iter
            precision
            
            AbsTol = 1e-3;
            testCase.verifyEqual(paramsMatrixR, paramsMatrix, 'AbsTol', AbsTol);
        end

        
        function test_WFI(testCase)
            ModelParams.FatModel.relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016]; % 1;
            ModelParams.FatModel.freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83]; % -3.4; 
            
            nSpecies = 8;
            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = ModelParams.FatModel.relAmps;
            ModelParams.C_rho = C_rho;
            C_phi = double(logical(C_rho));
            % C_phi = zeros([nSpecies, nSpecies]);
            % C_phi(:, 1) = 1;
            ModelParams.C_phi = C_phi;
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = -1;
            ModelParams.C_omega = C_omega;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            ModelParams.C_r = C_r;
            
            fieldmap_Hz = 10;
            matrixSize = [10, 10, 10];
            nVoxels = prod(matrixSize);
            ModelParams.fieldmap_Hz = fieldmap_Hz * ones(matrixSize);
            R2s_Hz = 2;
            ModelParams.R2s_Hz = R2s_Hz * ones(matrixSize);
            
            centerFreq_Hz = 42.58e6 * 3;

            paramsMatrix = zeros([nSpecies, 4]);
            paramsMatrix(1, 1) = 50;
            paramsMatrix(1, 2) = 0;
            paramsMatrix(1, 3) = fieldmap_Hz;
            paramsMatrix(1, 4) = R2s_Hz;
            paramsMatrix(2:end, 1) = (100-paramsMatrix(1, 1)) .* ModelParams.FatModel.relAmps';
            paramsMatrix(2:end, 2) = 0;
            paramsMatrix(2:end, 3) = paramsMatrix(1, 3) + ModelParams.FatModel.freqs_ppm' * 1e-6 * centerFreq_Hz;
            paramsMatrix(2:end, 4) = 2;
            
            paramsMatrix

            ImDataParams.TE_s = testCase.TE_s;
            ImDataParams.centerFreq_Hz = centerFreq_Hz;
            signalVoxel = build_signal(ImDataParams.TE_s, paramsMatrix);
            signal = squeeze(permute(repmat(signalVoxel, [1, 1, nVoxels]), [3, 1, 2]));
            size(signal)
            nTE = length(ImDataParams.TE_s);
            ImDataParams.signal = reshape(signal, [matrixSize, nTE]);
            
            WFIparams = separate_waterFat(ImDataParams, ModelParams)
            
            abs(WFIparams.water(1, 1))
            angle(WFIparams.fat(1, 1))
        end

        
        function test_invivo(testCase)
            warning('off', 'all')
            
            load 20160216_133028_0202_ImDataParams.mat % ./20160216_134233_0402_ImDataParams.mat;
            ImDataParams
            
            ImDataParams.signal = ImDataParams.signal(:, :, 5:6, :);
            
            imagine(ImDataParams.signal)

            ModelParams.FatModel.relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016]; % 1;
            ModelParams.FatModel.freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83]; % -3.4; 
            
            nSpecies = 8;
            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = ModelParams.FatModel.relAmps;
            ModelParams.C_rho = C_rho;
            C_phi = double(logical(C_rho));
            % C_phi = zeros([nSpecies, nSpecies]);
            % C_phi(:, 1) = 1;
            ModelParams.C_phi = C_phi;
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = 1;
            ModelParams.C_omega = C_omega;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            ModelParams.C_r = C_r;
            ModelParams
            
            warning('off', 'all');
            WFIparams = separate_waterFat(ImDataParams, ModelParams);
            warning('on', 'all');
            
            imagine(WFIparams.water, 'Window', [0, 1000], ...
                    WFIparams.fat, 'Window', [0, 1000], ...
                    WFIparams.fieldmap_Hz, 'Window', [-300, 300], ...
                    WFIparams.R2s_Hz, 'Window', [0, 300], ...
                    WFIparams.residual, ...
                    WFIparams.iterations);
            
            AlgoParams.FatModel = ModelParams.FatModel
            AlgoParams.precision = 1e-4;
            AlgoParams.iterMax = 100;
            AlgoParams.fieldmap_Hz = 0;
            AlgoParams.R2s_Hz = 0;
            ImDataParams.signal = double(ImDataParams.signal);
            WFIparams = IDEALmex(ImDataParams, AlgoParams);
            
            imagine(WFIparams.water, 'Window', [0, 1000], ...
                    WFIparams.fat, 'Window', [0, 1000], ...
                    WFIparams.fieldmap_Hz, 'Window', [-300, 300], ...
                    WFIparams.R2s_Hz, 'Window', [0, 300], ...
                    WFIparams.residual, ...
                    WFIparams.iterations);
            
            warning('off', 'all')
        end

        
        function compare_Amatrix(testCase)
            
            TE_s = (1:3) * 1e-3;
            centerFreq_Hz = 12.8e6;
            relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016];
            freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83];
            freqs_Hz = freqs_ppm * 1e-6 * centerFreq_Hz;
            
            nSpecies = length(freqs_Hz) + 1;
            
            paramsMatrix = zeros(nSpecies, 4);
            paramsMatrix(2:end, 3) = freqs_Hz;

            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = relAmps;
            unknowRhoMask = logical(diag(C_rho));


            A = get_CSS_Amatrix(TE_s, paramsMatrix);
            P = diag(squeeze(exp(1j .* paramsMatrix(:, 2))));
            Ar = (A * P) * C_rho(:, unknowRhoMask)
            
            ImDataParams = add_vars2struct(struct(), TE_s, centerFreq_Hz);
            FatModel = add_vars2struct(struct(), relAmps, freqs_ppm);
            IDEAL_A = get_Amatrix(ImDataParams, FatModel)
            
            testCase.verifyEqual(Ar, IDEAL_A, 'AbsTol', 1e-10);
            
        end

        
        function compare_Bmatrix(testCase)
            TE_s = (1:3) * 1e-3;
            centerFreq_Hz = 12.8e6;
            relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016];
            freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83];
            freqs_Hz = freqs_ppm * 1e-6 * centerFreq_Hz;
            
            nSpecies = length(freqs_Hz) + 1;
            
            paramsMatrix = zeros(nSpecies, 4);
            paramsMatrix(1, 1) = 60;
            paramsMatrix(1, 3) = 110;
            paramsMatrix(1, 4) = 10;
            paramsMatrix(2:end, 1) = 40 * relAmps;
            paramsMatrix(2:end, 3) = 110 + freqs_Hz;
            paramsMatrix(2:end, 4) = 10;
            paramsMatrix

            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = relAmps;
            unknowRhoMask = logical(diag(C_rho));
            
            C_phi = double(logical(C_rho))
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = 1;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);

            signal = build_signal(TE_s, paramsMatrix);
            
            paramsMatrix(1, 3) = 0;
            paramsMatrix(2:end, 3) = 0 + freqs_Hz;
            paramsMatrix(:, 4) = 0;
            A = get_CSS_Amatrix(TE_s, paramsMatrix);
            P = diag(squeeze(exp(1j .* paramsMatrix(:, 2))));
            Ar = (A * P) * C_rho(:, unknowRhoMask);
            
            rho = Ar \ signal
            
            paramsMatrix(1, 1) = abs(rho(1));
            paramsMatrix(1, 2) = angle(rho(1));
            paramsMatrix(2:end, 1) = relAmps .* abs(rho(2));
            paramsMatrix(2:end, 2) = angle(rho(2));

            J = compute_multiSpeciesJacobian(TE_s, paramsMatrix, constrainMatrices);
            unknownParamsMask = logical(cat(1, zeros(nSpecies, 1), diag(C_phi), diag(C_omega), diag(C_r)));
            unknownParamsMask = logical(cat(1, zeros(2*nSpecies, 1), diag(C_omega), diag(C_r)));
            J = J(:, unknownParamsMask)
            
            s = build_signal(TE_s, paramsMatrix);
            ds = signal - s;
            updates  = split_ReIm(J) \ split_ReIm(ds)
            
        %% ideal
            ImDataParams = add_vars2struct(struct(), TE_s, centerFreq_Hz);
            FatModel = add_vars2struct(struct(), relAmps, freqs_ppm);
            A = get_Amatrix(ImDataParams, FatModel);
            ModelParams.fieldmap_Hz = 0;
            ModelParams.R2s_Hz = 0;
            M = update_modelMatrix(A, TE_s, ModelParams);
            rho = M \ signal
            B = get_Bmatrix(TE_s, M, rho, ModelParams)
            
            signalModel = M * rho;
            dS = split_ReIm(signal - signalModel);
            B = split_ReIm(get_Bmatrix(TE_s, M, rho, ModelParams));
            updatesArray = B(:, 3:end) \ dS
            
            split_ReIm(J)
            split_ReIm(B)
        end


        function compare_updates(testCase)
            TE_s = (1:3) * 1e-3;
            centerFreq_Hz = 12.8e6;
            relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016];
            freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83];
            freqs_Hz = freqs_ppm * 1e-6 * centerFreq_Hz;
            
            nSpecies = length(freqs_Hz) + 1;
            
            paramsMatrix = zeros(nSpecies, 4);
            paramsMatrix(1, 1) = 60;
            paramsMatrix(1, 3) = 110;
            paramsMatrix(1, 4) = 10;
            paramsMatrix(2:end, 1) = 40 * relAmps;
            paramsMatrix(2:end, 3) = 110 + freqs_Hz;
            paramsMatrix(2:end, 4) = 10;
            paramsMatrix

            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = relAmps;
            unknowRhoMask = logical(diag(C_rho));
            
            C_phi = double(logical(C_rho));
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = -1;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);

            signal = build_signal(TE_s, paramsMatrix);
            
            paramsMatrix(1, 3) = 0;
            paramsMatrix(2:end, 3) = 0 + freqs_Hz;
            paramsMatrix(:, 4) = 0;
            A = get_CSS_Amatrix(TE_s, paramsMatrix);
            P = diag(squeeze(exp(1j .* paramsMatrix(:, 2))));
            Ar = A * C_rho(:, unknowRhoMask);
            
            rho = Ar \ signal
            
            paramsMatrix
            paramsMatrix(1, 1) = abs(rho(1));
            paramsMatrix(1, 2) = angle(rho(1));
            paramsMatrix(2:end, 1) = relAmps .* abs(rho(2));
            paramsMatrix(2:end, 2) = angle(rho(2));
            paramsMatrix

            s = build_signal(TE_s, paramsMatrix);
            ds = signal - s;
            nRho = sum(logical(diag(C_rho)));
            nPhi = sum(logical(diag(C_phi)));

            J = compute_multiSpeciesJacobian(TE_s, paramsMatrix, constrainMatrices);
            unknownParamsMask = logical(cat(1, zeros(2*nSpecies, 1), diag(C_omega), diag(C_r)));
            unknownParamsMask(unknownParamsMask~=0)
            J = J(:, unknownParamsMask(unknownParamsMask~=0))
            updates = J \ ds
            abs(updates)
            % updates(updates~=0)
            
        %% ideal
            
            ImDataParams = add_vars2struct(struct(), TE_s, centerFreq_Hz);
            FatModel = add_vars2struct(struct(), relAmps, freqs_ppm);
            A = get_Amatrix(ImDataParams, FatModel);
            ModelParams.fieldmap_Hz = 0;
            ModelParams.R2s_Hz = 0;
            M = update_modelMatrix(A, TE_s, ModelParams);
            rho = M \ signal
            signalModel = M * rho;
            dS = split_ReIm(signal - signalModel);
            B = split_ReIm(get_Bmatrix(TE_s, M, rho, ModelParams));
            updatesArray = B \ dS
            
        end


        function compare2IDEAL(testCase)
            warning('off', 'all')
            

            load 20160216_133028_0202_ImDataParams.mat
            matrixSize = get_matrixSize(ImDataParams);
            ix = 46; % 60;
            iy = 56; % 60;
            iz = 5;
            TE_s = ImDataParams.TE_s;
            signal = squeeze(ImDataParams.signal(ix, iy, iz, :));
            centerFreq_Hz = ImDataParams.centerFreq_Hz;

            relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016];
            freqs_ppm = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83];
            freqs_Hz = freqs_ppm * 1e-6 * centerFreq_Hz;
            
            nSpecies = length(freqs_Hz) + 1;
            
            paramsMatrix = zeros(nSpecies, 4);
            paramsMatrix(2:end, 3) = freqs_Hz;

            C_rho = zeros([nSpecies, nSpecies]);
            C_rho(1, 1) = 1;
            C_rho(2:end, 2) = relAmps;
            C_phi = double(logical(C_rho));
            C_omega = zeros([nSpecies, nSpecies]);
            C_omega(1, 1) = 1;
            C_omega(2:end, 1) = 1;
            C_r = zeros([nSpecies, nSpecies]);
            C_r(1, 1) = 1;
            C_r(2:end, 1) = 1;
            constrainMatrices = cat(3, C_rho, C_phi, C_omega, C_r);

            tol = 1e-4;
            iterMax = 100;
            
            magn = norm(signal)

            A = get_CSS_Amatrix(TE_s, paramsMatrix);
            unknowRhoMask = logical(diag(C_rho));
            Ar = A * C_rho(:, unknowRhoMask)
            [paramsMatrix, residual, iter, precision] = ...
                CSSvoxel(signal, TE_s, paramsMatrix, constrainMatrices, tol, iterMax);
            paramsMatrix
            (paramsMatrix(:, 1)' * logical(C_rho(:, 1:2)))'
            iter
            residual
            
            ModelParams.fieldmap_Hz = 0;
            ModelParams.R2s_Hz = 0;
            AlgoParamsVoxel.precision = tol;
            AlgoParamsVoxel.iterMax = iterMax;
            A = get_Amatrix(ImDataParams, add_vars2struct(struct(), relAmps, freqs_ppm))
            
            OutParamsVoxel = IDEALvoxel(signal, TE_s, A, ModelParams, AlgoParamsVoxel);
            
            wm = abs(OutParamsVoxel.water)
            wp = angle(OutParamsVoxel.water)
            fm = abs(OutParamsVoxel.fat)
            fp = angle(OutParamsVoxel.fat)
            fieldmap_Hz = OutParamsVoxel.fieldmap_Hz
            R2s_Hz = OutParamsVoxel.R2s_Hz
            iterations = OutParamsVoxel.iterations
            delta = OutParamsVoxel.delta
            residual = OutParamsVoxel.residual
            
            warning('on', 'all')
        end
        
        
        function test_apply_linearUpdates(testCase)
            rhoComplex = [1 + 2j; 3 + 4j]
            
            abs(rhoComplex)
            angle(rhoComplex)

            C_rho = [1, 0, 0; ...
                     0, 0.5, 0; ...
                     0, 0.5, 0;];
            C_phi = [1, 0, 0; ...
                     0, 1, 0; ...
                     0, 1, 0;];
            
            paramsMatrix = zeros(3, 4)
            
            paramsMatrix = apply_linearUpdates(paramsMatrix, rhoComplex, C_rho, C_phi)
        end

        
        function test_apply_nonlinearUpdates(testCase)
            updates = [-1; -2; -10; -20; -3; 0.1];
            
            C_rho = [1, 0, 0; ...
                     0, 0.5, 0; ...
                     0, 0.5, 0;];
            C_phi = [1, 0, 0; ...
                     0, 1, 0; ...
                     0, 1, 0;];
            C_omega = [1, 0, 0; ...
                       1, 0, 0; ...
                       1, 0, 0];
            C_r = [1, 0, 0; ...
                   1, 0, 0; ...
                   1, 0, 0];
            
            paramsMatrix = ones(3, 4);
            
            paramsMatrix = apply_nonlinearUpdates(paramsMatrix, updates, C_rho, C_phi, C_omega, C_r)
        end

    end

end