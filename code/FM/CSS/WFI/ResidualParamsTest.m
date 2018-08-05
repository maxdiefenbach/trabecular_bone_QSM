classdef ResidualParamsTest < matlab.unittest.TestCase
    
    methods (Test)
        
        function test_get_ResidualParams(testCase)
            
            signalTrue = 1.0e+08 * ...
                [1.3269 - 1.0434i  -0.5881 - 1.0531i  -0.6935 + 0.7569i  -0.6252 + 0.0697i ...
                 -0.0013 + 0.6820i  -0.2109 + 0.4374i   0.3545 + 0.3302i   0.2259 + 0.3788i ...
                 0.3527 - 0.1086i   0.3636 + 0.0380i  -0.2244 - 0.1827i  -0.2268 + 0.1676i ...
                 -0.2256 + 0.0406i  -0.0206 + 0.2090i  -0.0812 + 0.1580i   0.1027 + 0.1158i ...
                 0.0444 + 0.1336i   0.1133 + 0.0096i   0.0905 + 0.0480i   0.0109 + 0.0829i];
            ImDataParams.TE_s = [0.0008 0.0019 0.0029 0.0040 0.0050 0.0061 0.0072 0.0082 0.0093 0.0104 0.0114 0.0125 0.0136 0.0146 0.0157 0.0168 0.0178 0.0189 0.0199 0.0210];
            ImDataParams.centerFreq_Hz = 127754321;
            FatModel.freqs_ppm = [-3.8000   -3.4000   -3.1000   -2.6800   -2.4600   -1.9500   -0.5000 0.4900    0.5900]
            FatModel.relAmps = [0.0880    0.6420    0.0580    0.0620    0.0580    0.0060    0.0390 0.0100    0.0370];
            FatModel.freqs_ppm = -3.4
            FatModel.relAmps = 1

            ResidualParams = get_ResidualParams(signalTrue, ImDataParams, FatModel)
            mysurf(squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :)), ...
                   ResidualParams.fieldmapRange_Hz, ...
                   ResidualParams.waterR2sRange_Hz, ...
                   200, 25)
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))
        end
        
        
        function test_simSignal(testCase)
            F = FatModel;
            fatFraction_percent = 30;
            F.add_waterPeak(fatFraction_percent)
            F.build_signal
            F.plot_signal
            
            ResidualParams = get_ResidualParams(signalTrue, ImDataParams, FatModel)
            mysurf(squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :)), ...
                   ResidualParams.fieldmapRange_Hz, ...
                   ResidualParams.waterR2sRange_Hz, ...
                   200, 25)
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))
        end

        
        function test_residual(testCase)
            GYRO = 42.58e6;                % [Hz/T]
            nTE = 3;
            dTE_s = 1.0e-3;
            TEmin_s = 1.0e-3;
            TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            centerFreq_Hz = GYRO * fieldStrength_T;
            
            ModelParamsTrue.Species(1).magnitude = 30;
            ModelParamsTrue.Species(1).phase = pi/4;
            ModelParamsTrue.Species(1).freq_Hz = 0;
            ModelParamsTrue.Species(1).R2s_Hz = 50;
            ModelParamsTrue.Species(2).magnitude = 70;
            ModelParamsTrue.Species(2).phase = pi/4;
            ModelParamsTrue.Species(2).freq_Hz = 0 - 3.4e-6 * centerFreq_Hz;
            ModelParamsTrue.Species(2).R2s_Hz = 50;

            ModelParamsVary.Species(1).magnitude = 30;
            ModelParamsVary.Species(1).phase = pi/4;
            ModelParamsVary.Species(1).freq_Hz = (0:0.01:1)/(dTE_s);
            ModelParamsVary.Species(1).R2s_Hz = 0:10:300;
            ModelParamsVary.Species(2).magnitude = 70;
            ModelParamsVary.Species(2).phase = pi/4;
            ModelParamsVary.Species(2).freq_Hz = ModelParamsVary.Species(1).freq_Hz - 3.4e-6 * centerFreq_Hz;
            ModelParamsVary.Species(2).R2s_Hz = 0:10:300;
            
            skipStepMat = [0, 0, 0, 0; ...
                           1, 1, 1, 1];
            
            additionalPhase_rad = zeros(size(TE_s));
            residualArray = compute_residual(TE_s, ModelParamsTrue, ModelParamsVary, skipStepMat, additionalPhase_rad);
            [countMatrix, nCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            shape = countMatrix(countMatrix~=1)';
            residualArray = reshape(residualArray, shape);
            
            close all;
            x = ModelParamsVary.Species(1).freq_Hz;
            y = ModelParamsVary.Species(2).R2s_Hz;
            mysurf(residualArray, x, y, 200, 25)
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            colormap viridis

        end
        
        
        function test_singlePeak(testCase)
            F = FatModel('single peak');
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1.0e-3;
            F.TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            GYRO = 42.58e6;
            F.centerFreq_Hz = GYRO * fieldStrength_T;
            F.fatFraction_percent = 70;
            F.R2s_Hz = 50;
            F.add_waterPeak;
            F.freqs_ppm
            F.build_signal
            F.plot_signal
            FM = F.get_FatModel;
            
            ImDataParams.TE_s = F.TE_s;
            ImDataParams.centerFreq_Hz = F.centerFreq_Hz;
            
            ResidualParams = get_ResidualParams(F.signal, ImDataParams, FM)
            mysurf(squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :)), ...
                   ResidualParams.fieldmapRange_Hz, ...
                   ResidualParams.waterR2sRange_Hz, ...
                   200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))

            figure;
            xlabel('fieldmap [Hz]')
            ylabel('residual [a.u.]')
            plot(ResidualParams.fieldmapRange_Hz(ResidualParams.subMin(2)), squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, ResidualParams.subMin(3))));
            
            figure;
            xlabel('R2s [1/s]');
            ylabel('residual [a.u.]')
            plot(ResidualParams.waterR2sRange_Hz, squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), ResidualParams.subMin(2), :)));

            align_allFigures
        end

        
        function test_compute_residual(testCase)
            GYRO = 42.58e6;                % [Hz/T]
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1e-3;
            TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            centerFreq_Hz = GYRO * fieldStrength_T;
            
            ModelParamsTrue.Species(1).magnitude = 30;
            ModelParamsTrue.Species(1).phase = 0;
            ModelParamsTrue.Species(1).freq_Hz = 0;
            ModelParamsTrue.Species(1).R2s_Hz = 50;
            ModelParamsTrue.Species(2).magnitude = 70;
            ModelParamsTrue.Species(2).phase = 0;
            ModelParamsTrue.Species(2).freq_Hz = 0 - 3.4e-6 * centerFreq_Hz;
            ModelParamsTrue.Species(2).R2s_Hz = 50;

            ModelParamsVary.Species(1).magnitude = 30;
            ModelParamsVary.Species(1).phase = 0;
            ModelParamsVary.Species(1).freq_Hz = -500:1:500;
            ModelParamsVary.Species(1).R2s_Hz = 0:100;
            ModelParamsVary.Species(2).magnitude = 70;
            ModelParamsVary.Species(2).phase = 0;
            ModelParamsVary.Species(2).freq_Hz = (-500:1:500) - 3.4e-6 * centerFreq_Hz;
            ModelParamsVary.Species(2).R2s_Hz = 0:100;
            
            skipStepMat = [0, 0, 0, 0; ...
                           1, 1, 1, 1];
            
            additionalPhase_rad = zeros(size(TE_s));
            residualArray = compute_residual(TE_s, ModelParamsTrue, ModelParamsVary, skipStepMat, additionalPhase_rad);
            [countMatrix, nCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            shape = countMatrix(countMatrix~=1)';
            residualArray = reshape(residualArray, shape);
            
            x = ModelParamsVary.Species(1).freq_Hz;
            y = ModelParamsVary.Species(2).R2s_Hz;
            mysurf(residualArray, x, y, 200, 25)
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            colormap viridis
            
            figure;
            xlabel('fieldmap [Hz]')
            ylabel('residual [a.u.]')
            plot(ModelParamsVary.Species(1).freq_Hz, residualArray(:, 50));
            
            figure;
            xlabel('R2s [1/s]');
            ylabel('residual [a.u.]')
            plot(ModelParamsVary.Species(1).R2s_Hz, residualArray(500, :));

            align_allFigures
        end

        
        function compare_residual12(testCase)
            GYRO = 42.58e6;                % [Hz/T]
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1e-3;
            TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            centerFreq_Hz = GYRO * fieldStrength_T;
            
            %%
            F = FatModel('single peak');
            F.TE_s = TE_s;
            fieldStrength_T = 3;
            F.centerFreq_Hz = GYRO * fieldStrength_T;
            F.fatFraction_percent = 70;
            F.fieldmap_Hz = 200;
            F.R2s_Hz = 150;
            F.add_waterPeak;
            F.build_signal
            FM = F.get_FatModel;
            
            ImDataParams.TE_s = F.TE_s;
            ImDataParams.centerFreq_Hz = F.centerFreq_Hz;
            
            ResidualParams = get_ResidualParams(F.signal, ImDataParams, FM)
            mysurf(squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :)), ...
                   ResidualParams.fieldmapRange_Hz, ...
                   ResidualParams.waterR2sRange_Hz, ...
                   200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))

            
            %%            
            ModelParamsTrue.Species(1).magnitude = 100 - F.fatFraction_percent;
            ModelParamsTrue.Species(1).phase = 0;
            ModelParamsTrue.Species(1).freq_Hz = F.fieldmap_Hz;
            ModelParamsTrue.Species(1).R2s_Hz = F.R2s_Hz;
            ModelParamsTrue.Species(2).magnitude = F.fatFraction_percent;
            ModelParamsTrue.Species(2).phase = 0;
            ModelParamsTrue.Species(2).freq_Hz = ModelParamsTrue.Species(1).freq_Hz - 3.4e-6 * centerFreq_Hz;
            ModelParamsTrue.Species(2).R2s_Hz = F.R2s_Hz;

            ModelParamsVary.Species(1).magnitude = 100 - F.fatFraction_percent;
            ModelParamsVary.Species(1).phase = 0;
            ModelParamsVary.Species(1).freq_Hz = (0:0.01:2) / dTE_s;
            ModelParamsVary.Species(1).R2s_Hz = 0:10:300;
            ModelParamsVary.Species(2).magnitude = F.fatFraction_percent;
            ModelParamsVary.Species(2).phase = 0;
            ModelParamsVary.Species(2).freq_Hz = ModelParamsVary.Species(1).freq_Hz -3.4e-6 * centerFreq_Hz;
            ModelParamsVary.Species(2).R2s_Hz = 0:10:300;
            
            signalTrue = build_signal(TE_s, ModelParams2paramsMatrix(ModelParamsTrue));
            
            testCase.verifyEqual(signalTrue(:).', F.signal);

            skipStepMat = [0, 0, 0, 0; ...
                           1, 1, 1, 1];
            
            additionalPhase_rad = zeros(size(TE_s));
            residualArray = compute_residual(TE_s, ModelParamsTrue, ModelParamsVary, skipStepMat, additionalPhase_rad);
            [countMatrix, nCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            shape = [1, countMatrix(countMatrix~=1)'];
            residualArray = reshape(residualArray, shape);
            
            x = ModelParamsVary.Species(1).freq_Hz;
            y = ModelParamsVary.Species(2).R2s_Hz;
            mysurf(squeeze(residualArray(1, :, :)), x, y, 200, 25)
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            colormap viridis
            
            testCase.verifyEqual(residualArray, ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :), 'AbsTol', 1e-9)

            align_allFigures
        end

        
        function test_build_ModelParams(testCase)
            GYRO = 42.58e6;                % [Hz/T]
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1e-3;
            TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            centerFreq_Hz = GYRO * fieldStrength_T;
            
            %%
            F = FatModel('single peak');
            F.TE_s = TE_s;
            fieldStrength_T = 3;
            GYRO = 42.58e6;
            F.centerFreq_Hz = GYRO * fieldStrength_T;
            F.fatFraction_percent = 70;
            F.fieldmap_Hz = 200;
            F.R2s_Hz = 150;
            F.add_waterPeak;
            F.build_signal
            FM = F.get_FatModel

            fatFraction_percent = 70;
            freqWrapStep_Hz = 1/dTE_s;
            fieldmap_Hz = (0:0.01:2) * freqWrapStep_Hz; % compute in range [0, 2pi] (phase_rad pi <=> fieldmap_Hz 1/2/dTE_s)
            waterR2s_Hz = 0:10:300;
            fatR2s_Hz = waterR2s_Hz;
            ModelParamsVary1 = build_ModelParams(centerFreq_Hz, FM, fatFraction_percent, fieldmap_Hz, waterR2s_Hz, fatR2s_Hz);
            
            
            ModelParamsTrue.Species(1).magnitude = 100 - F.fatFraction_percent;
            ModelParamsTrue.Species(1).phase = 0;
            ModelParamsTrue.Species(1).freq_Hz = 0;
            ModelParamsTrue.Species(1).R2s_Hz = F.R2s_Hz;
            ModelParamsTrue.Species(2).magnitude = F.fatFraction_percent;
            ModelParamsTrue.Species(2).phase = 0;
            ModelParamsTrue.Species(2).freq_Hz = ModelParamsTrue.Species(1).freq_Hz - 3.4e-6 * centerFreq_Hz;
            ModelParamsTrue.Species(2).R2s_Hz = F.R2s_Hz;

            ModelParamsVary2.Species(1).magnitude = 100 - F.fatFraction_percent;
            ModelParamsVary2.Species(1).phase = 0;
            ModelParamsVary2.Species(1).freq_Hz = (0:0.01:4) / dTE_s;
            ModelParamsVary2.Species(1).R2s_Hz = 0:10:300;
            ModelParamsVary2.Species(2).magnitude = F.fatFraction_percent;
            ModelParamsVary2.Species(2).phase = 0;
            ModelParamsVary2.Species(2).freq_Hz = ModelParamsVary2.Species(1).freq_Hz- 3.4e-6 * centerFreq_Hz;
            ModelParamsVary2.Species(2).R2s_Hz = 0:10:300;

            ModelParamsVary1
            ModelParamsVary1.Species(1)
            ModelParamsVary1.Species(2)
            ModelParamsVary2
            ModelParamsVary2.Species(1)
            ModelParamsVary2.Species(2)
            
            testCase.verifyEqual(ModelParamsVary1, ModelParamsVary2, 'AbsTol', 1e-10)
        end
        
        
        function test_compute_residual2(testCase)
            GYRO = 42.58e6;                % [Hz/T]
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1e-3;
            TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            fieldStrength_T = 3;
            centerFreq_Hz = GYRO * fieldStrength_T;
            
            F = FatModel('single peak');
            F.TE_s = TE_s;
            fieldStrength_T = 3;
            F.centerFreq_Hz = GYRO * fieldStrength_T;
            F.fatFraction_percent = 70;
            F.fieldmap_Hz = 200;
            F.R2s_Hz = 150;
            F.add_waterPeak;
            F.build_signal
            
            ImDataParams.TE_s = F.TE_s;
            ImDataParams.centerFreq_Hz = F.centerFreq_Hz;
            
            ResidualParams = get_ResidualParams(F.signal, ImDataParams, F.get_FatModel)
            h = mysurf(squeeze(ResidualParams.residualLandscape(ResidualParams.subMin(1), :, :)), ...
                       ResidualParams.fieldmapRange_Hz, ...
                       ResidualParams.waterR2sRange_Hz, ...
                       200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('residual [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))

            
            h = mysurf(squeeze(ResidualParams.paramsMatrixArray(ResidualParams.subMin(1), :, :, 1, 1)), ...
                       ResidualParams.fieldmapRange_Hz, ...
                       ResidualParams.waterR2sRange_Hz, ...
                       200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('water [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))

            h = mysurf(squeeze(ResidualParams.paramsMatrixArray(ResidualParams.subMin(1), :, :, 2, 1)), ...
                       ResidualParams.fieldmapRange_Hz, ...
                       ResidualParams.waterR2sRange_Hz, ...
                       200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('fat [a.u.]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))


            ffLandscape = 100 * ...
                squeeze(ResidualParams.paramsMatrixArray(ResidualParams.subMin(1), :, :, 2, 1) ./ ...
                        (ResidualParams.paramsMatrixArray(ResidualParams.subMin(1), :, :, 1, 1) + ResidualParams.paramsMatrixArray(ResidualParams.subMin(1), :, :, 2, 1)));
            h = mysurf(ffLandscape, ...
                       ResidualParams.fieldmapRange_Hz, ...
                       ResidualParams.waterR2sRange_Hz, ...
                       200, 25);
            xlabel('fieldmap [Hz]')
            ylabel('R2s [1/s]')
            zlabel('fat fraction [%]')
            title(sprintf('fat fraction %d', ResidualParams.fatFraction_percent))

            align_allFigures
        end
        
        
        function test_varpro_residual(testCase)
            
            F = FatModel('single peak');
            F.TE_s = (1:3) * 1e-3;
            F.fieldmap_Hz = 100;
            F.R2s_Hz = 100;
            F.add_waterPeak(50);
            F.build_signal;
            F.plot_signal;
            
            TE_s = F.TE_s;
            signal = F.signal;
            ModelParamsVary = build_ModelParams(42.58e6 * 3, ...
                                                F.get_FatModel, ...
                                                F.fatFraction_percent, ...
                                                -1000:10:1000, ...
                                                0:25:250, ...
                                                0:25:250, ...
                                                0, ...
                                                0);
            nSpecies = numel(ModelParamsVary.Species)
            skipStepMat = ones(nSpecies, 4);
            skipStepMat(1, :) = 0;
            additionalPhase_rad = zeros(size(signal));
            
            tic
            [residualArray, paramsMatrixArraySubMin, subMin, residualMin] = compute_varpro_residual_mex(TE_s, double(signal(:)).', ModelParamsVary, ...
                                                              skipStepMat, additionalPhase_rad(:)');
            toc
            
            [countMatrix, nParamCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            residualSize = countMatrix(countMatrix~=0 & countMatrix~=1)'
            residualLandscape = reshape(residualArray, residualSize);
            
            figure
            plot(1:size(residualLandscape, 1), residualLandscape(:, 10))
            yyaxis right
            grid on
            plot(1:size(residualLandscape, 1)-1, diff(residualLandscape(:, 10)))
            hold on
            plot(1:size(residualLandscape, 1)-2, diff(residualLandscape(:, 10), 2))
            legend('varpro residual', '1st derivative slope', '2nd derivative curvature')

        end

        
        function test_find_locMinima(testCase)
            F = FatModel('single peak');
            F.TE_s = (1:3) * 1e-3;
            F.fieldmap_Hz = 100;
            F.R2s_Hz = 100;
            F.add_waterPeak(50);
            F.build_signal;
            
            TE_s = F.TE_s;
            signal = F.signal;
            ModelParamsVary = build_ModelParams(42.58e6 * 3, ...
                                                F.get_FatModel, ...
                                                F.fatFraction_percent, ...
                                                -1000:10:1000, ...
                                                0:25:250, ...
                                                0:25:250, ...
                                                0, ...
                                                0);
            nSpecies = numel(ModelParamsVary.Species)
            skipStepMat = ones(nSpecies, 4);
            skipStepMat(1, :) = 0;
            additionalPhase_rad = zeros(size(signal));
            
            tic
            [residualArray, paramsMatrixArraySubMin, subMin, residualMin] = compute_varpro_residual_mex(TE_s, double(signal(:)).', ModelParamsVary, ...
                                                              skipStepMat, additionalPhase_rad(:)');
            toc
            
            [countMatrix, nParamCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            residualSize = countMatrix(countMatrix~=0 & countMatrix~=1)'
            residualLandscape = reshape(residualArray, residualSize);
            
            [maxtab, mintab] = peakdet(residualLandscape(:, 10), 1, -1000:10:1000)
            [maxtab, mintab] = peakdet(residualLandscape(:, 10), 1)
            
            R = residualLandscape(:, 10).'
            Extrema = find_Extrema(R, 1)
            Extrema.Minima
            Extrema.Maxima
                        
            figure
            plot(R)
            hold on
            plot(Extrema.Maxima.ind, Extrema.Maxima.val, 'o')
            plot(Extrema.Minima.ind, Extrema.Minima.val, 'o')
        end
        
        
        function test_find_Extrema(testCase)
            F = FatModel('single peak');
            F.TE_s = (1:10) * 1e-3;
            F.fieldmap_Hz = 100;
            F.R2s_Hz = 100;
            F.add_waterPeak(50);
            F.build_signal;
            
            TE_s = F.TE_s;
            signal = F.signal;
            ModelParamsVary = build_ModelParams(42.58e6 * 3, ...
                                                F.get_FatModel, ...
                                                F.fatFraction_percent, ...
                                                -1000:10:1000, ...
                                                0:25:250, ...
                                                0:25:250, ...
                                                0, ...
                                                0);
            nSpecies = numel(ModelParamsVary.Species)
            skipStepMat = ones(nSpecies, 4);
            skipStepMat(1, :) = 0;
            additionalPhase_rad = zeros(size(signal));
            
            tic
            [residualArray, paramsMatrixArraySubMin, subMin, residualMin] = compute_varpro_residual_mex(TE_s, double(signal(:)).', ModelParamsVary, ...
                                                              skipStepMat, additionalPhase_rad(:)');
            toc
            
            [countMatrix, nParamCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            residualSize = countMatrix(countMatrix~=0 & countMatrix~=1)';
            residualLandscape = reshape(residualArray, residualSize);
            
            Extrema = find_Extrema(residualLandscape, 1, 1);
            
            close all;
            h = mysurf(residualLandscape, -1000:10:1000, 0:25:250, 500, 50);
            hold on
            scatter3(Extrema.Minima.sub(:, 1), Extrema.Minima.sub(:, 2), Extrema.Minima.val, 100, 'red', 'filled')
            scatter3(Extrema.Maxima.sub(:, 1), Extrema.Maxima.sub(:, 2), Extrema.Maxima.val, 100, 'green', 'filled')

        end

        
        function test_shift_array_1D(testCase)
            array = 1:3;
            testCase.verifyEqual(shift_array(array, 2, 'forward', 'periodic'), [2, 3, 1])
            testCase.verifyEqual(shift_array(array, 2, 'backward', 'periodic'), [3, 1, 2])
            testCase.verifyEqual(shift_array(array, 2, 'forward', 'Dirichlet'), [2, 3, 0])
            testCase.verifyEqual(shift_array(array, 2, 'backward', 'Dirichlet'), [0, 1, 2])
            testCase.verifyEqual(shift_array(array, 2, 'forward', 'replicate'), [2, 3, 3])
            testCase.verifyEqual(shift_array(array, 2, 'backward', 'replicate'), [1, 1, 2])
        end

        
        function test_shift_array_2D(testCase)
            array = repmat((1:3).', [1, 3]);
            testCase.verifyEqual(shift_array(array, 1, 'forward', 'periodic'), repmat([2, 3, 1].', [1, 3]))
            testCase.verifyEqual(shift_array(array, 1, 'backward', 'periodic'), repmat([3, 1, 2].', [1, 3]))
            testCase.verifyEqual(shift_array(array, 1, 'forward', 'Dirichlet'), repmat([2, 3, 0].', [1, 3]))
            testCase.verifyEqual(shift_array(array, 1, 'backward', 'Dirichlet'), repmat([0, 1, 2].', [1, 3]))
            testCase.verifyEqual(shift_array(array, 1, 'forward', 'replicate'), repmat([2, 3, 3].', [1, 3]))
            testCase.verifyEqual(shift_array(array, 1, 'backward', 'replicate'), repmat([1, 1, 2].', [1, 3]))
        end

        
        function test_finite_differences_1D(testCase)
            array = 1:3;
            testCase.verifyEqual(finite_difference(array, 2, 'forward', 'Dirichlet'), [1, 1, -3])
            testCase.verifyEqual(finite_difference(array, 2, 'backward', 'Dirichlet'), [-1, -1, -1])
            testCase.verifyEqual(finite_difference(array, 2, 'forward', 'periodic'), [1, 1, -2])
            testCase.verifyEqual(finite_difference(array, 2, 'backward', 'periodic'), [2, -1, -1])
            testCase.verifyEqual(finite_difference(array, 2, 'forward', 'replicate'), [1, 1, 0])
            testCase.verifyEqual(finite_difference(array, 2, 'backward', 'replicate'), [0, -1, -1])
        end

        
        function test_finite_differences_2D(testCase)
            array = repmat((1:3).', [1, 3]);
            testCase.verifyEqual(finite_difference(array, 1, 'forward', 'Dirichlet'), repmat([1, 1, -3].', [1, 3]))
            testCase.verifyEqual(finite_difference(array, 1, 'backward', 'Dirichlet'), repmat([-1, -1, -1].', [1, 3]))
            testCase.verifyEqual(finite_difference(array, 1, 'forward', 'periodic'), repmat([1, 1, -2].', [1, 3]))
            testCase.verifyEqual(finite_difference(array, 1, 'backward', 'periodic'), repmat([2, -1, -1].', [1, 3]))
            testCase.verifyEqual(finite_difference(array, 1, 'forward', 'replicate'), repmat([1, 1, 0].', [1, 3]))
            testCase.verifyEqual(finite_difference(array, 1, 'backward', 'replicate'), repmat([0, -1, -1].', [1, 3]))
        end

        
        function test_HernandoResidual(testCase)
        % Function: computeResidual
        %
        % Description: compute fit residual for water/fat imaging
        % 
        % Parameters:
        %% Input: structures imDataParams and algoParams
        %%   - imDataParams.images: acquired images, array of size[nx,ny,1,ncoils,nTE]
        %%   - imDataParams.TEs: echo times (in seconds)
        %%   - imDataParams.fieldStrength: (in Tesla)
        %%
        %%   - algoParams.species(ii).frequency = frequency shift in ppm of each peak within species ii
        %%   - algoParams.species(ii).relAmps = relative amplitude (sum normalized to 1) of each peak within species ii
        %%   Example
        %%      - algoParams.species(1).name = 'water' % Water
        %%      - algoParams.species(1).frequency = [0] 
        %%      - algoParams.species(1).relAmps = [1]   
        %%      - algoParams.species(2).name = 'fat' % Fat
        %%      - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
        %%      - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
        %% 
        %%   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
        %%   - algoParams.range_r2star = [0 0]; % Range of R2* values
        %%   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
        %%   - algoParams.range_fm = [-400 400]; % Range of field map values
        %%   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
        %%   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
        %%   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
        %%   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
        %%   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
        %%   - algoParams.lambda = 0.05; % Regularization parameter
        %%   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
        %%   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
        %%   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 
        %
        % Returns: 
        %  - residual: the residual, of size NUM_FMS X sx X sy
            
            I = ImDataParamsBMRR('/Users/maxdiefenbach/programs/BMRR/tests/testdata/20160216WaterBottlePhantom/2016_02_16/20160216_133028_0202_ImDataParams.mat')
            imDataParams = I.get_imDataParams4toolbox;
            imDataParams.images = imDataParams.images(:, :, ceil((size(imDataParams.images, 3)+1)/2), :, :)
            
            I.set_FatModel('single peak')
            algoParams = I.get_algoParams4toolbox;
            algoParams.range_fm = [-600, 600];
            algoParams.NUM_FMS = 1200;
            algoParams.range_r2star = [0, 250];
            algoParams.NUM_R2STARS = 250;
            
            imDataParams
            algoParams
            
            tic
            residual = computeResidual(imDataParams, algoParams);
            toc

            figure; 
            plot(1:1200, squeeze(residual(:, 60, 45)));
            
            TE_s = imDataParams.TE;
            signal = squeeze(imDataParams.images(60, 45, 1, 1, :));
            ModelParamsVary = build_ModelParams(I.ImDataParams.centerFreq_Hz, ...
                                                I.AlgoParams.FatModel, ...
                                                50, ...
                                                linspace(algoParams.range_fm(1), algoParams.range_fm(2), algoParams.NUM_FMS), ...
                                                linspace(algoParams.range_r2star(1), algoParams.range_r2star(2), algoParams.NUM_R2STARS), ...
                                                linspace(algoParams.range_r2star(1), algoParams.range_r2star(2), algoParams.NUM_R2STARS), ...
                                                0, 0);
            nSpecies = numel(ModelParamsVary.Species)
            skipStepMat = ones(nSpecies, 4);
            skipStepMat(1, :) = 0;
            additionalPhase_rad = zeros(size(signal));
            
            tic
            [residualArray] = compute_varpro_residual(TE_s, double(signal(:)).', ModelParamsVary, ...
                                                      skipStepMat, additionalPhase_rad(:)');
            toc
            tic
            [residualArray, paramsMatrixArraySubMin, subMin, residualMin] = compute_varpro_residual_mex(TE_s, double(signal(:)).', ModelParamsVary, ...
                                                              skipStepMat, additionalPhase_rad(:)');
            toc
            
            paramsMatrixArraySubMin

            [countMatrix, nParamCounts] = get_ModelParamsCount(ModelParamsVary, skipStepMat);
            residualSize = countMatrix(countMatrix~=0 & countMatrix~=1)'
            residualLandscape = reshape(residualArray, residualSize);
            size(residualLandscape)
            
            figure
            plot(1:size(residualLandscape, 1), residualLandscape(:, 10))
            
            align_allFigures
        end

    end % methods (Test)
    
end % classdef