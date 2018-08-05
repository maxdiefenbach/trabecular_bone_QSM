classdef FatModelTest < matlab.unittest.TestCase
    
    methods(Test)

        function test_FatModel(testCase)
            close all;
            F = FatModel('Hamilton liver')
            F.get_paramsMatrix;
            F.build_signal;
            F.sort_peaks('freqs_ppm')
            F.plot_spectrum;
            h = F.plot_signal;
            ylim([0, 110])
            hold all;
            
            F.set_FatModel('single peak');
            F.build_signal;
            plot(F.TE_s * 1e3, abs(F.signal), 'LineWidth', 2)

            F.add_waterPeak(70)
            F.build_signal;
            plot(F.TE_s * 1e3, abs(F.signal), 'LineWidth', 2)

            F.plot_speciesVectors
            align_allFigures
        end
        
        
        function test_buildSignal(testCase)
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

            signalTrue = build_signal(TE_s, ModelParams2paramsMatrix(ModelParamsTrue));
            signalTrue = signalTrue(:).'
            
            F = FatModel('single peak');
            F.fieldmap_Hz = 0;
            F.R2s_Hz = 50;
            F.TE_s = TE_s;
            F.centerFreq_Hz = centerFreq_Hz;
            F.add_waterPeak(70);
            F.get_paramsMatrix
            F.build_signal;
            F.signal
            
        % flipped
        % testCase.verifyEqual(F.get_paramsMatrix, ModelParams2paramsMatrix(ModelParamsTrue));
            testCase.verifyEqual(F.signal, signalTrue);
        end

        
        function test_paramsMatrix(testCase)
            nTE = 3;
            dTE_s = 1.6e-3;
            TEmin_s = 1e-3;
            
            %%
            F = FatModel('single peak');
            F.TE_s = TEmin_s + (0:(nTE-1)) * dTE_s;
            F.fieldStrength_T = 3;
            GYRO = 42.58e6;
            F.centerFreq_Hz = GYRO * F.fieldStrength_T;
            F.fatFraction_percent = 70;
            F.R2s_Hz = 50;
            F.add_waterPeak;
            F.build_signal;

            fatFraction_percent = 70;
            freqWrapStep_Hz = 1/dTE_s;
            fieldmap_Hz = (0:0.01:4) * freqWrapStep_Hz; % compute in range [0, 2pi] (phase_rad pi <=> fieldmap_Hz 1/2/dTE_s)
            waterR2s_Hz = 0:10:300;
            fatR2s_Hz = waterR2s_Hz;
            ModelParamsVary1 = build_ModelParams(F.centerFreq_Hz, F.get_FatModel, fatFraction_percent, fieldmap_Hz, waterR2s_Hz, fatR2s_Hz);
            
            ModelParamsTrue.Species(1).magnitude = 100 - F.fatFraction_percent;
            ModelParamsTrue.Species(1).phase = 0;
            ModelParamsTrue.Species(1).freq_Hz = 0;
            ModelParamsTrue.Species(1).R2s_Hz = F.R2s_Hz;
            ModelParamsTrue.Species(2).magnitude = F.fatFraction_percent;
            ModelParamsTrue.Species(2).phase = 0;
            ModelParamsTrue.Species(2).freq_Hz = ModelParamsTrue.Species(1).freq_Hz - 3.4e-6 * F.centerFreq_Hz;
            ModelParamsTrue.Species(2).R2s_Hz = F.R2s_Hz;
            
            testCase.verifyEqual(F.get_paramsMatrix, ModelParams2paramsMatrix(ModelParamsTrue));
        end
        
        
        function test_set_characteristicParameters(testCase)
            F = FatModel('scanner')
            F.set_characteristicParameters
            F
        end

    end % methods

end % classdef
