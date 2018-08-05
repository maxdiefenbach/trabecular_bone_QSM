classdef FatModel < matlab.mixin.Copyable
% class to structure different fat models by number of peaks, 
% their locations and their amplitudes
    
    properties
        modelCSV
        modelTable
        name
        centerFreq_Hz
        relAmps_percent
        freqs_ppm
        freqs_Hz
        ModelParams
        paramsMatrixArray
        Characteristics
        waterPeakLocation_ppm
        nFatPeaks
        TE_s
        fieldStrength_T
        fatFraction_percent
        phase_rad
        fieldmap_Hz
        R2s_Hz
        signal
    end

    
    methods
        
        %% constructor
        function this = FatModel(name)
            this.fieldStrength_T = 3;
            this.set_centerFreq_Hz;
            this.fatFraction_percent = 100;
            this.waterPeakLocation_ppm = 4.7;
            this.fieldmap_Hz = 0;
            this.R2s_Hz = 0;
            this.phase_rad = 0;
            this.modelCSV = fullfile(...
                fileparts([mfilename('fullpath') '.m']), ...
                'fatModel.csv');
            this.TE_s = (0:0.01:10) * 1e-3;
            this.signal = complex(zeros(size(this.TE_s)));
            this.paramsMatrixArray = zeros(11+1, 4);
            if nargin == 1 & ~isempty(name)
                this.set_FatModel(name);
            end
        end

        
        function set_FatModel(this, name)
            this.name = name;
            T = readtable(this.modelCSV);
            this.modelTable = T;
            chemshiftRow = strcmp(T.modelName, name) & strcmp(T.parameterName, 'chemical shift');
            relAmpsRow = strcmp(T.modelName, name) & strcmp(T.parameterName, 'relative amplitude');
            abc2j = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
            freqs_ppm = cellfun(@str2num, T{chemshiftRow, abc2j});
            relAmps_percent = cellfun(@str2num, T{relAmpsRow, abc2j});
            this.freqs_ppm = freqs_ppm(:).';
            this.freqs_Hz = this.centerFreq_Hz * 1e-6 * freqs_ppm;
            this.relAmps_percent = relAmps_percent(:).';
            this.normalize_amplitudes;
            this.nFatPeaks = length(relAmps_percent(relAmps_percent ~= 0));
            this.set_ModelParams;
            this.set_Characteristics;
        end
        
        
        function normalize_amplitudes(this)
            this.relAmps_percent = 100 .* this.relAmps_percent / sum(this.relAmps_percent);
        end


        
        function set_centerFreq_Hz(this, centerFreq_Hz)
            GYRO = 42.58e6;
            if nargin == 2
                this.centerFreq_Hz = centerFreq_Hz;
            else
                this.centerFreq_Hz = GYRO * this.fieldStrength_T;
            end
        end
        
        
        function [Crho, Cphi, Cf, Cr] = get_Cmats(this)
            peakNonzero = this.relAmps_percent ~= 0;
            relAmps_percent = this.relAmps_percent(peakNonzero);
            M = length(relAmps_percent) + 1;
            Crho = eye(M, 2);
            Crho(2:end, 2) = relAmps_percent / 100;
            Cphi = double(logical(Crho));
            Cf = ones(M, 1);
            Cr = ones(M, 1);
        end

        
        function set_ModelParams(this, fatFraction_percent)
            if nargin < 2
                fatFraction_percent = this.fatFraction_percent;
            else
                this.fatFraction_percent = fatFraction_percent;
            end
            ModelParams = struct();
            ModelParams.Species(1).name = 'water';
            ModelParams.Species(1).magnitude = 100 - fatFraction_percent;
            ModelParams.Species(1).phase = this.phase_rad;
            ModelParams.Species(1).chemShift_Hz = 0;
            ModelParams.Species(1).chemShift_ppm = 0;
            ModelParams.Species(1).freq_Hz = this.fieldmap_Hz;
            ModelParams.Species(1).freq_ppm = 1e6 * this.fieldmap_Hz / this.centerFreq_Hz;
            ModelParams.Species(1).R2s_Hz = this.R2s_Hz;
            peakNonzero = this.relAmps_percent ~= 0;
            relAmps_percent = this.relAmps_percent(peakNonzero);
            freqs_ppm = this.freqs_ppm(peakNonzero);
            abc2j = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
            abc2j = abc2j(peakNonzero);
            for ip = 1:length(relAmps_percent)
                ModelParams.Species(ip+1).name = ['Fat peak ' abc2j{ip}];
                ModelParams.Species(ip+1).magnitude = relAmps_percent(ip)/100 * this.fatFraction_percent;
                ModelParams.Species(ip+1).phase = this.phase_rad;
                ModelParams.Species(ip+1).chemShift_ppm = (freqs_ppm(ip) - this.waterPeakLocation_ppm);
                ModelParams.Species(ip+1).chemShift_Hz = (freqs_ppm(ip) - this.waterPeakLocation_ppm) * 1e-6 * this.centerFreq_Hz;
                ModelParams.Species(ip+1).freq_ppm = 1e6 * this.fieldmap_Hz / this.centerFreq_Hz + (freqs_ppm(ip) - this.waterPeakLocation_ppm);
                ModelParams.Species(ip+1).freq_Hz = this.fieldmap_Hz + (freqs_ppm(ip) - this.waterPeakLocation_ppm) * 1e-6 * this.centerFreq_Hz;
                ModelParams.Species(ip+1).R2s_Hz = this.R2s_Hz;
            end
            [Crho, Cphi, Cf, Cr] = this.get_Cmats;
            ModelParams.Crho = Crho;
            ModelParams.Cphi = Cphi;
            ModelParams.Cf = Cf;
            ModelParams.Cr = Cr;
            this.ModelParams = ModelParams;
        end
        
        
        function set_Characteristics(this)
        % Berglund, J., Ahlström, H., & Kullberg, J. (2012). 
        % Model-based mapping of fat unsaturation and chain length 
        % by chemical shift imaging-phantom validation and in vivo feasibility. 
        % Magnetic Resonance in Medicine, 68(6), 1815–1827. http://doi.org/10.1002/mrm.24196
            
            ABCDEFGHIJ = this.freqs_ppm;
            F1 = 9 * ABCDEFGHIJ(1) + ...
                 6 * ABCDEFGHIJ(3) + ...
                 6 * ABCDEFGHIJ(5) + ...
                 2 * ABCDEFGHIJ(7) + ...
                 2 * ABCDEFGHIJ(8) + ...
                 ABCDEFGHIJ(9);
            F2 = 2 * ABCDEFGHIJ(2);
            F3 = 4 * ABCDEFGHIJ(4) + 2 * ABCDEFGHIJ(10);
            F4 = 2 * ABCDEFGHIJ(6) + 2 * ABCDEFGHIJ(10);
            
            CL = 4 + abs((F2 + 4 * F3 + 3 * F4) / (3 * F1));
            UD = abs((F3 + F4) / F1);
            PUD = abs(F4 / F1);
            UF = abs(F3 / (3 * F1));
            SF = 1 - UF;
            PUF = PUD / 3;
            MUF = UF - PUF;
            
            this.Characteristics.Fs = [F1, F2, F3, F4];
            this.Characteristics.chainLength = CL;
            this.Characteristics.unsaturationDegree = UD;
            this.Characteristics.polyUnsaturationDegree = PUD;
            this.Characteristics.unsaturationFraction = UF;
            this.Characteristics.saturationFraction = SF;
            this.Characteristics.polyUnsaturationFraction = PUF;
            this.Characteristics.monoUnsaturationFraction = MUF;
        end

        
        function compute_FatModel(this, cl, ndb, nmidb, spectrumName)
            if nargin < 5
                spectrumName = 'Berglund 10 peaks';
            end
            switch spectrumName
              case 'Berglund 10 peaks'
                ABC2J = [9, 
                         6 * (cl - 4) - 8 * ndb + 2 * nmidb, 
                         6, 
                         4 * (ndb - nmidb), 
                         6, 
                         2 * nmidb, 
                         2, 
                         2, 
                         1, 
                         2 * ndb];
              case 'Hamilton 9 peaks'
                ABC2J = [9, 
                         6 * (cl - 4) - 8 * ndb + 2 * nmidb, 
                         6, 
                         4 * (ndb - nmidb), 
                         6, 
                         2 * nmidb, 
                         4, 
                         0, 
                         1, 
                         2 * ndb];
              case 'Peterson 8 peaks'
                ABC2J = [9, 
                         6 * (cl - 4) - 8 * ndb + 2 * nmidb, 
                         6, 
                         4 * (ndb - nmidb), 
                         6, 
                         2 * nmidb, 
                         4, 
                         0, 
                         2 * ndb + 1, 
                         0];
              case 'Hamilton 6 peaks'
                ABC2J = [9, 
                         6 * (cl - 4) - 8 * ndb + 2 * nmidb + 6, 
                         0, 
                         4 * (ndb - nmidb) + 6, 
                         0, 
                         2 * nmidb, 
                         4, 
                         0, 
                         0,
                         2 * ndb + 1];
            end
            relAmps_percent = 100 * ABC2J / sum(ABC2J);
            
            this.name = 'computed';
            T = readtable(this.modelCSV);
            this.modelTable = T;
            chemshiftRow = strcmp(T.modelName, spectrumName) & strcmp(T.parameterName, 'chemical shift');
            abc2j = {'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J'};
            this.freqs_ppm = cellfun(@str2num, T{chemshiftRow, abc2j});
            this.relAmps_percent = relAmps_percent;
            this.nFatPeaks = sum(relAmps_percent > 0);
            this.set_ModelParams;
            this.set_Characteristics;
            this.Characteristics.spectrumName = spectrumName;
            this.Characteristics.chainLength = cl;
            this.Characteristics.numberOfDoubleBonds = ndb;
            this.Characteristics.numberOfMethyleneInterruptedDoubleBonds = nmidb;
        end

        
        function pmN = get_paramsMatrixArray(this)
            pmN = get_paramsMatrixArray(this.ModelParams);
        end

        
        function pm = get_paramsMatrix(this, sub)
            if nargin < 2
                sub = ones(this.nFatPeaks+1, 4);
            end
            pm = get_paramsMatrix(this.ModelParams, sub);
        end


        function signal = build_signal(this, TE_s)
            if nargin ~= 2
                TE_s = this.TE_s;
            end
            sub = ones(this.nFatPeaks+1, 4);
            pm = get_paramsMatrix(this.ModelParams, sub);
            signal = build_signal(TE_s, pm);
            this.signal = signal;
        end
        
        
        function h = plot_signal(this)
            h = figure;
            hold on
            subplot(1, 2, 1)
            p = plot(this.TE_s * 1e3, abs(this.signal));
            set(p, 'LineWidth', 2);
            grid on;
            xlabel('t / [ms]')
            ylabel('magnitude [a.u.]')
            
            subplot(1, 2, 2)
            p = plot(this.TE_s * 1e3, angle(this.signal)/pi);
            set(p, 'LineWidth', 2);
            grid on;
            hold on
            p = plot(this.TE_s * 1e3, unwrap(angle(this.signal))/pi);
            set(p, 'LineWidth', 2);
            xlabel('t / [ms]')
            ylabel('phase [\pi]')
        end


        function h = plot_spectrum(this)
            FatSpectrum = this.get_FatSpectrum;
            h = plot_FatModel(FatSpectrum);
            grid on
            title(this.name);
        end
        

        function FatSpectrum = get_FatSpectrum(this)
            nonzero = this.freqs_ppm ~= 0;
            freqs_ppm = this.freqs_ppm(nonzero) - this.waterPeakLocation_ppm;
            relAmps = this.relAmps_percent(nonzero) / sum(this.relAmps_percent);
            
            FatSpectrum.name = this.name;
            FatSpectrum.freqs_ppm = freqs_ppm(:).';
            FatSpectrum.relAmps = relAmps(:).';
        end
        
        
        function FatModel = get_FatModel(this)
        % for backward compatibility
            FatModel = this.get_FatSpectrum;
        end

        
        function h = plot_speciesVectors(this, TE_s)
            if nargin < 2
                TE_s = (1:6) * 1e-3;
            end
            nFatPeaks = this.nFatPeaks;
            nTE = length(TE_s);

            h = figure;
            set(h, 'Position', [100 1000 1000 2400])
            colors = get(0, 'DefaultAxesColorOrder');
            colors = cat(1, colors, colors);
            colorCell = num2cell(colors, 2);
            for iTE = 1:nTE;
                for iPeak = 1:nFatPeaks
                    subplot(1, nTE, iTE)
                    spinVector = 100 * this.relAmps_percent(iPeak) .* ...
                        exp(2j * pi * (this.fieldmap_Hz + this.freqs_ppm(iPeak) * 1e-6 * this.centerFreq_Hz) * ...
                            TE_s(iTE));
                    p = compass(spinVector);
                    hold on;
                    set(p, 'LineWidth', 2, 'Color', colorCell{iPeak});
                end
            end
        end
        
    end % methods

end % classdef
