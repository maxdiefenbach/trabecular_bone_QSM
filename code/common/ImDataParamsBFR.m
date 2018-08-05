classdef ImDataParamsBFR < ImDataParamsWFI
    
    properties
        BFRparams
    end
    
    methods
        
        %% Constructor
        function this = ImDataParamsBFR(filename)
            this = this@ImDataParamsWFI(filename);
        end
        
        
        function save_BFRparams(this, filename)
            if nargin < 2
                filename = [this.fileID '_BFRparams.mat'];
            end
            this.save_property('BFRparams', filename);        
        end
        
        
        function load_BFRparams(this, filename)
            if nargin < 2
                filename = fullfile(fileparts(this.filename), [this.fileID '_BFRparams.mat']);
            end
            this.load_propertyFromFile('BFRparams', filename);
        end


        function set_BFRparams(this)
            if isfield(this.WFIparams, 'voxelSize_mm_interp')
                this.BFRparams.voxelSize_mm = this.WFIparams.voxelSize_mm_interp;
                disp('Use interpolated voxelSize_mm.')
            else
                this.BFRparams.voxelSize_mm = this.ImDataParams.voxelSize_mm;
            end
            this.BFRparams.B0dir = this.ImDataParams.B0dir;
            this.copy_fieldmap2BFRparams;
            this.set_RDF_ppm;
        end


        function set_RDF_ppm(this)      % relative difference field
            if ~isfield(this.BFRparams, 'fieldmap_Hz')
                this.copy_fieldmap2BFRparams;
            end
            fieldmap_Hz = this.BFRparams.fieldmap_Hz;
            centerFreq_Hz = this.ImDataParams.centerFreq_Hz;
            
            RDF_ppm = (fieldmap_Hz ./ centerFreq_Hz) * 1e6;
            
            this.BFRparams.totalRDF_ppm = RDF_ppm;
        end

        
        function copy_fieldmap2BFRparams(this)
            this.BFRparams.fieldmap_Hz = this.WFIparams.fieldmap_Hz_unwrap;
        end
        

        function tissueMask = get_tissueMask(this, airSignalThreshold_percent)
            if nargin < 2
                airSignalThreshold_percent = this.Masks.airSignalThreshold_percent
            end
            tissueMask = get_tissueMask(this.ImDataParams.signal, airSignalThreshold_percent);
        end

        
        function set_BFRmask(this)      % (boundary) mask for background field removal
            airSignalThreshold_percent = this.AlgoParams.airSignalThreshold_percent;
            orientation = this.ImDataParams.orientation;
            signal = this.ImDataParams.signal;
            
            tissueMask = get_tissueMask(signal, airSignalThreshold_percent);
            switch orientation
              case 'cor'
                openDimension = 1;
              case 'sag'
                openDimension = 2;
              case 'tra'
                openDimension = 3;
            end
            tissueMaskFilled = fill_mask3D(tissueMask, openDimension);
            
            this.BFRparams.BFRmask = tissueMaskFilled;
        end
        
        
        function remove_backgroundField_LBV(this, Options)
            if nargin < 2
                Options = this.BFRparams;
            end
            if ~isfield(Options, 'BFRmask')
                Options.BFRmask = this.BFRparams.BFRmask;
            end
            OutParams = remove_backgroundField_LBV(this.BFRparams, Options);
            this.BFRparams = merge_Structs(OutParams, this.BFRparams);
        end


        %% Laplacian Unwrapping
        function do_LaplacianUnwrapping(this)
            voxelSize_mm = this.ImDataParams.voxelSize_mm;
            P = get_phase(this.ImDataParams);
            Punwrap = do_LaplacianUnwrapping(P, voxelSize_mm);
            M = get_magnitude(this.ImDataParams);
            this.ImDataParams.signal = M .* exp(1j .* Punwrap);
        end


        %% background field removal
        function remove_backgroundField_PDF(this, Options)
        % Background field removal using Projection onto Dipole Fields
        %%%% NMR Biomed 2011;24(9):1129-36.
        %%%% MRM 2010;63(1):194-206
            if nargin < 2
                Options = this.BFRparams;
            end
            if ~isfield(Options, 'BFRmask')
                Options.BFRmask = this.BFRparams.BFRmask;
            end
            BFRparams = remove_backgroundField_PDF(this.BFRparams, Options)
            this.BFRparams = merge_Structs(this.BFRparams, BFRparams);
        end


        function remove_backgroundField_PDF_MEDI(this)        
        % Background field removal using Projection onto Dipole Fields
        %%%% NMR Biomed 2011;24(9):1129-36.
        %%%% MRM 2010;63(1):194-206
            RDF_ppm = this.BFRparams.RDF_ppm;
            matrixSize = size(RDF_ppm);
            Nstd = ones(matrixSize);
            BFRmask = this.BFRparams.BFRmask;
            voxelSize_mm = this.ImDataParams.voxelSize_mm;
            B0dir = this.ImDataParams.B0dir;

            localRDF_ppm = PDF(RDF_ppm, Nstd, BFRmask, matrixSize, voxelSize_mm, B0dir);
            
            this.BFRparams.localRDF_ppm = localRDF_ppm;
            this.BFRparams.backgroundFieldRemovalMethod = 'PDF';
        end
        
        
        function set_backgroundRDF_ppm(this)
            this.BFRparams.backgroundRDF_ppm = this.BFRparams.totalRDF_ppm - this.BFRparams.localRDF_ppm;
        end
        
        
        function h = plot_BFRparams(this)
            h = imagine(this.BFRparams.totalRDF_ppm, 'Name', 'total RDF [ppm]', ...
                        this.BFRparams.BFRmask, 'Name', 'BFR mask', ...
                        this.BFRparams.localRDF_ppm, 'Name', 'local RDF [ppm]', 'Window', [-1, 1], ...
                        this.BFRparams.backgroundRDF_ppm, 'Name', 'background RDF [ppm]');
        end

    end % methods

end % EOF