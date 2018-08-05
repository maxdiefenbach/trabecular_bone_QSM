classdef ImDataParamsBase < MyHandle
    
    properties
        filename
        ImDataParams
    end
    
    methods
        
        %% Constructor
        function this = ImDataParamsBase(filename)
            this = this@MyHandle;
            if ~isempty(filename)
                this.filename = filename;
                this.load_ImDataParams;
                this.set_fileID;
            end
        end
        
        
        function set_fileID(this)
            [~, fname, ~] = fileparts(this.filename);
            this.fileID = fname(1:20);
        end

        
        function load_ImDataParams(this, filename)
            if nargin < 2
                filename = this.filename;
            end
            if isdir(filename)
                filenameCell = regexpdir(filename, '^*ImDataParams.mat');
                filename = filenameCell{1};
                this.filename = filename;
                this.set_fileID;
            end
            this.load_propertyFromFile('ImDataParams', filename);
        end
        
        
        function save_ImDataParams(this, filename)
            if nargin < 2
                filename = [this.fileID '_ImDataParams.mat'];
            end
            this.save_property('ImDataParams', filename);
        end

        
        function h = plot(this, varargin)
            S = this.ImDataParams.signal;
            nTE = size(S, ndims(S));
            h = imagine(S, varargin{:});
            set(h, 'Name', ['signal ' this.fileID]);
            if nTE > 16
                h = imagine(S(:, :, :, 17:end));
                set(h, 'Name', ['signal ' this.fileID]);
            end
        end


        function set_B0dir(this)
            B0dir = get_B0dir(this.ImDataParams);
            this.ImDataParams.B0dir = B0dir;
        end

        
        function set_matrixSize(this)
            this.ImDataParams.matrixSize = get_matrixSize(this.ImDataParams);
        end
        
        
        function matrixSize = get_matrixSize(this)
            matrixSize = get_matrixSize(this.ImDataParams);
        end

        
        function M = magnitude(this)
            M = get_magnitude(this.ImDataParams);
        end


        function P = phase(this)
            P = get_phase(this.ImDataParams);
        end

        
        function load_field(this, fileName, fieldName)
            tmp = load(fileName);
            this.(fieldName) = tmp.(fieldName);
        end


        %% save
        function save_array2nii(this, array, filename, voxelSize, isoCenter, datatype, description)
            if nargin < 7 | isempty(description)
                description = this.fileID;
            end
            if nargin < 6 | isempty(datatype)
                datatype = [];
            end
            if nargin < 5 | isempty(isoCenter)
                isoCenter = this.ImDataParams.isoCenter_REC;
            end
            if nargin < 4 | isempty(voxelSize)
                voxelSize = this.ImDataParams.voxelSize_mm;
            end
            if nargin < 3 | isempty(filename)
                filename = fullfile(fileparts(this.filename), [this.fileID '_' inputname(2) '.nii.gz']);
            end
            nii = make_nii(real(array), voxelSize, isoCenter, datatype, description);
            dicom2nifti = get_permMat_convention('dicom', 'nifti');
            mat2RAI = get_permMat_mat2RAI(this.ImDataParams.orientation);
        % use orientation method 3 of nifti standard
        % give affine transform matrix of MATLAB's coordinates mat to nifti's xyz
            T = dicom2nifti * mat2RAI;
            nii.hdr.hist.sform_code = 1; 
            nii.hdr.hist.srow_x = [T(1, :) 0];
            nii.hdr.hist.srow_y = [T(2, :) 0];
            nii.hdr.hist.srow_z = [T(3, :) 0];
            save_nii(nii, filename);
            fprintf('Wrote %s.\n', GetFullPath(filename));
        end
        
        
        function interp(this, factors)
            this.ImDataParams.signal = interp_array3D(this.ImDataParams.signal, factors);
            this.ImDataParams.voxelSize_mm = this.ImDataParams.voxelSize_mm ./ factors;
        end
        
        
        function Prop = get_interp_property(this, propName, factors)
            prop = this.get_property(propName);
            Prop = interp_array3D(prop, factors);
        end

        
        function set_interp_property(this, propName, factors)
            prop = this.get_property(propName);
            Prop = interp_array3D(prop, factors);
            propNameCell = strsplit(propName, '.');
            nLevels = numel(propNameCell);
            propChain = '';
            for iLevel = 1:nLevels
                propChain = [propChain sprintf('.(propNameCell{%d})', iLevel)];
            end
            cmd = ['this' propChain ' = Prop;'];
            eval(cmd);
        end

        
        function interp_propertyStruct(this, structName, factors)
            Struct = this.(structName)
            fieldNameCell = fieldnames(Struct);
            for iField = 1:numel(fieldNameCell)
                fieldName = fieldNameCell{iField};
                field = Struct.(fieldName);
                if ndims(field) == 3 | ndims(field) == 4
                    Struct.(fieldName) = interp_array3D(field, factors);
                    fprintf('Interpolate %s by factors [%d, %d, %d].\n', fieldName, factors(1), factors(2), factors(3));
                end
            end
            Struct.voxelSize_mm = this.ImDataParams.voxelSize_mm;
            Struct.voxelSize_mm_interp = this.ImDataParams.voxelSize_mm ./ factors;
            this.(structName) = Struct;
        end
        
        
        function echoMIP = get_echoMIP(this)
            echoMIP = get_echoMIP(this.ImDataParams.signal);
        end

    end % methods

    end % EOF