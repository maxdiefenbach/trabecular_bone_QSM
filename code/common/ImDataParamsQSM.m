classdef ImDataParamsQSM < ImDataParamsWFI & ...
        ImDataParamsBFR
    
    properties
        QSMparams
    end
    
    methods
        
        %% Constructor
        function this = ImDataParamsQSM(filename)
            this@ImDataParamsBFR('');
            this = this@ImDataParamsWFI(filename);
            
            this.set_QSMalgoParams;
        end
        
        
        function set_QSMparams(this)
            this.QSMparams = []
            this.set_QSMalgoParams;
            this.copy_localRDF2QSMparams;
            this.QSMparams.voxelSize_mm = this.BFRparams.voxelSize_mm;
            this.QSMparams.B0dir = this.BFRparams.B0dir;
        end
        
        
        function set_QSMalgoParams(this)
            AlgoParams.precision = 1e-2;
            AlgoParams.iterMax = 70;
            AlgoParams.verbose = 1;
            this.QSMparams.AlgoParams = AlgoParams;
        end
        
        
        function copy_localRDF2QSMparams(this)
            this.QSMparams.localRDF_ppm = this.BFRparams.localRDF_ppm;
        end

        
        function set_scaledEchoMIP_dataWeighting(this)
            scaledEchoMIP = get_echoMIP(this.ImDataParams.signal);
            dataWeighting = scale_array2interval(scaledEchoMIP, [0, 1]);
            this.QSMparams.AlgoParams.dataWeighting;
        end

        
        %% QSMparams
        function load_QSMparams(this, filename)
            if nargin < 2
                filenameCell = get_filenameCell(pwd, [this.fileID '_QSMparams*.mat']);
                filename = filenameCell{1};
            end
            tmp = load(filename, 'QSMparams');
            this.QSMparams = tmp.QSMparams;
            fprintf('loaded %s\n', filename);
        end

        
        function save_QSMparams(this, filename)
            assert(~isempty(this.QSMparams))
            method = this.QSMparams.method;
            if nargin < 2
                filename = [this.fileID '_QSMparams_' method '.mat'];
            end
            if isdir(filename)
                filename = fullfile(filename, ...
                                    [this.fileID '_QSMparams_' method '.mat']);
            end
            this.save_property('QSMparams', filename)
        end


        function closedFormL2(this)
            DataParams.RDF_ppm = this.QSMparams.localRDF_ppm;
            DataParams.voxelSize_mm = this.QSMparams.voxelSize_mm;
            DataParams.B0dir = this.QSMparams.B0dir;
            
            AlgoParams = this.QSMparams.AlgoParams

            QSMparams = closedFormL2(DataParams, AlgoParams);
            
            this.QSMparams = merge_Structs(this.QSMparams, QSMparams)
        end


        function MEDI_l2(this)
            DataParams.RDF_ppm = this.QSMparams.localRDF_ppm;
            DataParams.voxelSize_mm = this.QSMparams.voxelSize_mm;
            DataParams.B0dir = this.QSMparams.B0dir;
            
            AlgoParams = this.QSMparams.AlgoParams

            QSMparams = MEDI_l2_imagespace(DataParams, AlgoParams)
            
            this.QSMparams = merge_Structs(this.QSMparams, QSMparams)
        end
        
        
        function MEDI_l1(this)
            DataParams.RDF_ppm = this.QSMparams.localRDF_ppm;
            DataParams.voxelSize_mm = this.QSMparams.voxelSize_mm;
            DataParams.B0dir = this.QSMparams.B0dir;
            
            AlgoParams = this.QSMparams.AlgoParams

            QSMparams = MEDI_l1(DataParams, AlgoParams)
            
            this.QSMparams = merge_Structs(this.QSMparams, QSMparams)
        end


        function set_dataWeighting_BFRmask(this)
            this.QSMparams.AlgoParams.dataWeighting = this.BFRparams.BFRmask;;
        end

        
        function set_dataWeighting_scaledEchoMIP(this)
            echoMIP = get_echoMIP(this.ImDataParams.signal);
            scaledEchoMIP = scale_array2interval(echoMIP, [0, 1]);
            
            this.QSMparams.AlgoParams.dataWeighting = scaledEchoMIP;
            this.QSMparams.AlgoParams.dataWeightingStr = 'scaledEchoMIP';
        end


        function set_gradWeighting_edgePenalty(this)
            tissueMaskFilled = this.Masks.tissueMaskFilled;
            edgeMask = this.Masks.edgeMask;
            objectBoundary = this.Masks.objectBoundary;
            
            edgePenalty = ~(edgeMask | objectBoundary);
%             edgePenalty = tissueMaskFilled .* ~edgeMask | objectBoundary;

            this.QSMparams.AlgoParams.gradWeighting = edgePenalty;
            this.QSMparams.AlgoParams.gradWeightingStr = 'edgePenalty';
%             this.QSMparams.AlgoParams.gradWeighting = tissueMaskFilled .* edgePenalty;
        end

        
        function set_gradWeighting_subcutaneousFatEnhancement(this)
            tissueMask = this.Masks.tissueMask;
            tissueMaskFilled = this.Masks.tissueMaskFilled;
            subcutaneousFatMask = this.Masks.subcutaneousFatMask;
            edgeMask = this.Masks.edgeMask;
            objectBoundary = this.Masks.objectBoundary;
            edgePenalty = ~(edgeMask | objectBoundary);
            
            subcutaneousFatEnhancement = edgePenalty .* (subcutaneousFatMask + tissueMask) + ~tissueMaskFilled;
            
            this.QSMparams.AlgoParams.gradWeighting = subcutaneousFatEnhancement;
            this.QSMparams.AlgoParams.gradWeightingStr = 'subcutaneousFatEnhancement';
        end

        
        function set_gradWeighting_scMIPsubFatEnh(this)
            echoMIP = get_echoMIP(this.ImDataParams.signal);
            scaledEchoMIP = scale_array2interval(echoMIP, [0, 1]);

            tissueMask = this.Masks.tissueMask;
            tissueMaskFilled = this.Masks.tissueMaskFilled;
            subcutaneousFatMask = this.Masks.subcutaneousFatMask;
            edgeMask = this.Masks.edgeMask;
            objectBoundary = this.Masks.objectBoundary;
            edgePenalty = ~(edgeMask | objectBoundary);
            
            subcutaneousFatEnhancement = edgePenalty .* (subcutaneousFatMask + tissueMask);
            
            this.QSMparams.AlgoParams.gradWeighting = scaledEchoMIP .* subcutaneousFatEnhancement + ~tissueMaskFilled;
            this.QSMparams.AlgoParams.gradWeightingStr = 'scMIPsubFatEnh';
        end
        
        
        function set_initialChimap(this)
            tissueMaskFilled = this.Masks.tissueMaskFilled;
            boneMarrowMask = this.Masks.boneMarrowMask;
            FF = this.OutParams.fatFraction_percent / 100;
            CHI_H2O_ppm = -0.61;
            CHI_bone_ppm = -0.8;
            chiInit_ppm = (1 - FF) .* CHI_H2O_ppm;
%             chiInit_ppm(~tissueMaskFilled) = 0;
            chiInit_ppm(~tissueMaskFilled) = mean(chiInit_ppm(tissueMaskFilled));
            chiInit_ppm(boneMarrowMask) = CHI_bone_ppm;
            
            this.QSMparams.AlgoParams.x0 = chiInit_ppm;
        end
        

        function reference_toMask(this, mask)
            if nargin < 2
                mask = this.Masks.subcutaneousFatMask;
            end
            chimap_ppm = this.QSMparams.chimap_ppm;
            reference_ppm = mean(chimap_ppm(mask));
            chimapRelative_ppm = chimap_ppm - reference_ppm;
            
            this.QSMparams.chimapRelative_ppm = chimapRelative_ppm;
            this.QSMparams.reference_ppm = reference_ppm;
        end

        
        function analyse_ROIs(this)
            this.reference_toMask(this.Masks.subcutaneousFatMask);

            chi = this.QSMparams.chimapRelative_ppm;
            ROIs.meanFat_ppm = mean(chi(this.Masks.subcutaneousFatMask));
            ROIs.stdFat_ppm = std(chi(this.Masks.subcutaneousFatMask));
            ROIs.meanWater_ppm = mean(chi(this.Masks.waterMask));
            ROIs.stdWater_ppm = std(chi(this.Masks.waterMask));
            ROIs.meanBone_ppm = mean(chi(this.Masks.boneMarrowMask));
            ROIs.stdBone_ppm = std(chi(this.Masks.boneMarrowMask));
            ROIs

            this.QSMparams.ROIs = ROIs;
        end
        
        
        function h = plot_QSMparams(this, mask)
            if nargin < 2
                mask = ones(size(this.QSMparams.chimap_ppm));
            end
            if isfield(this.QSMparams.AlgoParams, 'dataWeighting') & isfield(this.QSMparams.AlgoParams, 'gradWeighting')
                h = imagine(mask .* this.QSMparams.localRDF_ppm, 'Name', 'local RDF [ppm]', 'Window', [-2, 2], ...
                            mask .* this.QSMparams.chimap_ppm, 'Name', 'chi [ppm]', 'Window', [-2, 2], ...
                            mask .* this.QSMparams.AlgoParams.dataWeighting(:, :, :, 1), 'Name', 'data weighting (x)',  ...
                            mask .* this.QSMparams.AlgoParams.gradWeighting(:, :, :, 1), 'Name', 'grad weighting (x)');
            else
                h = imagine(mask .* this.QSMparams.localRDF_ppm, 'Window', [-2, 2], 'Name', 'local RDF [ppm]', ...
                            mask .* this.QSMparams.chimap_ppm, 'Window', [-2, 2]), 'Name', 'chi [ppm]';
            end
        end

    end % methods

    end % EOF