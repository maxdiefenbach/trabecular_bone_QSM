classdef ImDataParamsBMRR < ImDataParamsBase & ...
        ImDataParamsWFI & ...
        ImDataParamsBFR & ...
        ImDataParamsQSM

        
    methods
        
        %% Constructor
        function this = ImDataParamsBMRR(filename)
            this = this@ImDataParamsBase(filename);
            this@ImDataParamsWFI('');
            this@ImDataParamsBFR('');
            this@ImDataParamsQSM('');
        end
        
        
        function DataParams = get_DataParams(this)
            fieldmap_Hz = this.WFIparams.fieldmap_Hz;
            centerFreq_Hz = this.ImDataParams.centerFreq_Hz;
            
            DataParams.RDF_ppm = (fieldmap_Hz ./ centerFreq_Hz) * 1e6;
            DataParams.voxelSize_mm = this.ImDataParams.voxelSize_mm;
            DataParams.B0dir = this.ImDataParams.B0dir;            
        end

    end % methods

end % EOF