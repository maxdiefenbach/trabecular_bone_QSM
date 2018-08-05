classdef fwtoolboxTest < matlab.unittest.TestCase
    
    properties
        imDataParams
        algoParams
    end
    
    methods (TestClassSetup)
        function add2path(testCase)
            p = path;
            testCase.addTeardown(@path, p);
            addpath(genpath('../../'));
        end
        function set_imDataParams(testCase)
            testfile = '../../testdata/fwtoolbox_v1_data/kellman_data/PKdata1.mat';
            tmp = load(testfile)
            imDataParams = tmp.data;
            testCase.imDataParams = imDataParams;
        end
        function set_algoParams(testCase)
            algoParams.species(1).name = 'water';
            algoParams.species(1).frequency = 0;
            algoParams.species(1).relAmps = 1;
            algoParams.species(2).name = 'fat (7 peaks)';
            algoParams.species(2).frequency = -[3.30, 2.57, -0.71, 3.70, 3.01, 2.35, 1.83];
            algoParams.species(2).relAmps = [0.625, 0.095, 0.042, 0.085, 0.071, 0.066, 0.016];
            testCase.algoParams = algoParams;
        end
    end % methods (TestClassSetup)



    methods (Test)
        
        function test_imDataParams(testCase)
            testCase.imDataParams
        end
        
        
        function test_hIDEAL(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            % Algorithm-specific parameters
            algoParams.Verbose = 1;
            algoParams.AlwaysShowGUI = 0;
            algoParams.Visualize = 0;
            algoParams.Visualize_FatMapMultipler = 1.;
            algoParams.CorrectAmpForT2star = 0;
            %   - algoParams.Verbose - (optional) 1 = show info, 0 = no info
            %   - algoParams.AlwaysShowGUI - (optional) 1 = always show GUI to verify input
            %              parameters, 0 = no GUI if not needed. Default: 1
            %   - algoParams.Visualize - (optional) 1 = show graphics, 0 = no graphics
            %   - algoParams.Visualize_FatMapMultipler - (optional) Multiply fat fraction
            %              by this value for visualization (default 1.5)
            %   - algoParams.MinFractSizeToDivide - (optional) Minimum size of the image
            %              to stop dividing into smaller regions. (default 0.05)
            %   - algoParams.MaxNumDiv - (optional) Maximum number of hierarchical
            %              sub-divisions (default 6)
            %   - algoParams.AssumeSinglePeakAsWater - (optional) 1 = if there is only
            %              one peak, assume that it is water. (default 0)
            %   - algoParams.SnrToAssumeSinglePeak - (optional) Minimum signal to
            %              noise to assess when it is just noise. (default 2.5)
            %   - algoParams.CorrectAmpForT2star - (optional) 1 = provide water and fat
            %       maps that have been corrected for T2* decay. (default 1)
            %   - algoParams.MaxR2star - (optional) maximum acceptable R2* in case of 
            %     erroneously large values from bad fit (default 250 s^-1)

            outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams);
        end
        
        
        function test_graphcut(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;

            %   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
            %   - algoParams.range_r2star = [0 0]; % Range of R2* values
            %   - algoParams.NUM_R2STARS = 1; % Numbre of R2* values for quantization
            %   - algoParams.range_fm = [-400 400]; % Range of field map values
            %   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
            %   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
            %   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
            %   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
            %   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
            %   - algoParams.lambda = 0.05; % Regularization parameter
            %   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
            %   - algoParams.TRY_PERIODIC_RESIDUAL = 0; % Take advantage of periodic residual if uniform TEs (will change range_fm)  
            %   - algoParams.residual: in case we pre-computed the fit residual (mostly for testing) 

            outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams);            
        end

        
        function test_mixed_fitting(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams);
            algoParams.fieldmap = outParams.fieldmap;
            algoParams.r2starmap = outParams.r2starmap;
            
            %   - algoParams.range_r2 = [0 0]; % Range of R2* values           
            %   - algoParams.range_fm = [-400 400]; % Range of field map values
            %   - algoParams.NUM_ITERS = 40; % Number of descent iterations
            %   - algoParams.fieldmap = []; % initial fieldmap
            %   - algoParams.r2starmap = []; % initial fieldmap
            
            % Always add these default values.
            algoParams.NUM_MAGN = 1;
            algoParams.THRESHOLD = 0.04;
            algoParams.range_r2star = [0 200];

            outParams = fw_i2xm1c_3pluspoint_hernando_mixedfit(imDataParams, algoParams);            
        end


        function test_multiSeedRG(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            %   - algoParams.c1 = 0.75; % Threshold on magnitude weight for seed points
            %   - algoParams.c2 = 0.25; % Threshold on |log(W/F)| for seed points
            
            outParams = fw_i3cm0i_3point_berglund(imDataParams, algoParams);
        end
        
        
        function test_analytical_2point_method(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            outParams = fw_i3cm0c_2flexiblepoint_berglund(imDataParams, algoParams);
        end


        function test_Bydder(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            %   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
            %   - algoParams.range_r2star = [0 100]; % Range of R2* values
            %   - algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
            %   - algoParams.range_fm = [-400 400]; % Range of field map values
            %   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
            %   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
            %   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
            %   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
            %   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
            %   - algoParams.lambda = 0.05; % Regularization parameter
            %   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
            %   - algoParams.TRY_PERIODIC_RESIDUAL = 0;
            
            outParams = fw_i2cm1i_3pluspoint_hernando_graphcut(imDataParams, algoParams);
            
            outParamsCommonPhase = decomposeGivenFieldMapAndDampings_commonPhase( imDataParams,algoParams,outParams.fieldmap,outParams.r2starmap,outParams.r2starmap );
        end
        
        
        function test_goldeSectionSearch(testCase)
            imDataParams = testCase.imDataParams;
            algoParams = testCase.algoParams;
            
            %   - algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
            %   - algoParams.range_r2star = [0 120]; % Range of R2* values
            %   - algoParams.NUM_R2STARS = 11; % Numbre of R2* values for quantization
            %   - algoParams.range_fm = [-400 400]; % Range of field map values
            %   - algoParams.NUM_FMS = 301; % Number of field map values to discretize
            %   - algoParams.NUM_ITERS = 40; % Number of graph cut iterations
            %   - algoParams.SUBSAMPLE = 2; % Spatial subsampling for field map estimation (for speed)
            %   - algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
            %   - algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
            %   - algoParams.lambda = 0.05; % Regularization parameter
            %   - algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
            %   - algoParams.TRY_PERIODIC_RESIDUAL = 0;
            
            outParams = fw_3point_wm_goldSect(imDataParams, algoParams);
        end

    end % methods (Test)

end % EOF
