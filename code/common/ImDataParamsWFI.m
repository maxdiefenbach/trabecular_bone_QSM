classdef ImDataParamsWFI < ImDataParamsBase
        
    properties
            AlgoParams
            WFIparams
            VARPROparams
        end
        
        methods
            
            %% Constructor
function this = ImDataParamsWFI(filename)
    this = this@ImDataParamsBase(filename);
    this.set_FatModel;
end


function set_FatModel(this, name)
    if nargin < 2
        name = 'Ren marrow';
    end
    F = FatModel(name);
    this.WFIparams.FatModel = F.get_FatModel;
end


%% WFIparams
function load_VARPROparams(this, filename)
    if nargin < 2
        folder = fileparts(this.filename);
        filename = fullfile(folder, [this.fileID '_VARPROparams.mat'])
    end
    this.load_propertyFromFile('VARPROparams', filename);
end

function save_VARPROparams(this, filename)
    assert(~isempty(this.WFIparams))
    if nargin < 2
        filename = fullfile(pwd, [this.fileID '_VARPROparams' '.mat']);
    elseif isdir(filename)
        filename = fullfile(filename, [this.fileID '_VARPROparams' '.mat']);
    end
    this.save_property('VARPROparams', filename);
end

function load_WFIparams(this, filename)
    if nargin < 2
        folder = fileparts(this.filename);
        filename = fullfile(folder, [this.fileID '_WFIparams.mat'])
    end
    this.load_propertyFromFile('WFIparams', filename);
end


function save_WFIparams(this, filename)
    assert(~isempty(this.WFIparams))
    this.WFIparams
    initialization = '';            
    if isfield(this.WFIparams, 'initialization')
        initialization = ['_' this.WFIparams.initialization];
    end
    if nargin < 2
        filename = fullfile(pwd, [this.fileID ...
                            '_WFIparams_' this.WFIparams.method ...
                            initialization '.mat']);
    elseif isdir(filename)
        filename = fullfile(filename, ...
                            [this.fileID ...
                            '_WFIparams_' this.WFIparams.method ...
                            initialization '.mat']);
    end
    this.save_property('WFIparams', filename);
end


function set_WFIparams(this, outParamsFromToolbox)
    FatModel = this.WFIparams.FatModel;
    this.WFIparams = [];
    this.WFIparams.FatModel = FatModel;

    water = outParamsFromToolbox.species(1).amps;
    this.WFIparams.water = water;
    
    fat = outParamsFromToolbox.species(2).amps;
    this.WFIparams.fat = fat;
    
    if isfield(outParamsFromToolbox, 'fieldmap');
        fieldmap_Hz = outParamsFromToolbox.fieldmap;
        this.WFIparams.fieldmap_Hz = fieldmap_Hz;
    end
    
    if isfield(outParamsFromToolbox, 'r2starmap');
        r2starmap_Hz = outParamsFromToolbox.r2starmap;
        this.WFIparams.R2s_Hz = r2starmap_Hz;
    end
    
    if isfield(outParamsFromToolbox, 'phasemap');
        phasemap = outParamsFromToolbox.phasemap;
        this.WFIparams.phasor = phasemap;
    end
    
    if isfield(outParamsFromToolbox, 'phasor');
        phasor = outParamsFromToolbox.phasor;
        this.WFIparams.phasor = phasor;
    end
end


function unwrap_fieldmap(this)
% V. Fortier and I. R. Levesque, 
% Phase Processing for Quantitative Susceptibility Mapping of 
% Regions With Large Susceptibility and Lack of Signal, 
% Magn. Reson. Med., DOI 10.1002/mrm.26989, 2017.
    airSignalThreshold_percent = set_option(this.AlgoParams, ...
                                            'airSignalThreshold_percent', 3);
    tissueMask = this.get_tissueMask(airSignalThreshold_percent);

    dTE_s = diff(this.ImDataParams.TE_s(1:2));
    qualityCutoff = 3.5;
    phase = 2 * pi * this.WFIparams.fieldmap_Hz * dTE_s;
    phase_unwrapped = qualityGuidedUnwrapping(phase, ...
                                              tissueMask, ...
                                              qualityCutoff);
    this.WFIparams.fieldmap_Hz_unwrap = phase_unwrapped / (2 * pi * dTE_s);
end


function equalize_fieldmap_periods(this)
% if slices are in different varpro residual periods
% this function bring every slice of the field map 
% in the period around 0
    this.WFIparams.fieldmap_Hz = equalize_fieldmap_periods(this.WFIparams.fieldmap_Hz, ...
                                                      this.ImDataParams.TE_s);
end


%% Water-Fat Separation
function IDEAL(this)
    ImDataParams = this.ImDataParams;
    AlgoParams = this.AlgoParams;
    AlgoParams.FatModel = this.WFIparams.FatModel;
    if ~isfield(AlgoParams, 'airSignalThreshold_percent')
        AlgoParams.airSignalThreshold_percent = 5;
    end
    if ~isfield(AlgoParams, 'precision')
        AlgoParams.precision = 1e-4;
    end
    if ~isfield(AlgoParams, 'iterMax')
        AlgoParams.iterMax = 80;
    end
    if isfield(AlgoParams, 'iSlice')
        iSlice = AlgoParams.iSlice;
        ImDataParams.signal = ImDataParams.signal(:, :, iSlice, :);
        if isfield(AlgoParams, 'fieldmap_Hz')
            AlgoParams.fieldmap_Hz = AlgoParams.fieldmap_Hz(:, :, iSlice);
        end
    end
    if isfield(AlgoParams, 'iTE')
        iTE = AlgoParams.iTE;
        ImDataParams.signal = ImDataParams.signal(:, :, :, iTE);
        ImDataParams.TE_s = ImDataParams.TE_s(iTE);
    end
    if ~isfield(AlgoParams, 'mask')
        tissueMask = get_tissueMask(ImDataParams.signal, AlgoParams.airSignalThreshold_percent);
        AlgoParams.mask = logical(tissueMask);
    end
    
    ImDataParams.signal = double(ImDataParams.signal);
    
    WFIparams = IDEALmex(ImDataParams, AlgoParams);
    
    this.WFIparams = WFIparams;
    this.AlgoParams = AlgoParams;
end

function CSS(this)
    ModelParams.FatModel = this.WFIparams.FatModel;
    ModelParams.FatModel
    this.WFIparams.initialization = '';
    if isfield(this.WFIparams, 'fieldmap_Hz')
        ModelParams.fieldmap_Hz = this.WFIparams.fieldmap_Hz;
        fprintf('Initialize fieldmap with %s result. Done.\n', this.WFIparams.method);
        this.WFIparams.initialization = this.WFIparams.method;
        this.WFIparams.initFieldmap_Hz = this.WFIparams.fieldmap_Hz;
    end
    WFIparams = separate_waterFat(this.ImDataParams, ModelParams, this.AlgoParams);
    this.WFIparams = merge_Structs(this.WFIparams, WFIparams);
    this.WFIparams.FatModel = ModelParams.FatModel;
    this.WFIparams.method = 'CSS';
end

%% GANDALFs from here
function calculateVARPRO(this)
    imDataParams = this.get_imDataParams4toolbox;
    
    if ~isfield(this.AlgoParams, 'range_r2star')
        this.AlgoParams.range_r2star = [0 500];
        fprintf('choose range of r2star values (default = [0 500], found to be optimal by Cui et al.)\n');
    end
    
    if ~isfield(this.AlgoParams, 'NUM_R2STARS')
        this.AlgoParams.NUM_R2STARS = 30;
        fprintf('choose number of r2star values (default = 30, values >= 30 found to be optimal by Cui et al.)\n');
    end
    
    if ~isfield(this.AlgoParams, 'do_fill_mask3D')
        this.AlgoParams.do_fill_mask3D = 1;
        fprintf('If a closed object has low signal regions inside, they can be in-(1) or ex(0)-cluded from the graph.  (AlgoParams.do_fill_mask3D, default = 1\n');
    end
    
    if ~isfield(this.AlgoParams, 'sampling_stepsize')
        this.AlgoParams.sampling_stepsize = 2; %Hz
        fprintf('No sampling stepsize chosen (AlgoParams.sampling_stepsize, default = 2Hz)\n');
    end
    
    if ~isfield(this.AlgoParams, 'nSamplingPeriods')
        this.AlgoParams.nSamplingPeriods = 1;
        fprintf('Number of periods to be sampled not specified (AlgoParams.nSamplingPeriods, default = 1)\n2 periods should be enough for most applications\n');
    end
    algoParams = this.get_algoParams4toolbox;
    
    start = tic;
    VARPROparams = calculateVARPRO(imDataParams, algoParams);
    elapsedTime_s = toc(start);
    VARPROparams.t_elapsed = elapsedTime_s;
    VARPROparams.do_fill_mask3D = this.AlgoParams.do_fill_mask3D;
    this.VARPROparams = VARPROparams;
    this.save_VARPROparams;
end

function GANDALF2D(this)
    imDataParams = this.get_imDataParams4toolbox;
    
    if isfield(this.AlgoParams, 'iSlice')
        iSlice = this.AlgoParams.iSlice;
        imDataParams.images = imDataParams.images(:, :, iSlice, :, :);
    end
    
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    
    if ~isfield(this.AlgoParams, 'sampling_stepsize')
        this.AlgoParams.sampling_stepsize = 2; %Hz
        fprintf('No sampling stepsize chosen (AlgoParams.sampling_stepsize, default = 2Hz)\n');
    end

    algoParams = this.get_algoParams4toolbox;
    
    start = tic;
    outParams = GANDALF2D(imDataParams, algoParams);
    elapsedTime_s = toc(start);
    
    this.WFIparams.fieldmap_Hz = outParams.fieldmap;
    this.WFIparams.water = outParams.water;
    this.WFIparams.fat = outParams.fat;
    this.WFIparams.r2starmap = outParams.r2starmap;
    this.WFIparams.method = 'GANDALF2D';
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.sampling_stepsize = this.AlgoParams.sampling_stepsize;
    
    this.set_fatFraction_percent;
    this.set_T2s_ms;
end


function GANDALF3D(this)
    imDataParams = this.get_imDataParams4toolbox;
    
    if isfield(this.AlgoParams, 'iSlice')
        iSlice = this.AlgoParams.iSlice;
        imDataParams.images = imDataParams.images(:, :, iSlice, :, :);
    end
    
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    
    if ~isfield(this.AlgoParams, 'sampling_stepsize')
        this.AlgoParams.sampling_stepsize = 2; %Hz
        fprintf('No sampling stepsize chosen (AlgoParams.sampling_stepsize, default = 2Hz)\n');
    end

    algoParams = this.get_algoParams4toolbox;
    
    start = tic;
    outParams = GANDALF3D(imDataParams, algoParams);
    elapsedTime_s = toc(start);
    
    this.WFIparams.fieldmap_Hz = outParams.fieldmap;
    this.WFIparams.water = outParams.water;
    this.WFIparams.fat = outParams.fat;
    this.WFIparams.r2starmap = outParams.r2starmap;
    this.WFIparams.method = 'GANDALF3D';
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.sampling_stepsize = this.AlgoParams.sampling_stepsize;
    
    this.set_fatFraction_percent;
    this.set_T2s_ms;
end

function GANDALF2D_VL(this)
    imDataParams = this.get_imDataParams4toolbox;
    
    if isfield(this.VARPROparams, 'costLocalMinimaRescale')
        imDataParams.VARPROparams = this.VARPROparams;
        this.AlgoParams.VARPROparamsexist = 1;
        this.AlgoParams.nSamplingPeriods = imDataParams.VARPROparams.nSamplingPeriods;
        
        if isfield(this.AlgoParams, 'iSlice')
            iSlice = this.AlgoParams.iSlice;
            imDataParams.VARPROparams.costLocalMinimaRescale = imDataParams.VARPROparams.costLocalMinimaRescale(:, :, iSlice, :);
            imDataParams.VARPROparams.indexLocalMinimaRescale = imDataParams.VARPROparams.indexLocalMinimaRescale(:, :, iSlice, :);
            imDataParams.VARPROparams.masksignal = imDataParams.VARPROparams.masksignal(:, :, iSlice);
            imDataParams.VARPROparams.nMinimaPerVoxel = imDataParams.VARPROparams.nMinimaPerVoxel(:, :, iSlice);
        end
        fprintf('VARPROparams found\n');
        imDataParams.VARPROparams
    end
    
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    
    if ~isfield(this.AlgoParams, 'sampling_stepsize')
        this.AlgoParams.sampling_stepsize = 2; %Hz
        fprintf('No sampling stepsize chosen (AlgoParams.sampling_stepsize, default = 2Hz)\n');
    end
    
    if ~isfield(this.AlgoParams, 'nSamplingPeriods')
        this.AlgoParams.nSamplingPeriods = 1;
        fprintf('Number of periods to be sampled not specified (AlgoParams.nSamplingPeriods, default = 1)\n2 periods should be enough for most applications\n');
    end

    algoParams = this.get_algoParams4toolbox;
    
    start = tic;
    outParams = this.run_singleSliceAlgorithm(@GANDALF, imDataParams, algoParams);
    elapsedTime_s = toc(start);
    
    this.set_WFIparams(outParams);
    this.WFIparams.method = 'GANDALF2D_VL';
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.sampling_stepsize = this.AlgoParams.sampling_stepsize;
    this.WFIparams.nSamplingPeriods = this.AlgoParams.nSamplingPeriods;
    this.set_fatFraction_percent;
    this.set_T2s_ms;
end


function GANDALF3D_VL(this)
    fprintf(2,'Exhaustive amount of RAM is needed for big datasets with many slices.\nGANDALF does not check your machines capabilities, it may abort!\n\n');
    imDataParams = this.get_imDataParams4toolbox;
    imDataParams.fileID = this.fileID;
    
    if isfield(this.AlgoParams, 'iSlice')
        iSlice = this.AlgoParams.iSlice;
        imDataParams.images = imDataParams.images(:, :, iSlice, :, :);
    end
    
    if isfield(this.VARPROparams, 'costLocalMinimaRescale')
        imDataParams.VARPROparams = this.VARPROparams;
        this.AlgoParams.VARPROparamsexist = 1;
        this.AlgoParams.nSamplingPeriods = imDataParams.VARPROparams.nSamplingPeriods;
        
        if isfield(this.AlgoParams, 'iSlice')
            imDataParams.VARPROparams.costLocalMinimaRescale = imDataParams.VARPROparams.costLocalMinimaRescale(:, :, iSlice, :);
            imDataParams.VARPROparams.indexLocalMinimaRescale = imDataParams.VARPROparams.indexLocalMinimaRescale(:, :, iSlice, :);
            imDataParams.VARPROparams.masksignal = imDataParams.VARPROparams.masksignal(:, :, iSlice);
            imDataParams.VARPROparams.nMinimaPerVoxel = imDataParams.VARPROparams.nMinimaPerVoxel(:, :, iSlice);
        end
        fprintf('VARPROparams found\n');
        imDataParams.VARPROparams
    end

    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    
    if ~isfield(this.AlgoParams, 'sampling_stepsize')
        this.AlgoParams.sampling_stepsize = 2; %Hz
        fprintf('No sampling stepsize chosen (AlgoParams.sampling_stepsize, default = 2Hz)\n');
    end
    
    if ~isfield(this.AlgoParams, 'nSamplingPeriods')
        this.AlgoParams.nSamplingPeriods = 1;
        fprintf('Number of periods to be sampled not specified (AlgoParams.nSamplingPeriods, default = 1)\n2 periods should be enough for most applications\n');
    end
    

    
    algoParams = this.get_algoParams4toolbox;
    start = tic;
    outParams = GANDALF(imDataParams, algoParams);
    elapsedTime_s = toc(start);
    
    this.set_WFIparams(outParams);
    this.WFIparams.method = 'GANDALF3D_VL';
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.sampling_stepsize = this.AlgoParams.sampling_stepsize;
    this.WFIparams.nSamplingPeriods = this.AlgoParams.nSamplingPeriods;
    
    this.set_fatFraction_percent;
    this.set_T2s_ms;
end

%% ISMRM Water Fat Toolbox Methods
function imDataParams4toolbox = get_imDataParams4toolbox(this)
    imDataParams4toolbox = get_imDataParams4toolbox(this.ImDataParams);
end


function algoParams4toolbox = get_algoParams4toolbox(this)
% default multipeak fat model is the one used as the default in the
% mDixon quant sequence
% - algoParams.species(1).name = 'water' % Water
% - algoParams.species(1).frequency = [0] 
% - algoParams.species(1).relAmps = [1]   
% - algoParams.species(2).name = 'fat' % Fat
% - algoParams.species(2).frequency = [3.80, 3.40, 2.60, 1.94, 0.39, -0.60]
% - algoParams.species(2).relAmps = [0.087 0.693 0.128 0.004 0.039 0.048]
    algoParams.species(1).name = 'water';
    algoParams.species(1).frequency = 0;
    algoParams.species(1).relAmps = 1;
    algoParams.species(2).name = 'fat (7 peaks)';
% 7 peak fat model with same frequencies and relative amplitudes like in
% mDixon sequence
    algoParams.species(2).frequency = this.WFIparams.FatModel.freqs_ppm;
% [ppm] relative to water peak, one positive peak, rest is negative
    algoParams.species(2).relAmps = this.WFIparams.FatModel.relAmps;
    
    algoParams4toolbox = merge_Structs(algoParams, this.AlgoParams);
end


function hIDEAL(this)
    imDataParams = this.get_imDataParams4toolbox;
    if isfield(this.AlgoParams, 'iSlice')
        iSlice = this.AlgoParams.iSlice;
        imDataParams.images = imDataParams.images(:, :, iSlice, :, :);
    end
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end

    algoParams = this.get_algoParams4toolbox;
    
% Algorithm-specific parameters
    algoParams.Verbose = set_option(this.AlgoParams, 'Verbose', 1);
    algoParams.AlwaysShowGUI = set_option(this.AlgoParams, 'AlwaysShowGUI', 0);
    algoParams.Visualize = set_option(this.AlgoParams, 'Visualize', 0);
    algoParams.Visualize_FatMapMultipler = set_option(this.AlgoParams, 'Visualize_FatMapMultipler', 1.);
    algoParams.CorrectAmpForT2star = set_option(this.AlgoParams, 'CorrectAmpForT2star', 0);
    algoParams.MaxNumDiv = set_option(this.AlgoParams, 'MaxNumDiv', 6);
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
    algoParams            
    
    tic
    outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams);
    elapsedTime_s = toc;

    this.set_WFIparams(outParams);
    phasor = this.WFIparams.phasor;
    dTE_s = imDataParams.TE(2) - imDataParams.TE(1);
    fieldmap_Hz = angle(phasor) ./ (2 * pi * dTE_s);
    this.WFIparams.fieldmap_Hz = fieldmap_Hz;
    this.WFIparams.method = 'hIDEAL';
    if isfield(this.AlgoParams, 'iSlice')
        this.WFIparams.iSlice = this.AlgoParams.iSlice;
    end
    this.set_fatFraction_percent;
    this.set_T2s_ms;
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.algoParams = algoParams;
end


function multiseedRG(this)
% Description: Fat-water separation from three complex echoes with uniform
%              echo time spacing, using a multi-seeded region growing
%              scheme, as described in:
%
% Berglund J, Johansson L, Ahlstr??m H, Kullberg J. Three-point Dixon method enables whole-body 
% water and fat imaging of obese subjects. Magn Reson Med. 2010, Jun;63(6):1659-1668.
    
    imDataParams = this.get_imDataParams4toolbox;
    if isfield(this.AlgoParams, 'iSlice')
        iSlice = this.AlgoParams.iSlice;
        imDataParams.images = imDataParams.images(:, :, iSlice, :, :);
    end
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end

    algoParams = this.get_algoParams4toolbox;

    %%   - algoParams.c1 = 0.75; % Threshold on magnitude weight for seed points
    %%   - algoParams.c2 = 0.25; % Threshold on |log(W/F)| for seed points
    
% run multi-seed region growing
    tic
    outParams = fw_i3cm0i_3point_berglund(imDataParams, algoParams);
    elapsedTime_s = toc;
    
    outParams.fieldmap = -outParams.fieldmap;

% set method string
    this.set_WFIparams(outParams);
    this.WFIparams.method = 'multiseedRG';
    this.set_fatFraction_percent;
    this.set_T2s_ms;
    this.WFIparams.elapsedTime_s = elapsedTime_s;
end


%% single-slice algorithms
function outParams = run_singleSliceAlgorithm(this, functionHandle, imDataParams, algoParams)
    
    matrixSize = size(imDataParams.images);
    matrixSize = matrixSize(1:3);

    if ~isfield(this.AlgoParams, 'iSlice')
        this.AlgoParams.iSlice = 1:matrixSize(3);
    end
    
    nSlices = length(this.AlgoParams.iSlice);
    
    outParams.species(1).amps = complex(zeros([matrixSize(1:2), nSlices]));
    outParams.species(2).amps = complex(zeros([matrixSize(1:2), nSlices]));
    outParams.fieldmap = zeros([matrixSize(1:2), nSlices]);
    outParams.r2starmap = zeros([matrixSize(1:2), nSlices]);
    
    for i = 1:nSlices
        iSlice = this.AlgoParams.iSlice(i);
        fprintf(2,'Processing Slice %d of %d:\n', i, nSlices);
        imDataParamsSlice = imDataParams;
        
        if isfield(imDataParamsSlice.VARPROparams, 'costLocalMinimaRescale')
            imDataParamsSlice.VARPROparams.costLocalMinimaRescale = imDataParamsSlice.VARPROparams.costLocalMinimaRescale(:, :, i, :);
            imDataParamsSlice.VARPROparams.indexLocalMinimaRescale = imDataParamsSlice.VARPROparams.indexLocalMinimaRescale(:, :, i, :);
            imDataParamsSlice.VARPROparams.masksignal = imDataParamsSlice.VARPROparams.masksignal(:, :, i);
            imDataParamsSlice.VARPROparams.nMinimaPerVoxel = imDataParamsSlice.VARPROparams.nMinimaPerVoxel(:, :, i);
        end
        
        imDataParamsSlice.images = imDataParamsSlice.images(:, :, iSlice, :, :);
        algoParamsSlice = algoParams;
        if isfield(algoParamsSlice, 'fieldmap')
            algoParamsSlice.fieldmap = algoParamsSlice.fieldmap(:, :, iSlice);
        end
        
        outParamsSlice = functionHandle(imDataParamsSlice, algoParamsSlice);
        
        outParams.species(1).amps(:, :, i) = outParamsSlice.species(1).amps;
        outParams.species(2).amps(:, :, i) = outParamsSlice.species(2).amps;
        outParams.fieldmap(:, :, i) = outParamsSlice.fieldmap;
        if isfield(outParamsSlice, 'r2starmap')
            outParams.r2starmap(:, :, i) = outParamsSlice.r2starmap;
        end
    end
end


function graphcut(this)
% Description: Fat-water separation using regularized fieldmap
%              formulation and graph cut solution. 
    
% Hernando D, Kellman P, Haldar JP, Liang ZP. Robust water/fat separation in the presence of large 
% field inhomogeneities using a graph cut algorithm. Magn Reson Med. 2010 Jan;63(1):79-90.

    imDataParams = this.get_imDataParams4toolbox;
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    algoParams = this.get_algoParams4toolbox;
    
    dTE_s = diff(imDataParams.TE(1:2));
    maxFM = round(1/dTE_s);
    dFM = 3;
    nFM = set_option(this.AlgoParams, 'NUM_FMS', ceil(2*maxFM/dFM));
    rangeFM = [-maxFM, maxFM];

% Algorithm-specific parameters
    algoParams.size_clique = set_option(this.AlgoParams, 'size_clique', 1); % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
    algoParams.range_r2star = set_option(this.AlgoParams, 'range_r2star', [0 300]);% Range of R2* values
    algoParams.NUM_R2STARS = set_option(this.AlgoParams, 'NUM_R2STARS', 11);% Number of R2* values for quantization
    algoParams.range_fm = set_option(this.AlgoParams, 'rangeFM', rangeFM);% Range of field map values
    algoParams.NUM_FMS = set_option(this.AlgoParams, 'NUM_FMS', nFM);% Number of field map values to discretize
    algoParams.NUM_ITERS = set_option(this.AlgoParams, 'NUM_ITERS', 40);  % Number of graph cut iterations
    algoParams.SUBSAMPLE = set_option(this.AlgoParams, 'SUBSAMPLE', 2);   % Spatial subsampling for field map estimation (for speed)
    algoParams.DO_OT = set_option(this.AlgoParams, 'DO_OT', 0);% 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
    algoParams.LMAP_POWER = set_option(this.AlgoParams, 'LMAP_POWER', 2);  % Spatially-varying regularization (2 gives ~ uniformn resolution)
    algoParams.lambda = set_option(this.AlgoParams, 'lambda', 0.05);   % Regularization parameter
    algoParams.LMAP_EXTRA = set_option(this.AlgoParams, 'LMAP_EXTRA', 0.05);% More smoothing for low-signal regions
    algoParams.TRY_PERIODIC_RESIDUAL = set_option(this.AlgoParams, 'TRY_PERIODIC_RESIDUAL', 0);
    algoParams.THRESHOLD = set_option(this.AlgoParams, 'THRESHOLD', 0.01);
    algoParams
    
    tic
    outParams = this.run_singleSliceAlgorithm(@fw_i2cm1i_3pluspoint_hernando_graphcut, ...
                                              imDataParams, algoParams);
    elapsedTime_s = toc;

    this.set_WFIparams(outParams);
    this.WFIparams.method = 'graphcut';
    this.set_fatFraction_percent;
    this.set_T2s_ms;
    this.WFIparams.elapsedTime_s = elapsedTime_s;
    this.WFIparams.algoParams = algoParams;
end


function goldenSectionSearch(this)
% Description: Fat-water separation using golden section search

    imDataParams = this.get_imDataParams4toolbox;
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    algoParams = this.get_algoParams4toolbox;
    
    tic
    outParams = this.run_singleSliceAlgorithm(@fw_3point_wm_goldSect, ...
                                              imDataParams, algoParams);
    elapsedTime_s = toc;

    this.set_WFIparams(outParams);
    this.WFIparams.method = 'goldenSectionSearch';
    this.set_fatFraction_percent;
    this.set_T2s_ms;
    this.WFIparams.elapsedTime_s = elapsedTime_s;
end


function GOOSE(this)
    imDataParams = this.get_imDataParams4toolbox;
    if isfield(this.AlgoParams, 'iTE')
        iTE = this.AlgoParams.iTE;
        imDataParams.images = imDataParams.images(:, :, :, :, iTE);
        imDataParams.TE = imDataParams.TE(iTE);
    end
    algoParams = this.get_algoParams4toolbox;
    
% Algorithm-specific parameters
% algoParams.gridsize = round(gridspacing);
% algoParams.size_clique = 1; % Size of MRF neighborhood (1 uses an 8-neighborhood, common in 2D)
% algoParams.range_r2star = [0 30]; % Range of R2* values
% algoParams.NUM_R2STARS = 16; % Numbre of R2* values for quantization
% algoParams.range_fm = [-freqRange freqRange]; % Range of field map values
% algoParams.NUM_FMS = Numlayers; % Number of field map values to discretize
% algoParams.NUM_ITERS = 40; % Number of graph cut iterations
% algoParams.SUBSAMPLE = 1; % Spatial subsampling for field map estimation (for speed)
% algoParams.DO_OT = 1; % 0,1 flag to enable optimization transfer descent (final stage of field map estimation)
% algoParams.LMAP_POWER = 2; % Spatially-varying regularization (2 gives ~ uniformn resolution)
% algoParams.lambda = 0.05; % Regularization parameter
% algoParams.LMAP_EXTRA = 0.05; % More smoothing for low-signal regions
% algoParams.TRY_PERIODIC_RESIDUAL = 0;
    
    tic
%             outParams = GOOSE(imDataParams, algoParams);
    outParams = this.run_singleSliceAlgorithm(@GOOSE, imDataParams, algoParams);
    elapsedTime_s = toc;

    this.set_WFIparams(outParams);
    this.WFIparams.method = 'GOOSE';
    this.set_fatFraction_percent;
    this.set_T2s_ms;
    this.WFIparams.elapsedTime_s = elapsedTime_s;
end


%% T2s FF
function set_fatFraction_percent(this)
    this.WFIparams = get_fatFraction_percent(this.WFIparams);
    disp('Set fat fraction [%].');
end


function set_T2s_ms(this)
    if isfield(this.WFIparams, 'R2s_Hz')
        this.WFIparams = get_T2s_ms(this.WFIparams);
        disp('Set T2star [ms].')
    else
        disp('No R2star given => no T2star.')
    end
end


function set_inPhaseImg(this)
    this.WFIparams.inPhaseImg = abs(this.WFIparams.water + this.WFIparams.fat);
end


function set_outOfPhaseImg(this)
    this.WFIparams.outOfPhaseImg = abs(this.WFIparams.water - this.WFIparams.fat);
end


%% plotting
function h = plot_WFIparams(this)
    h = plot_WFIparams(this.WFIparams);
    set(h, 'Name', ['WFIparams ' this.WFIparams.method ' ' this.fileID])
end


function h = plot_WFIparamsSlice(this, iSlice)
    if nargin < 2
        iSlice = ceil((size(this.WFIparams.water, 3) + 1)/2);
    end
    h = plot_WFIparamsSlice(this.WFIparams, iSlice);
end

    end % methods

end % EOF