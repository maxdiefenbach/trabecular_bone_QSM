[pathname n e] = fileparts(which(mfilename));   % gives path of this mfile
BASEFOLDER = [pathname filesep];                % absolute path, not relative
addpath([BASEFOLDER 'gui']);
addpath([BASEFOLDER 'gui' filesep 'freezeColors']);
addpath([BASEFOLDER 'tsao_jiang']);
addpath([BASEFOLDER 'bydder']);
addpath([BASEFOLDER 'berglund']);
addpath([BASEFOLDER 'berglund' filesep 'Common']);
addpath([BASEFOLDER 'berglund' filesep 'FlexibleTwoPoint']);
addpath([BASEFOLDER 'berglund' filesep 'MultiSeedRegionGrowing']);
addpath([BASEFOLDER 'hernando' filesep 'graphcut']);
addpath([BASEFOLDER 'hernando' filesep 'common']);
addpath([BASEFOLDER 'hernando' filesep 'descent']);
addpath([BASEFOLDER 'hernando' filesep 'mixed_fitting']);
addpath([BASEFOLDER 'hernando' filesep 'matlab_bgl']);
addpath([BASEFOLDER 'eggers']);
addpath([BASEFOLDER 'lu']);
addpath([BASEFOLDER 'lu' filesep 'multiResSep']);
addpath([BASEFOLDER 'hu']);
addpath([BASEFOLDER 'hu' filesep 'Support Functions']);
addpath([BASEFOLDER 'sharma']);
addpath([BASEFOLDER 'sharma' filesep 'utils']);
addpath([BASEFOLDER 'pineda_bice']);


% Create the GUI
fw_gui;


