%%
% outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams)
% Description: Fat-water separation by hierarchical decomposition.
%              Version of method described, modified to handle multipeak 
%              and arbitrary TE
% Jiang Y, Tsao J. Fast and Robust Separation of Multiple Chemical Species 
% from Arbitrary Echo Times  with Complete Immunity to Phase Wrapping.
% In: Proceedings of the 20th Annual Meeting of ISMRM, Melbourne, Australia
% 2012
%
% Some properties:
%   - Image-space
%   - Any species
%   - Complex-fitting
%   - Multi-peak fat
%   - Single-R2*
%   - Independent water/fat phase
%   - Handles 3+ echoes at arbitrary echo times
%
% Input: structures imDataParams and algoParams
%   - imDataParams.images - images with dimensions (x,y,z,coils,TE)
%   - imDataParams.TE - TE in seconds
%   - imDataParams.FieldStrength - Field strength in Tesla
%   - algoParams.species - Species to be resolved.
%              Default: water and multi-peak fat
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
% Output: structure outParams
%   - outParams.water: water map with dimensions (x,y,z)
%   - outParams.fat: fat map with dimensions (x,y,z)
%   - outParams.fiterror: residual map with dimensions (x,y,z)
%   - outParams.phasemap: map of phase difference during delta TE, with
%                         dimensions (x,y,z). Equivalent to 
%                         exp(+i*2*pi*B0map*gyromagnetic*dTE + 2i/3*pi)
%   - outParams.r2starmap: R2* map with dimensions (x,y,z). (NOTE: R2* may
%                          be inaccurate due to the limited TE range)
%   - outParams.TE: TE in seconds
%   - outParams.FieldStrength - Field strength in Tesla
%   - outParams.WaterFatPpmDiff - ppm difference between water and fat peaks
%
% Authors: Jeffrey Tsao and Yun Jiang
% Date created: Oct 15, 2011
% Date last modified:

function outParams = fw_i2cm0c_3pluspoint_tsaojiang(imDataParams, algoParams)
    if nargin<1, help(mfilename); return; end
    if nargin<2, algoParams=[]; end

%%  % Default parameters
    %--------------------------------------------------------
    DefaultSpecies.water.name = 'water';
    DefaultSpecies.water.frequency = 0;
    DefaultSpecies.water.relAmps = 1;
    DefaultSpecies.fat.name = 'fat';
    DefaultSpecies.fat.frequency = [-3.80, -3.40, -2.60, -1.94, -0.39, 0.60];  % ppm
    DefaultSpecies.fat.relAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
 
    if ~isfield(algoParams,'species'), algoParams.species = []; end
    if ~isfield(algoParams,'Verbose'), algoParams.Verbose = []; end
    if ~isfield(algoParams,'AlwaysShowGUI'), algoParams.AlwaysShowGUI = []; end
    if ~isfield(algoParams,'Visualize'), algoParams.Visualize = []; end
    if ~isfield(algoParams,'Visualize_FatMapMultipler'), algoParams.Visualize_FatMapMultipler = []; end
    if ~isfield(algoParams,'MinFractSizeToDivide'), algoParams.MinFractSizeToDivide = []; end
    if ~isfield(algoParams,'MaxNumDiv'), algoParams.MaxNumDiv = []; end
    if ~isfield(algoParams,'AssumeSinglePeakAsWater'), algoParams.AssumeSinglePeakAsWater = []; end
    if ~isfield(algoParams,'SnrToAssumeSinglePeak'), algoParams.SnrToAssumeSinglePeak = []; end
    if ~isfield(algoParams,'CorrectAmpForT2star'), algoParams.CorrectAmpForT2star = []; end
    if ~isfield(algoParams,'MaxR2star'), algoParams.MaxR2star = []; end
    if isempty(algoParams.Verbose), algoParams.Verbose = 1; end
    if isempty(algoParams.AlwaysShowGUI), algoParams.AlwaysShowGUI = 0; end
    if isempty(algoParams.Visualize), algoParams.Visualize = 1; end
    if isempty(algoParams.Visualize_FatMapMultipler), algoParams.Visualize_FatMapMultipler = 1.5; end
    if isempty(algoParams.MinFractSizeToDivide), algoParams.MinFractSizeToDivide = 0.05; end
    if isempty(algoParams.MaxNumDiv), algoParams.MaxNumDiv = 6; end
    if isempty(algoParams.AssumeSinglePeakAsWater), algoParams.AssumeSinglePeakAsWater = 1; end
    if isempty(algoParams.SnrToAssumeSinglePeak), algoParams.SnrToAssumeSinglePeak = 2.5; end
    if isempty(algoParams.CorrectAmpForT2star), algoParams.CorrectAmpForT2star = 1; end
    if isempty(algoParams.MaxR2star), algoParams.MaxR2star = 250; end
    
    IDEALalgorithm = @hierarchical_IDEAL_engine_flexible;
    if algoParams.Verbose,
        fprintf('Algorithm: %s\n',IDEALalgorithm('name'));
    end

%%  % Parameter validation
    %--------------------------------------------------------
    MinNumIDEALechoes = 3;
    InputMissing_FieldStrength         = 0;
    InputMissing_PrecessionIsClockwise = 0;
    InputMissing_species               = 0;
       
    if ~isfield(imDataParams,'TE'),
        if isfield(imDataParams,'TEs'),
            imDataParams.TE = imDataParams.TEs;
            imDataParams = rmfield(imDataParams,'TEs');
        else
            error('TE field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'images'),
        if isfield(imDataParams,'image'),
            imDataParams.images = imDataParams.image;
            imDataParams = rmfield(imDataParams,'image');
        else
            error('images field is missing in imDataParams');
        end
    end
    if numel(imDataParams.TE)<MinNumIDEALechoes, error('Current implementation handles at least %d echoes only.',MinNumIDEALechoes); end
    if numel(imDataParams.TE)<numel(algoParams.species), error('There are fewer TE (%d) than the number of species (%d).',numel(imDataParams.TE),numel(algoParams.species)); end
    if size(imDataParams.images,5)~=numel(imDataParams.TE), error('The data have a number of echoes (%d) different from the number of TE (%d).', size(imDataParams.images,5),numel(imDataParams.TE)); end
    if ndims(imDataParams.images)>5, error('Current implementation handles image data with 5 dimensions only.'); end
    
%%  % Check manditory parameters
    %--------------------------------------------------------
    if  ~isfield(imDataParams,'FieldStrength') |  isempty(imDataParams.FieldStrength) | ~isnumeric(imDataParams.FieldStrength) | ~isreal(imDataParams.FieldStrength),
      if isfield(imDataParams,'fieldStrength') & ~isempty(imDataParams.fieldStrength) &  isnumeric(imDataParams.fieldStrength) &  isreal(imDataParams.fieldStrength),
        imDataParams.FieldStrength = imDataParams.fieldStrength;
        imDataParams = rmfield(imDataParams,'fieldStrength');
      else
        InputMissing_FieldStrength = 1;
        imDataParams.FieldStrength = [];
      end
    end  
    if ~isfield(imDataParams,'PrecessionIsClockwise') |  isempty(imDataParams.PrecessionIsClockwise) | ~isnumeric(imDataParams.PrecessionIsClockwise) | ~isreal(imDataParams.PrecessionIsClockwise),
      InputMissing_PrecessionIsClockwise =  1;
      imDataParams.PrecessionIsClockwise = [];
    end
    if ~isfield(algoParams,'species') | isempty(algoParams.species),
      InputMissing_species = 1;
      algoParams.species   = [];
    end

%%  % Get missing input
    %--------------------------------------------------------
    if algoParams.AlwaysShowGUI | InputMissing_species | InputMissing_PrecessionIsClockwise | InputMissing_FieldStrength,
      tmpParams = fw_inputparams(algoParams,imDataParams); drawnow;
      if isempty(tmpParams),  % Cancelled
        clear tmpParams;
        outParams = [];
        return;
      end
      if isfield(tmpParams,'species'              ), algoParams.species               = tmpParams.species; end
      if isfield(tmpParams,'PrecessionIsClockwise'), imDataParams.PrecessionIsClockwise = tmpParams.PrecessionIsClockwise; end
      if isfield(tmpParams,'FieldStrength'        ), imDataParams.FieldStrength         = tmpParams.FieldStrength; end
      clear tmpParams;
    end

%%  % Check if manditory parameters are now available and valid.
    % Otherwise, fill in defaults or ask for values
    %--------------------------------------------------------
%%  % Check algoParams.species
    %--------------------------------------------------------
    if isempty(algoParams.species),
      algoParams.species(1).name      = DefaultSpecies.water.name;
      algoParams.species(1).frequency = DefaultSpecies.water.frequency;
      algoParams.species(1).relAmps   = DefaultSpecies.water.relAmps;
      algoParams.species(2).name      = DefaultSpecies.fat.name;
      algoParams.species(2).frequency = DefaultSpecies.fat.frequency;
      algoParams.species(2).relAmps   = DefaultSpecies.fat.relAmps;
     else
      for n=1:numel(algoParams.species)
        if ~isfield(algoParams.species(n),'frequency') | isempty(algoParams.species(n).frequency), 
          if ~isfield(algoParams.species(n),'name')
            error('algoParams.species(%d).frequency is not defined.',n);
          end
          if strcmpi(algoParams.species(n).name,'water'),
            algoParams.species(n).frequency = DefaultSpecies.water.frequency;
            algoParams.species(n).relAmps   = DefaultSpecies.water.relAmps;
          elseif strcmpi(algoParams.species(n).name,'fat'),
            algoParams.species(n).frequency = DefaultSpecies.fat.frequency;
            algoParams.species(n).relAmps   = DefaultSpecies.fat.relAmps;
          end
        end
        if ~isfield(algoParams.species(n),'relAmps')
          algoParams.species(n).relAmps = ones(algoParams.species(n).frequency,1)/algoParams.species(n).frequency;
          fprintf('ASSUME that algoParams.species(%d).relAmps = [%.3f',n,algoParams.species(n).frequency(1));
          fprintf(' %.3f',algoParams.species(n).frequency(2:end));
          fprintf(']\n');
        end
      end; clear n;
    end
    if isempty(imDataParams.PrecessionIsClockwise),
      imDataParams.PrecessionIsClockwise = +1; %Default is +1
      fprintf('PrecessionIsClockwise field is missing in imDataParams\n');
      if imDataParams.PrecessionIsClockwise>0,
        fprintf('ASSUME that precession is clockwise.\n');
      else
        fprintf('ASSUME that precession is counter-clockwise.\n');
      end
    end
    if isempty(imDataParams.FieldStrength),
      while 1
        answer=inputdlg({'The field ''FieldStrength'' is missing in imDataParams.What is the field strength in Tesla?'},...
               'Warning: Parameter missing',1,{'1.5'});
        if isempty(answer), % Pressed cancel
          error('FieldStrength field is missing in imDataParams');
        end
        answer = sscanf(answer{1},'%f');
        if ~isempty(answer), imDataParams.FieldStrength = answer; break; end
      end; clear answer;
    end
       
%%  % Coil combination
    %--------------------------------------------------------
    tic;
    if algoParams.Verbose,
      fprintf('Phase correction');
      if size(imDataParams.images,4)>1, fprintf(' and %d-coil combination',size(imDataParams.images,4)); end
    end
    img = CombineCoilImg(imDataParams.images); clear relphase;
    imDataParams = rmfield(imDataParams,'images');   % Reduce memory
    if algoParams.Verbose, fprintf(' (%.2fs)\n',toc); end
    
%%  % TE step
    %--------------------------------------------------------
    tmpTE = sort(imDataParams.TE);
    dTE = median(tmpTE(2:end)-tmpTE(1:end-1)); 
    clear tmpTE;
    dTE_multiples = imDataParams.TE./dTE;

%%  % Conjugation
    %--------------------------------------------------------
    if imDataParams.PrecessionIsClockwise<=0,
      img = conj(img);
%      imDataParams.PrecessionIsClockwise = +1;  % Fixed
    end

%%  % Sort by TE
    %--------------------------------------------------------
    [tmpTE,tmpidx] = sort(imDataParams.TE);
    if ~isequal(tmpidx(:).',(1:numel(imDataParams.TE))), % Sort data according to TE
      imDataParams.TE = tmpTE;
      img = img(:,:,:,:,tmpidx);
    end
    clear tmpTE tmpidx;

%%
    % Figure out which slice to start
    % (Start from slice with most above-threshold pixels)
    tic;
    %level = setthreshold(img);           % Estimate intensity threshold
    %prof = zeros(size(img,3),1);
    %for z=1:size(img,3), % Find number of pixels above threshold
    %    prof(z) = sum(sum(sum(abs(img(:,:,z,:))>level*1.5)));
    %end; clear z;
    prof = mean(mean(mean(abs(img),1),2),5);   % Fixed bug 3 -> 5
    [maxval,maxz] = max(prof); clear maxval prof;
    if algoParams.Verbose, fprintf('Calculated central slice (%d out of %d) to start (%.2fs)\n',maxz,size(img,3),toc); end
    
%%  % Visualize
    if algoParams.Visualize,
        fig = figure;
        cmap = [gray(round(256/algoParams.Visualize_FatMapMultipler)); ones(256-round(256/algoParams.Visualize_FatMapMultipler),3)];
    end
    
%%
   % Create inversion matrix
   % NOTE:   Mat * species = data
   %         error = (Mat * pinv(Mat) - I) data
   %  error' error = data' (Mat * pinv(Mat) - I)' (Mat * pinv(Mat) - I) data
   %               = data' (I - Mat * pinv(Mat)) data
   %               = data' MatForMin data
   GyromagneticRatio = 42.576;         % MHz/T
   LarmorFreq = imDataParams.FieldStrength*GyromagneticRatio;
   MaxCond = 1000;
   Mat = zeros(size(img,5), numel(algoParams.species));   % TE x species   
   for n=1:numel(algoParams.species)
     freq = LarmorFreq * -algoParams.species(n).frequency;  % freq in Hz, LarmorFreq in MHz, (-algoParams.species(n).frequency) in ppm
     Mat(:,n) = exp(-1i*2*pi*(imDataParams.TE(:) * freq(:).')) * algoParams.species(n).relAmps(:);
   end; clear n freq LarmorFreq GyromagneticRatio;
   [tmpU, tmpS, tmpV] = svd(Mat,0);
   tmpS = diag(tmpS);
   tmpidx1 = find(tmpS> tmpS/MaxCond);
   tmpidx2 = find(tmpS<=tmpS/MaxCond);
   CondNum = tmpS(1)/min(tmpS(tmpidx1));
   tmpS(tmpidx1)=1./tmpS(tmpidx1);
   tmpS(tmpidx2)=0;
   InvMat = tmpV*diag(tmpS)*tmpU';
   MatForMin = eye(size(tmpU,1)) - tmpU*tmpU';
   if algoParams.Verbose,
     if isempty(tmpS(tmpidx2)),
       fprintf('Inversion matrix: condition number = %f\n',CondNum);
     else
       fprintf('WARNING: ill-conditioned inversion matrix (condition number limited to %f)\n',CondNum);
     end
   end
   clear tmpU tmpS tmpV tmpidx1 tmpidx2 MaxCond;
   
   % Extract lower triangle of MatForMin
   CoefForMin = zeros((size(MatForMin,1)+1)*size(MatForMin,1)/2,1);
   n = 0;
   for row = 1:size(MatForMin,1)
     CoefForMin(n+(1:row)) = MatForMin(row,(1:row));
     n=n+row;
   end; clear row n;
   
%%     
    for n=1:numel(algoParams.species)
      if isfield(algoParams.species(n),'name'     ), outParams.species(n).name      = algoParams.species(n).name;      end
      if isfield(algoParams.species(n),'frequency'), outParams.species(n).frequency = algoParams.species(n).frequency; end
      if isfield(algoParams.species(n),'relAmps'  ), outParams.species(n).relAmps   = algoParams.species(n).relAmps;   end
      algoParams.species(n).amps = zeros(size(img,1),size(img,2),size(img,3));
    end; clear n;
%    outParams.fiterror  = zeros(size(img,1),size(img,2),size(img,3));
    outParams.phasemap  = zeros(size(img,1),size(img,2),size(img,3));
    outParams.r2starmap = zeros(size(img,1),size(img,2),size(img,3));
    if isfield(imDataParams,'FieldStrength'),        outParams.FieldStrength        = imDataParams.FieldStrength;     end
    if isfield(imDataParams,'TE'           ),        outParams.TE                   = imDataParams.TE;                end

    if algoParams.AssumeSinglePeakAsWater
      % Figure out which species is which
      SpeciesIdx_water=[];
      SpeciesIdx_fat  =[];
      for n=1:numel(algoParams.species)
        if isequal(lower(algoParams.species(n).name),'water'), SpeciesIdx_water = n;
        elseif isequal(lower(algoParams.species(n).name),'fat'), SpeciesIdx_fat = n;
        end
      end
    end
    
    StartingAngleRangeGlobal = 180;
    StartingAngleRangeNarrow = 90;
    
%%  % Quick preprocess
    %--------------------------------------------------------
    QuickParams.MinFractSizeToDivide = 0.1;
    QuickParams.MaxNumDiv = 4;
    origPhaseMap=1;
    lastPhaseMap=1;
    AngleOffset = 0.0;
    zprocessed = 0;
    tic;
    fprintf('Preprocessing with low-res correction');
    StartingAngleRange = StartingAngleRangeGlobal;
    for z=[maxz:-1:1,maxz+1:size(img,3)]
        tmpimg = single(img(:,:,z,:)); % x,y,z,TE
        
%%      % Use last phase map as starting point for correction
        if z~=maxz,
            if z==maxz-1 || z==maxz+1, lastPhaseMap = origPhaseMap; end
            for TE=1:size(tmpimg,4)
              tmpimg(:,:,1,TE) = tmpimg(:,:,1,TE).*(lastPhaseMap.^dTE_multiples(TE));
            end; clear TE;
        end

%%      % Averaging neighboring pixels to improve SNR tmpimg = (x,y,z,coil,TE)
        for TE=1:size(tmpimg,4)
            tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg([2:end,1  ],:,1,TE)+tmpimg([end,1:(end-1)        ],:,1,TE)) ...
                + 0.0833*(tmpimg([3:end,1:2],:,1,TE)+tmpimg([(end-1):end,1:(end-2)],:,1,TE));
            tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg(:,[2:end,1  ],1,TE)+tmpimg(:,[end,1:(end-1)        ],1,TE)) ...
                + 0.0833*(tmpimg(:,[3:end,1:2],1,TE)+tmpimg(:,[(end-1):end,1:(end-2)],1,TE));
        end; clear TE;

%%      % Calculate phase
        [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, CoefForMin, dTE_multiples, QuickParams.MinFractSizeToDivide, StartingAngleRange, AngleOffset, QuickParams.MaxNumDiv);
        StartingAngleRange = StartingAngleRangeNarrow;
        if z==maxz,        % For first slice, apply correction and try again to accommodate very large shifts.
            lastPhaseMap = phasemap;
            for TE=1:size(tmpimg,4),
                tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(lastPhaseMap.^dTE_multiples(TE));
            end; clear TE tmppower;
            [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, CoefForMin, dTE_multiples, QuickParams.MinFractSizeToDivide,StartingAngleRange,AngleOffset, QuickParams.MaxNumDiv);
        end
        phasemap = phasemap.*lastPhaseMap; % cumulative
        lastPhaseMap = phasemap;

        if z==maxz, origPhaseMap = lastPhaseMap; end
    
%%      % Apply IDEAL correction
        tmpimg = double(img(:,:,z,:)); % x,y,z,TE
        for TE=1:size(tmpimg,4)
            img(:,:,z,TE) = img(:,:,z,TE).*(phasemap.^dTE_multiples(TE));
        end; clear TE;
        outParams.phasemap(:,:,z) = phasemap;
        zprocessed = zprocessed + 1;         
    end; clear z zprocessed phasemap t2correctmap;
    fprintf(' (%.2fs)\n',toc);
    clear QuickParams origPhaseMap lastPhaseMap StartingAngleRange AngleOffset;
    %-------------------------------------------------------------------       

    if algoParams.AssumeSinglePeakAsWater
      % Figure out which species is which
      SpeciesIdx_water=[];
      SpeciesIdx_fat  =[];
      for n=1:numel(algoParams.species)
        if isequal(lower(algoParams.species(n).name),'water'), SpeciesIdx_water = n;
        elseif isequal(lower(algoParams.species(n).name),'fat'), SpeciesIdx_fat = n;
        end
      end
    end
    
%%  % Slice by slice
    %--------------------------------------------------------
    origPhaseMap=1;
    lastPhaseMap=1;
    AngleOffset = 0.0;
    zprocessed = 0;
    tic;
    StartingAngleRange = StartingAngleRangeNarrow;
    for z=[maxz:-1:1,maxz+1:size(img,3)]
        tmpimg = single(img(:,:,z,:)); % x,y,z,TE

%%      % Use last phase map as starting point for correction
        if z~=maxz,
            if z==maxz-1 || z==maxz+1, lastPhaseMap = origPhaseMap; end
            for TE=1:size(tmpimg,4)
              tmpimg(:,:,1,TE) = tmpimg(:,:,1,TE).*(lastPhaseMap.^dTE_multiples(TE));
            end; clear TE;
        end

%%      % Averaging neighboring pixels to improve SNR tmpimg = (x,y,z,coil,TE)
        for TE=1:size(tmpimg,4)
            tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg([2:end,1  ],:,1,TE)+tmpimg([end,1:(end-1)        ],:,1,TE)) ...
                + 0.0833*(tmpimg([3:end,1:2],:,1,TE)+tmpimg([(end-1):end,1:(end-2)],:,1,TE));
            tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg(:,[2:end,1  ],1,TE)+tmpimg(:,[end,1:(end-1)        ],1,TE)) ...
                + 0.0833*(tmpimg(:,[3:end,1:2],1,TE)+tmpimg(:,[(end-1):end,1:(end-2)],1,TE));
        end; clear TE;

%%      % Calculate phase
        [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, CoefForMin, dTE_multiples, algoParams.MinFractSizeToDivide, StartingAngleRange, AngleOffset, algoParams.MaxNumDiv);
        StartingAngleRange = StartingAngleRangeNarrow;
        if z==maxz,        % For first slice, apply correction and try again to accommodate very large shifts.
          lastPhaseMap = phasemap;
          for TE=1:size(tmpimg,4),
            tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(lastPhaseMap.^dTE_multiples(TE));
          end; clear TE tmppower;
          [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, CoefForMin, dTE_multiples, algoParams.MinFractSizeToDivide,StartingAngleRange,AngleOffset, algoParams.MaxNumDiv);
        end
        
        phasemap = phasemap.*lastPhaseMap; % cumulative
        lastPhaseMap = phasemap;

%%      % Blur it
        if 1,
            lastPhaseMap = lastPhaseMap.*sqrt(sum(abs(tmpimg).^2,4));  % magnitude weighting
            lastPhaseMap = convn(lastPhaseMap,ones(max(ceil([size(lastPhaseMap,1),size(lastPhaseMap,2),size(lastPhaseMap,3)]*algoParams.MinFractSizeToDivide/2),1)),'same');
            lastPhaseMap = exp(1i*angle(lastPhaseMap));
        end
        if z==maxz, origPhaseMap = lastPhaseMap; end

        %%      % R2* map
        r2starmap = max(t2correctmap,1); % r2starmap(r2starmap<=0)=1;
        r2starmap = log(r2starmap)/dTE;
        if any(r2starmap(:)>algoParams.MaxR2star),
          r2starmap = min(r2starmap,algoParams.MaxR2star);
          t2correctmap = exp(r2starmap*dTE);
        end

%%      % Apply IDEAL correction
        tmpimg = double(img(:,:,z,:)); % x,y,z,TE
        if algoParams.CorrectAmpForT2star,   % Correct
          dMap = t2correctmap.*phasemap;
          for TE=1:size(tmpimg,4)
            tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(dMap.^dTE_multiples(TE));
          end; clear TE;
          clear dMap;
        else
          for TE=1:size(tmpimg,4)
            tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(phasemap.^dTE_multiples(TE));
          end; clear TE;
        end
        clear t2correctmap;
        % tmpimg=transformKspaceToImage(tmpimg,4) / size(tmpimg,4);

        tmpimg = reshape(tmpimg,[size(tmpimg,1)*size(tmpimg,2)*size(tmpimg,3),size(tmpimg,4)]);
        fiterror = sqrt(real(reshape(sum((tmpimg*(MatForMin.')) .* conj(tmpimg), 2),[size(img,1),size(img,2),1])));
        tmpimg = reshape(tmpimg*(InvMat.'),[size(img,1),size(img,2),1,size(InvMat,1)]);
        
%%      % Use heuristics: if there is only one peak, check if it is 
        % in fat position. If so, shift it.
        %   tmpimg(:,:,1,1) should be water
        %   tmpimg(:,:,1,2) should be fat
        %   tmpimg(:,:,1,3) should be blank
        % i.e. if 1. |fat| > SNR x max(|water|,|blank|)
        %         2. max(|water|,|blank|) < SNR x min(|water|,|blank|)
        if algoParams.AssumeSinglePeakAsWater & ~isempty(SpeciesIdx_water) & ~isempty(SpeciesIdx_fat)
            % Get rid of phase and then sum
            BinValue_fat   = tmpimg(:,:,:,SpeciesIdx_fat);
            BinValue_water = tmpimg(:,:,:,SpeciesIdx_water);
            BinValue_blank = fiterror;
            
            if size(BinValue_fat,1)>1
              BinValue_fat   = BinValue_fat(1:end-1,:,:,:).*exp(-1i*angle(BinValue_fat(2:end,:,:,:)));
              BinValue_water = BinValue_water(1:end-1,:,:,:).*exp(-1i*angle(BinValue_water(2:end,:,:,:)));
              BinValue_blank = BinValue_blank(1:end-1,:,:,:);
            end
            if size(BinValue_fat,2)>1
              BinValue_fat   = BinValue_fat(:,1:end-1,:,:).*exp(-1i*angle(BinValue_fat(:,2:end,:,:)));
              BinValue_water = BinValue_water(:,1:end-1,:,:).*exp(-1i*angle(BinValue_water(:,2:end,:,:)));
              BinValue_blank = BinValue_blank(:,1:end-1,:,:);
            end
            if size(BinValue_fat,3)>1
              BinValue_fat   = BinValue_fat(:,:,1:end-1,:).*exp(-1i*angle(BinValue_fat(:,:,2:end,:)));
              BinValue_water = BinValue_water(:,:,1:end-1,:).*exp(-1i*angle(BinValue_water(:,:,2:end,:)));
              BinValue_blank = BinValue_blank(:,:,1:end-1,:);
            end
            BinValue_fat   = abs(mean(mean(mean(BinValue_fat,1),2),3));
            BinValue_water = abs(mean(mean(mean(BinValue_water,1),2),3));
            BinValue_blank = abs(mean(mean(mean(BinValue_blank,1),2),3));

            maxwaterblank = max(BinValue_water,BinValue_blank);
            minwaterblank = min(BinValue_water,BinValue_blank);
            if BinValue_fat > maxwaterblank * algoParams.SnrToAssumeSinglePeak ...
            && maxwaterblank < minwaterblank * algoParams.SnrToAssumeSinglePeak ,
                if algoParams.Verbose, fprintf('S'); end
                tmperror = fiterror;
                fiterror = abs(tmp(:,:,:,SpeciesIdx_water));
                tmpimg(:,:,1,SpeciesIdx_water) = tmpimg(:,:,1,SpeciesIdx_fat);
                tmpimg(:,:,1,SpeciesIdx_fat) = tmperror; clear tmperror;
            else
              if algoParams.Verbose, fprintf('.'); end
            end
            clear BinValue_fat BinValue_water BinValue_blank maxwaterblank minwaterblank;
        else
          if algoParams.Verbose, fprintf('.'); end          
        end


%%      % Assign
        for n=1:numel(algoParams.species)
          outParams.species(n).amps(:,:,z) = tmpimg(:,:,1,n);
        end; clear n;
        outParams.fiterror(:,:,z) = fiterror;
        outParams.phasemap(:,:,z) = outParams.phasemap(:,:,z).*phasemap;
        outParams.r2starmap(:,:,z) = r2starmap;
        clear tmpimg wateridx fatidx blankidx;

%%      % Visualize
        if algoParams.Visualize,
          try,
            figure(fig); colormap(cmap);
            caxisval=[];
            for n=1:numel(algoParams.species)
              h = subplot(2,numel(algoParams.species),n,'Replace',fig);
              imagesc(abs(outParams.species(n).amps(:,:,z)).'   ,'Parent',h);
              if isfield(outParams.species(n),'name'),
                title(h,sprintf('%s (slice %d)',outParams.species(n).name,z)); 
              else
                title(h,sprintf('%d (slice %d)',n,z)); 
              end
              axis(h,'image'); xlabel(h,''); ylabel(h,'');
              if isempty(caxisval),
                caxisval=caxis; 
              else
                tmpcaxisval = caxis; 
                caxisval(1)=min(caxisval(1),tmpcaxisval(1)); 
                caxisval(2)=max(caxisval(2),tmpcaxisval(2)); 
                clear tmpcaxisval;
              end
            end; clear n h;
            h=subplot(2,3,4,'Replace',fig); imagesc(real(fiterror(:,:)).','Parent',h); title(h,sprintf('Fit error (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            tmpcaxisval = caxis;  caxisval(1)=min(caxisval(1),tmpcaxisval(1));  caxisval(2)=max(caxisval(2),tmpcaxisval(2));  clear tmpcaxisval;
            for n=1:numel(algoParams.species)
              h = subplot(2,numel(algoParams.species),n,'Parent',fig);
              caxis(h,caxisval);
            end; clear n h;
            h = subplot(2,3,4,'Parent',fig); caxis(h,caxisval);
            h=subplot(2,3,5,'Replace',fig); imagesc(permute(ComplexToRgb(phasemap(:,:),[],jet),[2,1,3]),'Parent',h); title(h,sprintf('\\DeltaPhase (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            h=subplot(2,3,6,'Replace',fig); imagesc(r2starmap(:,:).','Parent',h); title(h,sprintf('R_2* (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            drawnow; clear h caxisval;
          end
        end
        zprocessed = zprocessed + 1;         
        if algoParams.Verbose, 
          if mod(zprocessed,10)==0 || zprocessed==size(img,3),
            fprintf(' (%.1f%% done, %.2fs)\n',zprocessed/size(img,3)*100, toc);
          end
        end
    end; clear z zprocessed phasemap t2correctmap fiterror r2starmap ;
    clear Mat InvMat MatForMin;
    
%%  % Close figure
    if algoParams.Visualize,
        try delete(fig); catch, end
        clear fig cmap;
    end
    %--------------------------------------------------------
    clear img;
end


% Averaging neighboring pixels to improve SNR tmpimg = (x,y,z,coil,TE)
function tmpimg = avgneighbor(tmpimg)
  for TE=1:size(tmpimg,4)
     tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg([2:end,1  ],:,1,TE)+tmpimg([end,1:(end-1)        ],:,1,TE)) ...
                + 0.0833*(tmpimg([3:end,1:2],:,1,TE)+tmpimg([(end-1):end,1:(end-2)],:,1,TE));
     tmpimg(:,:,1,TE) = 0.3333* tmpimg(:,:,1,TE) ...
                + 0.2500*(tmpimg(:,[2:end,1  ],1,TE)+tmpimg(:,[end,1:(end-1)        ],1,TE)) ...
                + 0.0833*(tmpimg(:,[3:end,1:2],1,TE)+tmpimg(:,[(end-1):end,1:(end-2)],1,TE));
  end; clear TE;
end