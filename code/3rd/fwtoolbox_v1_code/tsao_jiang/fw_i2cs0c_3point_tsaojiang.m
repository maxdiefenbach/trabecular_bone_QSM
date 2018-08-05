%%
% outParams = fw_i2cs0c_3point_tsaojiang(imDataParams, algoParams)
% Description: Fat-water separation by hierarchical decomposition and
%              direct estimation of phase offset to locate signal null.
% Tsao J, Jiang Y. Hierarchical IDEAL: robust water–fat separation at high
% field by multiresolution field map estimation. In: Proceedings of the 18th
% Annual Meeting of ISMRM, Toronto, ON, Canada, 2008. p 653
%
% Some properties:
%   - Image-space
%   - 2 species (water-fat)
%   - Complex-fitting
%   - Single-peak fat
%   - Single-R2*
%   - Independent water/fat phase
%   - Requires 3 echoes at optimized echo times
%
% Input: structures imDataParams and algoParams
%   - imDataParams.images - images with dimensions (x,y,z,coils,TE)
%   - imDataParams.TE - TE in seconds
%   - imDataParams.FieldStrength - Field strength in Tesla
%   - algoParams.Verbose - (optional) 1 = show info, 0 = no info
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
%   - algoParams.CheckTE - (optional) 1 = Check TE to ensure that they are
%       are the optimized echo times. They are needed for the algorithm.
%       (default 1)
%   - algoParams.CheckTE_WaterFatPpmDiff - (optional) ppm difference
%       between water and fat peaks used in checking TE. (default 3.4)
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
% Date created: September 6, 2011
% Date last modified: September 7, 2011

function outParams = fw_i2cs0c_3point_tsaojiang(imDataParams, algoParams)
    if nargin<1, help(mfilename); return; end
    if nargin<2, algoParams=[]; end

    tic

%%  % Default parameters
    %--------------------------------------------------------
    if ~isfield(algoParams,'Verbose'), algoParams.Verbose = []; end
    if ~isfield(algoParams,'Visualize'), algoParams.Visualize = []; end
    if ~isfield(algoParams,'Visualize_FatMapMultipler'), algoParams.Visualize_FatMapMultipler = []; end
    if ~isfield(algoParams,'MinFractSizeToDivide'), algoParams.MinFractSizeToDivide = []; end
    if ~isfield(algoParams,'MaxNumDiv'), algoParams.MaxNumDiv = []; end
    if ~isfield(algoParams,'AssumeSinglePeakAsWater'), algoParams.AssumeSinglePeakAsWater = []; end
    if ~isfield(algoParams,'SnrToAssumeSinglePeak'), algoParams.SnrToAssumeSinglePeak = []; end
    if ~isfield(algoParams,'CheckTE'), algoParams.CheckTE = []; end
    if ~isfield(algoParams,'CheckTE_WaterFatPpmDiff'), algoParams.CheckTE_WaterFatPpmDiff = []; end
    if ~isfield(algoParams,'CorrectAmpForT2star'), algoParams.CorrectAmpForT2star = []; end
    if ~isfield(algoParams,'MaxR2star'), algoParams.MaxR2star = []; end
%    if ~isfield(algoParams,'PhaseIsInverted'), algoParams.PhaseIsInverted = []; end
    if isempty(algoParams.Verbose), algoParams.Verbose = 1; end
    if isempty(algoParams.Visualize), algoParams.Visualize = 1; end
    if isempty(algoParams.Visualize_FatMapMultipler), algoParams.Visualize_FatMapMultipler = 1.5; end
    if isempty(algoParams.MinFractSizeToDivide), algoParams.MinFractSizeToDivide = 0.05; end
    if isempty(algoParams.MaxNumDiv), algoParams.MaxNumDiv = 6; end
    if isempty(algoParams.AssumeSinglePeakAsWater), algoParams.AssumeSinglePeakAsWater = 0; end
    if isempty(algoParams.SnrToAssumeSinglePeak), algoParams.SnrToAssumeSinglePeak = 2.5; end
    if isempty(algoParams.CheckTE), algoParams.CheckTE = 1; end
    if isempty(algoParams.CheckTE_WaterFatPpmDiff), algoParams.CheckTE_WaterFatPpmDiff = 3.4; end
    algoParams.CheckTE_WaterFatPpmDiff = abs(algoParams.CheckTE_WaterFatPpmDiff);
    if isempty(algoParams.CorrectAmpForT2star), algoParams.CorrectAmpForT2star = 1; end
    if isempty(algoParams.MaxR2star), algoParams.MaxR2star = 250; end
%    if isempty(algoParams.PhaseIsInverted), algoParams.PhaseIsInverted = 0; end
    
    IDEALalgorithm = @hierarchical_IDEAL_engine;
    if algoParams.Verbose,
        fprintf('Algorithm: %s\n',IDEALalgorithm('name'));
    end

%%  % Parameter validation
    %--------------------------------------------------------
    NumIDEALechoes = 3;
    if ~isfield(imDataParams,'TE'),
        if isfield(imDataParams,'TEs'),
            imDataParams.TE = imDataParams.TEs;
            rmfield(imDataParams,'TEs');
        else
            error('TE field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'images'),
        if isfield(imDataParams,'image'),
            imDataParams.images = imDataParams.image;
            rmfield(imDataParams,'image');
        else
            error('images field is missing in imDataParams');
        end
    end
    if numel(imDataParams.TE)~=NumIDEALechoes, error('Current implementation handles %d echoes only.',NumIDEALechoes); end
    if size(imDataParams.images,5)~=NumIDEALechoes, error('Current implementation handles %d echoes only.', NumIDEALechoes); end
    if ndims(imDataParams.images)>5, error('Current implementation handles image data with 5 dimensions only.'); end
    if algoParams.CheckTE,
      if ~isfield(imDataParams,'FieldStrength'), 
        if isfield(imDataParams,'fieldStrength'),
          imDataParams.FieldStrength = imDataParams.fieldStrength;
          rmfield(imDataParams,'fieldStrength');
        else
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
      end
      status = VerifyEchoTimes(imDataParams.TE, imDataParams.FieldStrength, algoParams.CheckTE_WaterFatPpmDiff, algoParams.Verbose);
      if status>0,
        if status==1,  % Incorrect number of echoes
          error('Current implementation handles %d optimized echo times only,\naccording to Pineda et al. Magn Reson Med. 2005 Sep;54(3):625-35.\nNOTE: You can turn off echo-time checking with algoParams.CheckTE=0\n      The water-fat shift is set by algoParams.CheckTE_WaterFatPpmDiff (%.2f ppm)', NumIDEALechoes,algoParams.CheckTE_WaterFatPpmDiff);
        else % Non-optimal echo times
          str1 = 'Ignore';
          str2 = 'Stop';
          button = questdlg({'The current echo times are suboptimal according to Pineda et al.',...
                            'Magn Reson Med. 2005 Sep;54(3):625-35. This may be:',...
                            '  1. a real error,',...
                            '  2. incorrect setting for field strength',...
                    sprintf('      (imDataParams.FieldStrength = %.2fT), or',imDataParams.FieldStrength),...
                            '  3. incorrect setting for chemical shift between water and fat',...
                    sprintf('      (algoParams.CheckTE_WaterFatPpmDiff = %.2fppm)',algoParams.CheckTE_WaterFatPpmDiff),...
                            ' ',...
                            'You may ignore this error, which may produce incorrect images',...
                            'if the error is real, or you can stop.'},...
                            'Warning: Suboptimal TE',str1,str2,str2);
          if isequal(button,str2),
             error('Current implementation handles %d optimized echo times only,\naccording to Pineda et al. Magn Reson Med. 2005 Sep;54(3):625-35.\nNOTE: You can turn off echo-time checking with algoParams.CheckTE=0\n      The water-fat shift is set by algoParams.CheckTE_WaterFatPpmDiff (%.2f ppm)', NumIDEALechoes,algoParams.CheckTE_WaterFatPpmDiff);
          end
          clear str1 str2 button;
        end
      end
      if status==-1, % Swapped echo time
          imDataParams.images = imDataParams.images(:,:,:,:,[3,2,1]);
          fprintf('Swapped echo times\n');
      end
      clear status
    end
    if ~isfield(imDataParams,'PrecessionIsClockwise'),
      imDataParams.PrecessionIsClockwise = +1;  % -1;
      fprintf('PrecessionIsClockwise field is missing in imDataParams\n');
      if imDataParams.PrecessionIsClockwise>0,
        fprintf('ASSUME that precession is clockwise.\n');
      else
        fprintf('ASSUME that precession is counter-clockwise.\n');
      end
    end
    if imDataParams.PrecessionIsClockwise<=0,
      imDataParams.images = conj(imDataParams.images);
    end
    %--------------------------------------------------------

%%  % Coil combination
    %--------------------------------------------------------
    tic;
    if algoParams.Verbose,
      fprintf('Phase correction');
      if size(imDataParams.images,4)>1, fprintf(' and %d-coil combination',size(imDataParams.images,4)); end
    end
    img = CombineCoilImg(imDataParams.images); clear relphase;
    rmfield(imDataParams,'images');  % Reduce memory requirement
    if algoParams.Verbose, fprintf(' (%.2fs)\n',toc); end
%    if algoParams.PhaseIsInverted, img=conj(img); end
    
%%  %% Correct linear phase (i.e. shift in k-space)
    %%--------------------------------------------------------
    %img = CorrectLinearPhase(img,algoParams.Verbose);
    img = reshape(img,[size(img,1),size(img,2),size(img,3),size(img,5)]); % -> x,y,z,TE

%%
    % Figure out which slice to start
    % (Start from slice with most above-threshold pixels)
    tic;
    %level = setthreshold(img);           % Estimate intensity threshold
    %prof = zeros(size(img,3),1);
    %for z=1:size(img,3), % Find number of pixels above threshold
    %    prof(z) = sum(sum(sum(abs(img(:,:,z,:))>level*1.5)));
    %end; clear z;
    prof = mean(mean(mean(abs(img),1),2),4);   % Fixed bug 3 -> 4
    [maxval,maxz] = max(prof); clear maxval prof;
    if algoParams.Verbose, fprintf('Calculated central slice to start (%.2fs)\n',toc); end
    
%%  % Visualize
    if algoParams.Visualize,
        fig = figure;
        cmap = [gray(round(256/algoParams.Visualize_FatMapMultipler)); ones(256-round(256/algoParams.Visualize_FatMapMultipler),3)];
    end
    
%%
    outParams.water = zeros(size(img,1),size(img,2),size(img,3));
    outParams.fat   = zeros(size(img,1),size(img,2),size(img,3));
    outParams.fiterror = zeros(size(img,1),size(img,2),size(img,3));
    outParams.phasemap = zeros(size(img,1),size(img,2),size(img,3));
    outParams.r2starmap = zeros(size(img,1),size(img,2),size(img,3));
    outParams.TE = imDataParams.TE;
    if isfield(imDataParams,'FieldStrength'), 
      outParams.FieldStrength = imDataParams.FieldStrength;
    end
    outParams.WaterFatPpmDiff = algoParams.CheckTE_WaterFatPpmDiff;

%%  % Quick preprocess
    %--------------------------------------------------------
    QuickParams.MinFractSizeToDivide = 0.1;
    QuickParams.MaxNumDiv = 3;
    origPhaseMap=1;
    lastPhaseMap=1;
    StartingAngleRange = 180;
    AngleOffset = 0.0;
    zprocessed = 0;
    tic;
    fprintf('Preprocessing with low-res correction');
    for z=[maxz:-1:1,maxz+1:size(img,3)]
        %if algoParams.Verbose, fprintf('z=%d',z); end
        tmpimg = single(img(:,:,z,:)); % x,y,z,TE
        
%%      % Use last phase map as starting point for correction
        if z~=maxz,
            if z==maxz-1 || z==maxz+1, lastPhaseMap = origPhaseMap; end
            for TE=1:size(tmpimg,4)
                tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
                if tmppower~=0, tmpimg(:,:,1,TE) = tmpimg(:,:,1,TE).*(lastPhaseMap.^tmppower); end
            end; clear TE tmppower;
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
        [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, QuickParams.MinFractSizeToDivide, StartingAngleRange, AngleOffset, QuickParams.MaxNumDiv);
        if z==maxz,        % For first slice, apply correction and try again to accommodate very large shifts.
            lastPhaseMap = phasemap;
            for TE=1:size(tmpimg,4),
                tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
                if tmppower~=0, tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(lastPhaseMap.^tmppower); end
            end; clear TE tmppower;
            [phasemap,t2correctmap] = IDEALalgorithm(tmpimg,QuickParams.MinFractSizeToDivide,StartingAngleRange,AngleOffset, QuickParams.MaxNumDiv);
        end
        phasemap = phasemap.*lastPhaseMap; % cumulative
        lastPhaseMap = phasemap;

        if z==maxz, origPhaseMap = lastPhaseMap; end
    
%%      % Apply IDEAL correction
        for TE=1:size(tmpimg,4)
            tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
            if tmppower~=0, img(:,:,z,TE) = img(:,:,z,TE).*(phasemap.^tmppower); end
        end; clear TE tmppower;
        outParams.phasemap(:,:,z) = phasemap;
        zprocessed = zprocessed + 1;         
    end; clear z zprocessed phasemap t2correctmap;
    fprintf(' (%.2fs)\n',toc);
    clear QuickParams origPhaseMap lastPhaseMap StartingAngleRange AngleOffset;
    %-------------------------------------------------------------------   

%%  % Slice by slice
    %-------------------------------------------------------------------
    origPhaseMap=1;
    lastPhaseMap=1;
    StartingAngleRange = 180;
    AngleOffset = 0.0;
    zprocessed = 0;
    tic;
    for z=[maxz:-1:1,maxz+1:size(img,3)]
        %if algoParams.Verbose, fprintf('z=%d',z); end
        tmpimg = single(img(:,:,z,:)); % x,y,z,TE
        
%%      % Use last phase map as starting point for correction
        if z~=maxz,
            if z==maxz-1 || z==maxz+1, lastPhaseMap = origPhaseMap; end
            for TE=1:size(tmpimg,4)
                tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
                if tmppower~=0, tmpimg(:,:,1,TE) = tmpimg(:,:,1,TE).*(lastPhaseMap.^tmppower); end
            end; clear TE tmppower;
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
        [phasemap,t2correctmap] = IDEALalgorithm(tmpimg, algoParams.MinFractSizeToDivide, StartingAngleRange, AngleOffset, algoParams.MaxNumDiv);
        if z==maxz,        % For first slice, apply correction and try again to accommodate very large shifts.
            lastPhaseMap = phasemap;
            for TE=1:size(tmpimg,4),
                tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
                if tmppower~=0, tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(lastPhaseMap.^tmppower); end
            end; clear TE tmppower;
            [phasemap,t2correctmap] = IDEALalgorithm(tmpimg,algoParams.MinFractSizeToDivide,StartingAngleRange,AngleOffset, algoParams.MaxNumDiv);
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

%%      % Apply IDEAL correction
        tmpimg = double(img(:,:,z,:)); % x,y,z,TE
        for TE=1:size(tmpimg,4)
            tmppower = TE-(bitshift(size(tmpimg,4),-1)+1);
            if tmppower~=0, tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*(phasemap.^tmppower); end
        end; clear TE tmppower;
        r2starmap = t2correctmap;
        dTE = (imDataParams.TE(end)-imDataParams.TE(1))/(numel(imDataParams.TE)-1);
        r2starmap(r2starmap<=0)=1;
        r2starmap = log(r2starmap)/dTE;
        if any(r2starmap(:)>algoParams.MaxR2star),
          r2starmap = min(r2starmap,algoParams.MaxR2star);
          t2correctmap = exp(r2starmap*dTE);
        end
        if algoParams.CorrectAmpForT2star,   % Correct
          for TE=1:size(tmpimg,4)
            tmpimg(:,:,:,TE) = tmpimg(:,:,:,TE).*t2correctmap.^((TE-1)+imDataParams.TE(1)/dTE);
          end; clear TE;
        end
        clear t2correctmap dTE;
        tmpimg=transformKspaceToImage(tmpimg,4) / size(tmpimg,4);

%%      % Use heuristics: if there is only one peak, check if it is 
        % in fat position. If so, shift it.
        %   tmpimg(:,:,1,1) should be blank
        %   tmpimg(:,:,1,2) should be water
        %   tmpimg(:,:,1,3) should be fat
        % i.e. if 1. |fat| > SNR x max(|water|,|blank|)
        %         2. max(|water|,|blank|) < SNR x min(|water|,|blank|)
        blankidx = 1;
        wateridx = 2;
        fatidx = 3;
        if algoParams.AssumeSinglePeakAsWater
            % Get rid of phase and then sum
            if size(tmpimg,1)>1
                BinValues = tmpimg(1:end-1,:,:,:).*exp(-1i*angle(tmpimg(2:end,:,:,:)));
            else
                BinValues = tmpimg;
            end
            if size(BinValues,2)>1
                BinValues = BinValues(:,1:end-1,:,:).*exp(-1i*angle(BinValues(:,2:end,:,:)));
            end
            if size(BinValues,3)>1
                BinValues = BinValues(:,:,1:end-1,:).*exp(-1i*angle(BinValues(:,:,2:end,:)));
            end
            BinValues = permute(abs(mean(mean(mean(BinValues,1),2),3)),[4,1,2,3]);

            maxwaterblank = max(BinValues(wateridx),BinValues(blankidx));
            minwaterblank = min(BinValues(wateridx),BinValues(blankidx));
            if BinValues(fatidx) > maxwaterblank * algoParams.SnrToAssumeSinglePeak ...
            && maxwaterblank < minwaterblank * algoParams.SnrToAssumeSinglePeak ,
                if algoParams.Verbose, fprintf('S'); end
                tmpimg(:,:,1,[fatidx,wateridx,blankidx]) = tmpimg(:,:,1,[blankidx,fatidx,wateridx]);      % shift by one plane
            else
              if algoParams.Verbose, fprintf('.'); end
            end
            clear Bin_maxidx maxwaterblank minwaterblank;
        else
          if algoParams.Verbose, fprintf('.'); end          
        end

%%      % Assign
        outParams.water(:,:,z) = tmpimg(:,:,1,wateridx);
        outParams.fat(:,:,z) = tmpimg(:,:,1,fatidx);
        outParams.fiterror(:,:,z) = tmpimg(:,:,1,blankidx);
        outParams.phasemap(:,:,z) = outParams.phasemap(:,:,z).*phasemap;
        outParams.r2starmap(:,:,z) = r2starmap;
        clear tmpimg wateridx fatidx blankidx;

%%      % Visualize
        if algoParams.Visualize,
            figure(fig); colormap(cmap);
            h = subplot(2,3,1,'Replace',fig); imagesc(abs(outParams.water(:,:,z)).'   ,'Parent',h); title(h,sprintf('Water (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            caxisval=caxis; 
            h = subplot(2,3,2,'Replace',fig); imagesc(abs(outParams.fat(:,:,z)).'     ,'Parent',h); title(h,sprintf('Fat (slice %d)',z));   axis(h,'image'); xlabel(h,''); ylabel(h,'');
            if isempty(caxisval), caxisval=caxis; else tmpcaxisval = caxis; caxisval(1)=min(caxisval(1),tmpcaxisval(1)); caxisval(2)=max(caxisval(2),tmpcaxisval(2)); end; clear tmpcaxisval
            h = subplot(2,3,3,'Replace',fig); imagesc(abs(outParams.fiterror(:,:,z)).','Parent',h); title(h,sprintf('Blank (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            if isempty(caxisval), caxisval=caxis; else tmpcaxisval = caxis; caxisval(1)=min(caxisval(1),tmpcaxisval(1)); caxisval(2)=max(caxisval(2),tmpcaxisval(2)); end; clear tmpcaxisval
            
            h=subplot(2,3,1,'Parent',fig); caxis(h,caxisval); 
            h=subplot(2,3,2,'Parent',fig); caxis(h,caxisval); 
            h=subplot(2,3,3,'Parent',fig); caxis(h,caxisval); 
            h=subplot(2,2,3,'Replace',fig); imagesc(permute(ComplexToRgb(phasemap(:,:),[],jet),[2,1,3]),'Parent',h); title(h,sprintf('\\DeltaPhase (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            h=subplot(2,2,4,'Replace',fig); imagesc(r2starmap(:,:).','Parent',h); title(h,sprintf('R_2* (slice %d)',z)); axis(h,'image'); xlabel(h,''); ylabel(h,'');
            drawnow; clear h caxisval;
        end
        zprocessed = zprocessed + 1;         
        if algoParams.Verbose, 
          if mod(zprocessed,10)==0 || zprocessed==size(img,3),
            fprintf(' (%.1f%% done, %.2fs)\n',zprocessed/size(img,3)*100, toc);
          end
        end
    end; clear z zprocessed phasemap t2correctmap;

%%  % Close figure
    if algoParams.Visualize,
        try delete(fig); catch, end
        clear fig cmap;
    end
    %--------------------------------------------------------
    clear img;
end

