% outParams = ConvertSingleToMultipeakRecon(imDataParams,algoParams)
%
% 2011, Jeffrey Tsao

%
function outParams = ConvertSingleToMultipeakRecon(imDataParams,algoParams),
    if nargin<1, help(mfilename); return; end
    if nargin<2, algoParams=[]; end

    % Default parameters
    %--------------------------------------------------------
    GyromagneticRatio = 42.58;
    MaxPpmRange = 0.25;
    NumChemSpecies = 2;
    NumIDEALechoes = 3;
    DefaultWaterPpm  = 0;
    DefaultWaterRelAmps = 1;
    DefaultFatPpm  = [-3.80, -3.40, -2.60, -1.94, -0.39, +0.60];;
    DefaultFatRelAmps = [0.087 0.693 0.128 0.004 0.039 0.048];
    
    % Algorithm parameters
    %--------------------------------------------------------------
    if ~isfield(algoParams,'Verbose'), algoParams.Verbose = []; end
    if ~isfield(algoParams,'Visualize'), algoParams.Visualize = []; end
    if ~isfield(algoParams,'species'), algoParams.species = []; end
    if isempty(algoParams.Verbose), algoParams.Verbose = 1; end
    if isempty(algoParams.Visualize), algoParams.Visualize = 1; end
    if isempty(algoParams.species),
      algoParams.species(1).name = 'water'; % Water
      algoParams.species(1).frequency = DefaultWaterPpm; 
      algoParams.species(1).relAmps = DefaultWaterRelAmps;   
      algoParams.species(2).name = 'fat'; % Fat
      algoParams.species(2).frequency = DefaultFatPpm;
      algoParams.species(2).relAmps = DefaultFatRelAmps;
    end
    if numel(algoParams.species)~=NumChemSpecies,
      error('Current implementation handles %d species only.',NumChemSpecies);
    end
    for SpeciesNum = 1:numel(algoParams.species),
      if ~isfield(algoParams.species(SpeciesNum),'name'),
        error('name field does not exist for algoParams.species(%d).',SpeciesNum);
      end
    end; clear SpeciesNum;
    if     strcmpi(algoParams.species(1).name,'water'), wateridx = 1;
    elseif strcmpi(algoParams.species(2).name,'water'), wateridx = 2;
    else   error('"water" label not found in chemical species.');
    end
    if    strcmpi(algoParams.species(1).name,'fat'), fatidx = 1;
    elseif strcmpi(algoParams.species(2).name,'fat'), fatidx = 2;
    else   error('"fat" label not found in chemical species.');
    end
    if ~isfield(algoParams.species(wateridx),'frequency'),
      algoParams.species(wateridx).frequency = DefaultWaterPpm; 
      algoParams.species(wateridx).relAmps = DefaultWaterRelAmps;       
    end
    if ~isfield(algoParams.species(wateridx),'relAmps'),
      if numel(algoParams.species(wateridx).frequency)==1 || ~isequal(algoParams.species(wateridx).frequency, DefaultWaterPpm)
        algoParams.species(wateridx).relAmps = ones(size(algoParams.species(wateridx).frequency));
      else
        algoParams.species(wateridx).relAmps = DefaultWaterRelAmps;
      end
    end
    if numel(algoParams.species(wateridx).frequency)~=1,
      error('This implementation only expects one frequency for water in algoParams.species(%d).',wateridx);
    end
    if ~isfield(algoParams.species(fatidx),'frequency'),
      algoParams.species(fatidx).frequency = DefaultFatPpm; 
      algoParams.species(fatidx).relAmps = DefaultFatRelAmps;       
    end
    if ~isfield(algoParams.species(fatidx),'relAmps'),
      if numel(algoParams.species(fatidx).frequency)==1 || ~isequal(algoParams.species(fatidx).frequency, DefaultFatPpm)
        algoParams.species(fatidx).relAmps = ones(size(algoParams.species(fatidx).frequency));
      else
        algoParams.species(fatidx).relAmps = DefaultFatRelAmps;
      end
    end
    for SpeciesNum = 1:numel(algoParams.species),
      if numel(algoParams.species(SpeciesNum).frequency)~=numel(algoParams.species(SpeciesNum).relAmps)
        error('Number of frequencies and relative amplitudes do not match in algoParams.species(%d).',SpeciesNum);
      end
      if ~any(algoParams.species(SpeciesNum).relAmps~=0),
        error('Relative amplitudes should not be all zero for algoParams.species(%d).',SpeciesNum);
      end
    end; clear SpeciesNum;
    %--------------------------------------------------------------
    
    %--------------------------------------------------------------
    % Image parameters
    if ~isfield(imDataParams,'TE'),
        if isfield(imDataParams,'TEs'),
            imDataParams.TE = imDataParams.TEs;
            rmfield(imDataParams,'TEs');
        else
            error('TE field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'FieldStrength'),
        if isfield(imDataParams,'fieldStrength'),
            imDataParams.FieldStrength = imDataParams.fieldStrength;
            rmfield(imDataParams,'fieldStrength');
        else
            error('FieldStrength field is missing in imDataParams');
        end
    end
    if numel(imDataParams.TE)~=NumIDEALechoes,
      error('Current implementation handles %d echoes only.',NumIDEALechoes); 
    end
    if (imDataParams.TE(1)>imDataParams.TE(2) || imDataParams.TE(2)>imDataParams.TE(3)),
      error('Current implementation expects the %d TE to be increasing.',NumIDEALechoes);
    end
    gap12 = imDataParams.TE(2)-imDataParams.TE(1);
    gap23 = imDataParams.TE(3)-imDataParams.TE(2);
    if max(gap12,gap23)>min(gap12,gap23)*1.1,
      error('Current implementation expects the %d TE to be equally spaced.',NumIDEALechoes);
    end; clear gap12 gap23;
    if ~isfield(imDataParams,'water'),
        if isfield(imDataParams,'watermap'),
            imDataParams.water = imDataParams.watermap;
            rmfield(imDataParams,'watermap');
        else
            error('water field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'fat'),
        if isfield(imDataParams,'fatmap'),
            imDataParams.fat = imDataParams.fatmap;
            rmfield(imDataParams,'fatmap');
        else
            error('fat field is missing in imDataParams');
        end
    end
    if ~isfield(imDataParams,'fiterror'),
        if isfield(imDataParams,'fiterrormap'),
            imDataParams.fiterror = imDataParams.fiterrormap;
            rmfield(imDataParams,'fiterrormap');
        elseif isfield(imDataParams,'errormap'),
            imDataParams.fiterror = imDataParams.errormap;
            rmfield(imDataParams,'errormap');
        else
            imDataParams.fiterror = zeros(size(imDataParams.watermap));
        end
    end
    if ~isequal(size(imDataParams.water),size(imDataParams.fat)),
        error('water and fat maps have different matrix sizes.');
    end
    if ~isequal(size(imDataParams.water),size(imDataParams.fiterror)),
        error('fiterror map has a different matrix size.');
    end
    if ~isfield(imDataParams,'WaterFatPpmDiff'), imDataParams.WaterFatPpmDiff = []; end
    if isempty(imDataParams.WaterFatPpmDiff), imDataParams.WaterFatPpmDiff = 3.4; end
    imDataParams.WaterFatPpmDiff = abs(imDataParams.WaterFatPpmDiff);
    %--------------------------------------------------------------

    WaterPpm = algoParams.species(wateridx).frequency;
    WaterRelAmps = algoParams.species(wateridx).relAmps;
    FatPpm = algoParams.species(fatidx).frequency;
    FatRelAmps = algoParams.species(fatidx).relAmps;
    
    LarmorFreq = imDataParams.FieldStrength*GyromagneticRatio;
    blankppmInit = algoParams.species(wateridx).frequency(1) + imDataParams.WaterFatPpmDiff; clear maxval maxidx;
    
    % Full signal model for common R2
    % SignalAtTE = exp(-imDataParams.TE(:)*R2star) .* ... % T2 decay
    %  ( exp(-j*2*pi*imDataParams.TE(:)*LarmorFreq*-WaterPpm(:).')*WaterRelAmps(:) ...  % Water signal
    %  + exp(-j*2*pi*imDataParams.TE(:)*LarmorFreq*  -FatPpm(:).')*  FatRelAmps(:));    % Fat signal
    
    % Check where the signal blank is for water- or fat-only signal
    R2star = 1./(max(imDataParams.TE(:))*2);       % T2* set to twice longest echo time
    WaterSignalAtTE = exp(-j*2*pi*imDataParams.TE(:)*LarmorFreq*-WaterPpm(:).')*WaterRelAmps(:); % Water signal (no T2* decay)
    FatSignalAtTE   = exp(-j*2*pi*imDataParams.TE(:)*LarmorFreq*  -FatPpm(:).')*  FatRelAmps(:); % Fat signal (no T2* decay)
    WaterShortT2starSignalAtTE = exp(-imDataParams.TE(:)*R2star) .* WaterSignalAtTE;
    FatShortT2starSignalAtTE   = exp(-imDataParams.TE(:)*R2star) .* FatSignalAtTE;
    
    WaterBlankPpm            = fminsearch(@(blankppm) abs((           WaterSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-blankppm))/3).^2, blankppmInit);
    FatBlankPpm              = fminsearch(@(blankppm) abs((             FatSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-blankppm))/3).^2, blankppmInit);
    WaterShortT2starBlankPpm = fminsearch(@(blankppm) abs((WaterShortT2starSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-blankppm))/3).^2, blankppmInit);
    FatShortT2starBlankPpm   = fminsearch(@(blankppm) abs((  FatShortT2starSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-blankppm))/3).^2, blankppmInit);

    % Visualize results
    if algoParams.Visualize,
      ppmrange = (-5:0.05:5);
      fig = figure('name','Single- to multi-peak conversion');
      h=subplot(2,1,1,'Parent',fig); 
      s = (FatShortT2starSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-ppmrange))/3;
      plot(h,ppmrange,abs(s),'r:');
      hold(h,'on');
      s = (WaterShortT2starSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-ppmrange))/3;
      plot(h,ppmrange,abs(s),'b:'); 
      s = (FatSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-ppmrange))/3;
      plot(h,ppmrange,abs(s),'r-');
      s = (WaterSignalAtTE(:).'*exp(+j*2*pi*imDataParams.TE(:)*LarmorFreq*-ppmrange))/3;
      plot(h,ppmrange,abs(s),'b-'); 
      
      axis(h,'tight'); hold(h,'on');
      yaxisval = get(h,'YLim');
      set(h,'YLim',[0,yaxisval(2)],'xDir','reverse');
      xlabel(h,'ppm');
      plot(h,[FatShortT2starBlankPpm  ,FatShortT2starBlankPpm],[0,yaxisval(2)],'r:', ...
             [WaterShortT2starBlankPpm,WaterShortT2starBlankPpm],[0,yaxisval(2)],'b:',...
             [FatBlankPpm  ,FatBlankPpm],[0,yaxisval(2)],'r-', ...
             [WaterBlankPpm,WaterBlankPpm],[0,yaxisval(2)],'b-');
      hold(h,'off');
      clear h s yaxisval ppmrange
    end
    %clear WaterSignalAtTE FatSignalAtTE WaterShortT2starSignalAtTE FatShortT2starSignalAtTE;
    
    ppmRange = max(abs([FatBlankPpm, WaterShortT2starBlankPpm, FatShortT2starBlankPpm]-WaterBlankPpm));
    if ppmRange>MaxPpmRange,
      error('The ppm difference between the single and multipeak models\nis too large for conversion (%.3f)',ppmRange);
    end
    if algoParams.Verbose, 
        fprintf('Max frequency difference between the single and multipeak models is %.3f ppm\n',ppmRange);
    end
    clear ppmRange;
    
    % Conversion matrix (ignoring T2 decay)
    %   Multi-peak:    SignalsAtTE =  [FatSignalAtTE(:),WaterSignalAtTE(:)] [MP_Water] 
    %                                                                       [MP_Fat  ]
    %                              =  Matrix MP
    %
    %   Single-peak:   SignalsAtTE =  FT [SP_Water]
    %                                    [SP_Fat  ]
    %                                    [SP_Blank]
    %                              =  FT SP
    %  MP = (pinv(Matrix) FT) SP        
    tmp = sin(pi/3)*1i; %ft = [-0.5-tmp,1,-0.5+tmp; 1,1,1; -0.5+tmp,1,-0.5-tmp]; clear tmp;
    %   water    fat       blank  (water is at zero frequency, fat is positive freq)
    ft = [1,  -0.5+tmp, -0.5-tmp; 
          1,         1,        1;
          1,  -0.5-tmp, -0.5+tmp];
    clear tmp;
    outParams.SPtoMPmatrix = pinv( [WaterSignalAtTE(:),FatSignalAtTE(:)] ) * ft; clear tmp;
    
    if algoParams.Visualize,
      h=subplot(2,1,2,'Parent',fig);
      imagesc(abs(outParams.SPtoMPmatrix),'Parent',h); axis image;
      set(h,'XTick',[1,2,3],'XTickLabel',{'Water','Fat','Blank'},'YTick',[1,2],'YTickLabel',{'Water','Fat'});
      xlabel(h,'Single peak');
      ylabel(h,'Multi-peak');
      title(h,'Conversion matrix');
      clear h;
    end
    
    outParams.species = algoParams.species;
    tmp = outParams.SPtoMPmatrix*[imDataParams.water(:).';imDataParams.fat(:).';imDataParams.fiterror(:).'];
    outParams.species(wateridx).amps = reshape(tmp(1,:),size(imDataParams.water));
    outParams.species(fatidx  ).amps = reshape(tmp(2,:),size(imDataParams.fat  ));
    clear tmp;
    if isfield(imDataParams,'FieldStrength'), outParams.FieldStrength = imDataParams.FieldStrength; end
    if isfield(imDataParams,'TE'           ), outParams.TE = imDataParams.TE; end
    if isfield(imDataParams,'r2star'       ), outParams.r2star = imDataParams.r2star; end
    if isfield(imDataParams,'phasemap'     ), outParams.phasemap = imDataParams.phasemap; end
    if isfield(imDataParams,'fieldmap'     ), outParams.fieldmap = imDataParams.fieldmap; end
end