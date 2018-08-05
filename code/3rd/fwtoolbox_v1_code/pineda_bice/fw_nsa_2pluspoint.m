function outParamsNSA = fw_nsa_2pluspoint(imDataParams, algoParams, simParams, N)
%simParams - a structure expecting a value for "fieldmap" and "r2star"
%to use models without these parameters, leave value out of simParams, as
%a 0 value will use a model with that 0 value
%any misnamed inputs will be ignored
%note that R2* model must include field map

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency)*(imDataParams.FieldStrength)];

% Only single-peak supported (for now), so take the max fat peak for deltaF:
[maxAmp,imaxAmp] = max(abs(algoParams.species(2).relAmps(:)));
dF = gyro*(algoParams.species(2).frequency(imaxAmp)- algoParams.species(1).frequency)*(imDataParams.FieldStrength)
%dF = max(abs(deltaF))

%get fieldmap and r2star, if given
if(nargin < 3) simParams = []; end
if isfield(simParams,'fieldmap')
    fieldmap = simParams.fieldmap;
else
    disp('No ''fieldmap'' parameter found - using model without field map.');
end
if isfield(simParams,'r2star')
    r2star = simParams.r2star;
    if ~isfield(simParams,'fieldmap')
        disp('No supported model for R2* without field map - ignoring given R2*');
    end
else
    disp('No ''r2star'' parameter found - using model without R2*.');
end

% override to test certain thetas
%imDataParams.TE = [0 pi]./dF;
%imDataParams.TE = [-2*pi/3 0 2*pi/3]./dF;
%imDataParams.TE = [-pi/6 pi/2 7*pi/6]./dF;
%imDataParams.TE = [0 pi 2*pi]./dF;

if(nargin < 5) %specify other signal properties for CRB calculation
    fatPhase = pi/4; 
    waterPhase = pi/4;
    
    if(nargin < 4)
        N = 51;
    end
end

%number of statistics depends on model
pc = 0; %flag to use phase-constrained model
if(isfield(simParams,'fieldmap') && isfield(simParams,'r2star'))
    numstat = 6;
elseif(isfield(simParams,'fieldmap'))
    numstat = 5;
else
    numstat = 4;
end
NSA = zeros(numstat,N);

%generate N fat/water ratios to test
ENkv = 10.^(linspace(-2,2,N));
FF = zeros(1,N);
for k = 1:N
    waterAmp = 1/(1+ENkv(k));
    fatAmp = ENkv(k)/(1+ENkv(k))+1e-6;
    %fat fraction
    FF(k) = fatAmp/(waterAmp+fatAmp);
    
    if(isfield(simParams,'fieldmap') && isfield(simParams,'r2star')) %IDEAL-T2* model (R2* given)
        
        NSA(:,k) = length(imDataParams.TE)* ...
            (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star)./ ...
             getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star));
         
    elseif(isfield(simParams,'fieldmap')) %IDEAL model (no R2* supplied)
      
        NSA(:,k) = length(imDataParams.TE)* ...
            (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap)./ ...
             getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap));
        
    else %Dixon model
        NSA(:,k) = length(imDataParams.TE)* ...
            (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF)./ ...
             getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF));

    end

end

%put results into struct for output
outParamsNSA.fatfraction = FF;
outParamsNSA.nsa = NSA;

% plot nsa by fat fraction
ld = ['\rho_w';'\rho_f';'\phi_w';'\phi_f';'  \psi';' R_2^*'];
subrow = ceil((numstat)/2);
subcol = 2;
figure;
for l = 1:numstat
    subplot(subrow,subcol,l)
    % plot theoretical
    semilogx(ENkv,outParamsNSA.nsa(l,:),'-');
    %plot(outParamsNSA.fatfraction,outParamsNSA.nsa(l,:),'-');
    
    axis([0 100 0 length(imDataParams.TE)+0.5])
    %axis([0 1 0 length(imDataParams.TE)+0.5])
    title(ld(l,:),'FontSize',14,'FontWeight','b');
    if(l > (numstat-2))
        xlabel('log_{10}(Fat/Water)','FontSize',14);
        %xlabel('Fat Fraction','FontSize',14); 
    end
    set(gca,'xtick',[1e-2 1e0 1e2],'FontSize',14);
    %set(gca,'xtick',[0 0.5 1],'FontSize',14);
    if(mod(l,subcol) == 1)
        ylabel('NSA','FontSize',14);
    end
    set(gca,'ytick',0:length(imDataParams.TE),'FontSize',14);
end

end