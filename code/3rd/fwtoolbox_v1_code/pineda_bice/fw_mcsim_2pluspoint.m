function outParamsMC = fw_mcsim_2pluspoint(imDataParams, algoParams, simParams, N,SNR,INR)

if(nargin < 6) INR = 500; end % INR: Independent noise realizations per fat/water ratio
if(nargin < 5) SNR = 200; end
if(nargin < 4) N = 51; end

gyro = 42.58;
deltaF = [0 ; gyro*(algoParams.species(2).frequency(:) - algoParams.species(1).frequency)*(imDataParams.FieldStrength)];
% Only single-peak for now, so take the max fat peak for deltaF:
[maxAmp,imaxAmp] = max(abs(algoParams.species(2).relAmps(:)));
dF = gyro*(algoParams.species(2).frequency(imaxAmp)- algoParams.species(1).frequency)*(imDataParams.FieldStrength)
%dF = max(abs(deltaF));

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
        disp('No supported model for R2* without field map. Ignoring given R2*');
    end
else
    disp('No ''r2star'' parameter found - using model without R2*.');
end

% override to test certain thetas
% imDataParams.TE = [0 pi]./dF;

if(nargin < 7) %specify other signal properties for CRB calculation
    fatPhase = pi/4; 
    waterPhase = pi/4;
end

%number of statistics depends on model
if(isfield(simParams,'fieldmap') && isfield(simParams,'r2star'))
    numstat = 6;
elseif(isfield(simParams,'fieldmap'))
    numstat = 5;
else
    numstat = 4;
end

NSA = zeros(numstat,N);
NSA_0 = zeros(numstat,N);
%%MSE = zeros(numstat,N);
%%VAR = zeros(numstat,N);
ENkv = 10.^(linspace(-2,2,N));
FF = zeros(1,N);
for i = 1:N
	% generate water and fat signals
	waterAmp = 1/(1+ENkv(i));
    fatAmp = ENkv(i)/(1+ENkv(i))+1e-6;
    %fat fraction
    FF(i) = fatAmp/(waterAmp+fatAmp);
    if(isfield(simParams,'fieldmap') && isfield(simParams,'r2star')) %IDEAL-T2* model (R2* given)
        % get theoretical NSA
        NSA_0(:,i) = length(imDataParams.TE)* ...
                (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star)./ ...
                 getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star));
        % INR independent noise realizations for each ratio
        rhos = zeros(numstat,INR);
        for k = 1:INR
            % generate random gaussian noise
            noise = randn(2*length(imDataParams.TE),1) * (1/SNR);
            %noise=zeros(2*length(thetas),1);

            % create actual signal
            s = buildSignals(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,noise,fieldmap,r2star);
            % solve for x with least squares (nonlinear)
            rho = LSQsolve(imDataParams.TE,dF,s,...
                [waterAmp+1e-2;fatAmp+1e-2;waterPhase+1e-2;fatPhase+1e-2;fieldmap+1e-2;r2star+1e-2]);
            rhos(:,k) = rho;
        end
        % get monte carlo variances
        mc_var = var(rhos,0,2);
        NSA(:,i) = (1/SNR)^2.*length(imDataParams.TE)*...
            getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,dF,fieldmap,r2star)./mc_var;
        %%VAR(:,i) = mc_var;
        % get mean squared error
        %%MSE(:,i) = (1/INR)*sum((([waterAmp;fatAmp;waterPhase;fatPhase;fieldmap;r2star]*ones(1,INR))-rhos).^2,2);
        
    elseif(isfield(simParams,'fieldmap')) %IDEAL model (no R2* supplied)
        % get theoretical NSA
        NSA_0(:,i) = length(imDataParams.TE)* ...
                (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)),fieldmap)./ ...
                 getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)),fieldmap));
        % INR independent noise realizations for each ratio
        rhos = zeros(numstat,INR);
        for k = 1:INR
            % generate random gaussian noise
            noise = randn(2*length(imDataParams.TE),1) * (1/SNR);
            %noise=zeros(2*length(thetas),1);

            % create actual signal
            s = buildSignals(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)),noise,fieldmap);
            % solve for x with least squares (nonlinear)
            rho = LSQsolve(imDataParams.TE,abs(deltaF(2)),s,...
                  [waterAmp+1e-2;fatAmp+1e-2;waterPhase+1e-2;fatPhase+1e-2;fieldmap+1e-2]);
            rhos(:,k) = rho;
        end
        % get monte carlo variances
        mc_var = var(rhos,0,2);
        NSA(:,i) = (1/SNR)^2.*length(imDataParams.TE)*...
            getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)),fieldmap)./mc_var;
        %%VAR(:,i) = mc_var;
        % get mean squared error
        %%MSE(:,i) = (1/INR)*sum((([waterAmp;fatAmp;waterPhase;fatPhase;fieldmap]*ones(1,INR))-rhos).^2,2);
        
    else
        % get theoretical NSA
        NSA_0(:,i) = length(imDataParams.TE)* ...
                (getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)))./ ...
                 getCRB(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2))));
        % INR independent noise realizations for each ratio
        rhos = zeros(numstat,INR);
        for k = 1:INR
            % generate random gaussian noise
            noise = randn(2*length(imDataParams.TE),1) * (1/SNR);
            %noise=zeros(2*length(thetas),1);

            % create actual signal
            s = buildSignals(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)),noise);
            % solve for x with least squares (nonlinear)
            rho = LSQsolve(imDataParams.TE,abs(deltaF(2)),s,...
                [waterAmp+1e-2;fatAmp+1e-2;waterPhase+1e-2;fatPhase+1e-2],0);
            rhos(:,k) = rho;
        end
        % get monte carlo variances
        mc_var = var(rhos,0,2);
        NSA(:,i) = (1/SNR)^2.*length(imDataParams.TE)*...
            getCRBnorm(imDataParams.TE,waterAmp,fatAmp,waterPhase,fatPhase,abs(deltaF(2)))./mc_var;
        %%VAR(:,i) = mc_var;
        % get mean squared error
        %%MSE(:,i) = (1/INR)*sum((([waterAmp;fatAmp;waterPhase;fatPhase]*ones(1,INR))-rhos).^2,2);

    end
end

%put results into struct for output
outParamsMC.fatfraction = FF;
outParamsMC.nsaTheoretical = NSA_0;
outParamsMC.nsaSimulated = NSA;

% plot nsa by fat fraction
ld = ['\rho_w';'\rho_f';'\phi_w';'\phi_f';'  \psi';' R_2^*'];
subrow = ceil((numstat)/2);
subcol = 2;
figure;
for l = 1:numstat
    subplot(subrow,subcol,l)
    % plot theoretical
    semilogx(ENkv,outParamsMC.nsaTheoretical(l,:),'-');
    %plot(outParamsMC.fatfraction,outParamsMC.nsaTheoretical(l,:),'-');
    hold on;
    % plot lsq
    semilogx(ENkv,outParamsMC.nsaSimulated(l,:),'o');
    %plot(outParamsMC.fatfraction,outParamsMC.nsaSimulated(l,:),'o');
    hold off;
    
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

function rho = LSQsolve(echoes,df,s,x0)

rho = x0;

if(length(x0) == 6) %IDEAL-T2* model (R2* given)
    for j = 1:10
    wp = rho(3);
    fp = rho(4);
    psi = rho(5);
    r2s = rho(6);

    for i=1:length(echoes)
        A(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                        [cos(wp+psi*echoes(i)*2*pi), cos(fp+(df+psi)*echoes(i)*2*pi); ...
                         sin(wp+psi*echoes(i)*2*pi), sin(fp+(df+psi)*echoes(i)*2*pi)];

        A_wp(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                           [-sin(wp+psi*echoes(i)*2*pi), 0; ...
                             cos(wp+psi*echoes(i)*2*pi), 0];

        A_fp(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                           [0,-sin(fp+(df+psi)*echoes(i)*2*pi); ...
                            0, cos(fp+(df+psi)*echoes(i)*2*pi)];

        A_psi(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                            [-echoes(i)*2*pi*sin(wp+psi*echoes(i)*2*pi),-echoes(i)*2*pi*sin(fp+(df+psi)*echoes(i)*2*pi); ...
                              echoes(i)*2*pi*cos(wp+psi*echoes(i)*2*pi), echoes(i)*2*pi*cos(fp+(df+psi)*echoes(i)*2*pi)];

        A_r2s(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                        [-echoes(i)*cos(wp+psi*echoes(i)*2*pi), -echoes(i)*cos(fp+(df+psi)*echoes(i)*2*pi); ...
                         -echoes(i)*sin(wp+psi*echoes(i)*2*pi), -echoes(i)*sin(fp+(df+psi)*echoes(i)*2*pi)];

    end

    g = A*[rho(1);rho(2)]-s;
    g_prime = zeros(2*length(echoes),6);
    g_prime(:,1) = A(:,1);
    g_prime(:,2) = A(:,2);
    g_prime(:,3) = A_wp*[rho(1);rho(2)];
    g_prime(:,4) = A_fp*[rho(1);rho(2)];
    g_prime(:,5) = A_psi*[rho(1);rho(2)];
    g_prime(:,6) = A_r2s*[rho(1);rho(2)];

    %    rho = rho + (g_prime.'*g_prime)\(g_prime.'*(-g));
    
    if sum(isnan(g_prime(:)))>0 |  sum(isinf(g_prime(:)))>0
      g_prime
      wfamp = rho(1:2)
      wp
      fp
      psi
      r2s
      echoes
      x0
      s
    end

    rho = rho + pinv(g_prime.'*g_prime)*(g_prime.'*(-g));
    end

elseif(length(x0) == 5) %IDEAL model (no R2* supplied)
    for j = 1:10
    wp = rho(3);
    fp = rho(4);
    psi = rho(5);

    for i=1:length(echoes)
        A(2*i-1:2*i,:) = [cos(wp+psi*echoes(i)*2*pi), cos(fp+(df+psi)*echoes(i)*2*pi); ...
                          sin(wp+psi*echoes(i)*2*pi), sin(fp+(df+psi)*echoes(i)*2*pi)];

        A_wp(2*i-1:2*i,:) = [-sin(wp+psi*echoes(i)*2*pi), 0; ...
                              cos(wp+psi*echoes(i)*2*pi), 0];

        A_fp(2*i-1:2*i,:) = [0,-sin(fp+(df+psi)*echoes(i)*2*pi); ...
                             0, cos(fp+(df+psi)*echoes(i)*2*pi)];

        A_psi(2*i-1:2*i,:) = [-echoes(i)*2*pi*sin(wp+psi*echoes(i)*2*pi),-echoes(i)*2*pi*sin(fp+(df+psi)*echoes(i)*2*pi); ...
                               echoes(i)*2*pi*cos(wp+psi*echoes(i)*2*pi), echoes(i)*2*pi*cos(fp+(df+psi)*echoes(i)*2*pi)];

    end
    g = A*[rho(1);rho(2)]-s;
    g_prime = zeros(2*length(echoes),5);
    g_prime(:,1) = A(:,1);
    g_prime(:,2) = A(:,2);
    g_prime(:,3) = A_wp*[rho(1);rho(2)];
    g_prime(:,4) = A_fp*[rho(1);rho(2)];
    g_prime(:,5) = A_psi*[rho(1);rho(2)];

    rho = rho + (g_prime.'*g_prime)\(g_prime.'*(-g));
    end
    
elseif(length(x0) == 4) %Dixon
    for j = 1:10
    wp = rho(3);
    fp = rho(4);

    for i=1:length(echoes)
        A(2*i-1:2*i,:) = [cos(wp), cos(fp+df*echoes(i)*2*pi); ...
                          sin(wp), sin(fp+df*echoes(i)*2*pi)];

        A_wp(2*i-1:2*i,:) = [-sin(wp), 0; ...
                              cos(wp), 0];

        A_fp(2*i-1:2*i,:) = [0,-sin(fp+df*echoes(i)*2*pi); ...
                             0, cos(fp+df*echoes(i)*2*pi)];

    end
    g = A*[rho(1);rho(2)]-s;
    g_prime = zeros(2*length(echoes),4);
    g_prime(:,1) = A(:,1);
    g_prime(:,2) = A(:,2);
    g_prime(:,3) = A_wp*[rho(1);rho(2)];
    g_prime(:,4) = A_fp*[rho(1);rho(2)];

    rho = rho + (g_prime.'*g_prime)\(g_prime.'*(-g));
    end
    
else
    return
end

end
