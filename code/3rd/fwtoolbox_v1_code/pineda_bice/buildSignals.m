function s = buildSignals(echoes,wm,fm,wp,fp,df,noise,psi,r2s)

if(nargin < 7) noise = zeros(2*length(thetas),1); end
if(nargin < 6) df = 210; end
if(nargin < 5) fp = pi/4; end
if(nargin < 4) wp = pi/4; end

rho = [wm;fm];

if(nargin < 8) %2-pt model (no field map (or R2*) supplied)
    
    for i=1:length(echoes)
        A(2*i-1:2*i,:) = [cos(wp), cos(fp+df*echoes(i)*2*pi); ...
                          sin(wp), sin(fp+df*echoes(i)*2*pi)];
        
    end

elseif(nargin < 9) %IDEAL model (no R2* supplied)
    
    for i=1:length(echoes)
        A(2*i-1:2*i,:) = [cos(wp+psi*echoes(i)*2*pi), cos(fp+(df+psi)*echoes(i)*2*pi); ...
                          sin(wp+psi*echoes(i)*2*pi), sin(fp+(df+psi)*echoes(i)*2*pi)];
        
    end

else %IDEAL-T2* model (R2* given)
    
    for i=1:length(echoes)
        A(2*i-1:2*i,:) = exp(-echoes(i)*r2s).* ...
                        [cos(wp+psi*echoes(i)*2*pi), cos(fp+(df+psi)*echoes(i)*2*pi); ...
                         sin(wp+psi*echoes(i)*2*pi), sin(fp+(df+psi)*echoes(i)*2*pi)];
    end
    
end

s = A*rho + noise;

end