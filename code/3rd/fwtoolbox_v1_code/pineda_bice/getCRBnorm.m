function b = getCRBnorm(TE,wm,fm,wp,fp,df,psi,r2s)

rho = [wm;fm];
thetas = TE.*df*2*pi; % DH* Added 2PI factor

% make sure TE (& thetas) is a column vector
if(size(TE,1) < size(TE,2))
    TE = TE.';
    thetas = thetas.';
end

if(nargin < 7) %2-pt model (no field map (or R2*) supplied)
    
    b = ones(4,1);
    b(1:2,1) = 1./length(TE);
    b(3:4,1) = 1./(length(TE)*rho.^2);

elseif(nargin < 8) %IDEAL model (no R2* supplied)
    b = ones(5,1);
    b(1:2,1) = 1./length(TE);
    b(3:4,1) = 1./(length(TE)*rho.^2);
    b(5,1) = 1./((2 * pi)^2 * sum(TE.^2.*(sum(rho.^2)+2*wm*fm*cos(wp-fp-thetas))));
    % is the same as: 
    %s = make_signals_c(thetas,wm,fm,wp,fp,df,psi);
     %b(5,1) = 1./(sum(TE.^2.*(s.*conj(s))));

else %IDEAL-T2* model (R2* given)

    b = ones(6,1);
    b(1:2,1) = 1./(sum(exp(-2*TE*r2s)));
    b(3:4,1) = 1./(rho.^2*sum(exp(-2*TE*r2s)));
    b(5,1) = 1./((2 * pi)^2 * sum(TE.^2.*exp(-2*TE*r2s).*(wm^2+fm^2+2*wm*fm*cos(wp-fp-thetas))));
    b(6,1) = 1./(sum(TE.^2.*exp(-2*TE*r2s).*(wm^2+fm^2+2*wm*fm*cos(wp-fp-thetas))));
    % is the same as:
    %s = make_signals_c(thetas,wm,fm,wp,fp,df,psi,r2s);
    %b(5:6,1) = 1./(sum(TE.^2.*(s.*conj(s))));
    
end

end