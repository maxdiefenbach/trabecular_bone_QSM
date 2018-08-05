%  x = nonlin_cg(params)
%
% Nonlinear conjugate gradient algorithm for minimizing the objective
% function
% f(b) = ||FTx -y ||_2^2 +  ||Wx||_1 + TV(x)
%
%
% Inputs: 
% params structure containing data and reconstruction parameters
% 
% params.data   % undersampled k-space data
% params.image  % initial image
% 
% params.FT     % undersampled Fourier operator
% params.TV     % finite differences operator
% params.WT     % wavelet transform operator
% 
%
% Outputs: reconstructed image x 
% 
% Based on
% M.Lustig, D.Donoho and J.Pauly, "Sparse MRI: The application of compressed sensing for
% rapid MR imaging" Magn Reson Med 2007 58: 1182-95 
%
%==========================================================================
function [x,f_x]= nonlin_cg (params)
%Initialization============================================================
x  = params.image;             % the initial image
r  = -gradient (x,params);     % steepest direction 
dx = r;                        % initial search direction
k  = params.alpha;             % alpha (0,0.5)the angle of the backtracking line search
h  = params.beta;              % beta (0,1) factor to set the step size in the line search
t0 = 1;

if params.debug
    f_x = zeros(params.n,1);
end
%Body======================================================================
    for i=1:params.n
             
    t = t0;
    [FTx, FTdx, Dx, Ddx, WTx, WTdx] = preobjective(x, dx, params);
    c1       = objective(FTx, FTdx, Dx, Ddx, WTx, WTdx, 0, params);  
    
    g        = (gradient(x,params));
    c2       = (k*(g(:)'*dx(:)));
    
    % Perform backtracking line search
    itcount = 0;
    while   (objective(FTx, FTdx, Dx, Ddx, WTx, WTdx, t, params)>= (c1 + (t*c2)))&&(itcount < 100)
        t=t*h;
        itcount = itcount +1;
    end
    
    
    % Adjust the initial step size for the next iteration
    
    if (itcount >2)
        t0 = t0*h;
    end
    
    if (itcount == 0)
        t0 = t0/h;
    end
    
    
    % Update image
    x     = x + t*dx;           
    

    r_new = -gradient (x, params);  %new residual
    beta  = (r_new(:)'*r_new(:))/(r(:)'*r(:));      %Fletcher Reeves
   % beta  = (r_new(:)'*(r_new(:)-r(:)))/(r(:)'*r(:));%Polak Riberie
   % beta = max(beta,0);
    
    % Update search direction 
    dx     = r_new + (beta*dx);
    r      = r_new;
    
            
    % track value of the objective function for debug purposes        
    if params.debug        
    f_x(i)=c1;
    end
    
    end
    
  
function [FTx, FTdx, Dx, Ddx, WTx, WTdx] = preobjective(x, dx, params)
% precalculates Fourier, wavelet and finite differences transforms to make line search cheap

FTx = params.FT*x;
FTdx = params.FT*dx;

if params.TVWeight
    Dx = params.TV*x;
    Ddx = params.TV*dx;
else
    Dx = 0;
    Ddx = 0;
end

if params.WTWeight    
    WTx  = params.WT*x;
    WTdx = params.WT*dx;
else
    WTx  = 0;
    WTdx = 0;
end


function res = objective(FTx, FTdx, Dx, Ddx, WTx, WTdx, t, params)
%computes the objective function

p = params.pNorm;

obj = FTx + t*FTdx - params.data;
obj = obj(:)'*obj(:);

if params.TVWeight
    w = Dx(:) + t*Ddx(:);    
    sum(Dx(:));
    TV = (w.*conj(w)+params.l1Smooth).^(p/2); 
else
    TV = 0;
end

if params.WTWeight
   w = WTx(:) + t*WTdx(:);   
   WT = (w.*conj(w)+params.l1Smooth).^(p/2);
else
    WT=0;
end


TV = sum(TV.*params.TVWeight(:));
WT = sum(WT.*params.WTWeight(:));

res = obj + (TV) + (WT);


function grad = gradient(x,params)
% computes the gradient at position x
gradWT = 0;
gradTV  = 0;

gradObj = gOBJ(x,params);
if params.WTWeight
gradWT = gWT(x,params);
end
if params.TVWeight
gradTV = gTV(x,params);
end

grad = (gradObj +  params.WTWeight.*gradWT + params.TVWeight.*gradTV);


function gradObj = gOBJ(x,params)
% computes the gradient of the data consistency

	gradObj = params.FT'*(params.FT*x - params.data);

gradObj = 2*gradObj;

function grad = gWT(x,params)
% compute gradient of the L1 transform operator

p = params.pNorm;

WTx = params.WT*x;

G = p*WTx.*(WTx.*conj(WTx)+params.l1Smooth).^(p/2-1);
grad = params.WT'*G; 


function grad = gTV(x,params)
% compute gradient of TV operator

p = params.pNorm;

Dx = params.TV*x;

G = p*Dx.*(Dx.*conj(Dx) + params.l1Smooth).^(p/2-1);
grad = params.TV'*G;

