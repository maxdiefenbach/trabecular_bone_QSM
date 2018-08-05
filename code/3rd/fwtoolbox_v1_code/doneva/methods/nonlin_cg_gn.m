function x = nonlin_cg_gn(x0, param)
%
%
%  nonlinear CG to minimize
%  f(xn,dx)  = ||DF(xn)dx - (y-F(xn)||l2^2 + alphan||dx||l2^2 +
%  lambda1||Psi(xn + dx||l1 +  lambda2||Phi(xn + dxn)||l1
% 
% Psi operates on the water and fat images,
% Phi operates on the offresonance only
% 
% need current point param.xn


% set backtracking line search parameters
alpha    = param.alpha;
beta     = param.beta; 
t0       = 1;

% initialization 
x        = x0;

resr = param.data - param.DA*x;
resr = resr(:)'*resr(:);
resr0 = resr; 


g        = gradient(x0,param);
dx       = -g;
gtg      = g(:)'*g(:);                  
n        = 0;


% iterations 
while (n<param.maxIter) && ( resr/resr0 > param.threshold)

[DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw, WTf, WTdf, FD2o, FD2do] = preobjective(x,dx,param); 
f0 = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, 0, param);  
% 
% 
c2 = -abs(g(:).'*dx(:));   

t  = t0;
m  = 0;
% first calculation of f1 
f1 = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t, param);  
% backtraching line search          
while(f1 > f0 + alpha*t/beta*c2)&&(m<param.maxIterLS)
    
    m  = m+1;
    t  = t*beta;  
    f1 = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t, param); 
    
end


% adjust initial step size for faster search
if m>2
    
    t0=t0*beta;
end

if m<1
    t0 = t0/beta;
end

% update image
  x = x + t*dx;
  
% update gradient
  g = gradient(x,param);

% update search direction
  
  c     = g(:)'*g(:);
  gamma = g(:)'*g(:)/gtg;
  gtg   = c;                      % save g transposed g for the next iteration
  
  dx = -g + gamma*dx;    
  n = n + 1;
  
  resr = param.data - param.DA*x;
  resr = resr(:)'*resr(:);
end



 function [DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw, WTf, WTdf, FD2o, FD2do] = preobjective(x,dx,param)
 % finite differences or wavelet constraint for water and fat, for now with the same constant   
 % second order finite differences constraint for offresonance
 
 DAx    = param.DA*x;
 DAdx   = param.DA*dx;
 

 if param.FD1Weight 
      % water part
      FD1w  = param.FD1*(x(:,:,1) + param.xn(:,:,1));
      FD1dw = param.FD1*dx(:,:,1);
      % fat part
      FD1f  = param.FD1*(x(:,:,2) + param.xn(:,:,2));
      FD1df = param.FD1*dx(:,:,2);
 else
     FD1w  = 0;
     FD1dw = 0;
     FD1f  = 0;
     FD1df = 0;
 end  
 
 if param.FD2Weight 
     % offresonance part
      FD2o  = param.FD2*(x(:,:,3) + param.xn(:,:,3));
      FD2do = param.FD2*dx(:,:,3);
     
 else
     FD2o  = 0;
     FD2do = 0;
 end  
 
 
 
 if param.WTWeight 
      % water part
      WTw   = param.WT*(x(:,:,1)+ param.xn(:,:,1));
      WTdw  = param.WT*(dx(:,:,1));
      
      % fat part
      WTf   = param.WT*(x(:,:,2)+ param.xn(:,:,2));
      WTdf  = param.WT*(dx(:,:,2));      
      
 else
     WTw  = 0;
     WTdw = 0;
     WTf  = 0;
     WTdf = 0;
 end  
 
 
 
 
function res = objective(x, dx, DAx, DAdx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t, param)
        
 res = (DAx + t*DAdx - param.data);
 res = res(:)'*res(:);
 
 
 
 if param.FD1Weight
     p     = param.pnorm_im;
     temp1 = FD1w + t*FD1dw;
     temp2 = FD1f + t*FD1df;
     FD1   = (temp1.*conj(temp1) + param.mu).^(p/2) + (temp2.*conj(temp2) + param.mu).^(p/2);
 else
     FD1   = 0;
 end
  

 if param.FD2Weight  
     p    = param.pnorm_off;
     temp = FD2o + t*FD2do;
     FD2  = (temp.*conj(temp) + param.mu).^(p/2);
 else
     FD2=0;
 end
 
 if param.WTWeight     
   p    = param.pnorm_im; % 120426: This line was missing
   temp1 = WTw + t*WTdw;
   temp2 = WTf + t*WTdf;
   WT    = (temp1.*conj(temp1) + param.mu).^(p/2) + (temp2.*conj(temp2) + param.mu).^(p/2);
 else
     WT=0;
 end
 
 tmp = x + t*dx;
 tmp = tmp(:)'*tmp(:);

 res    = param.dataWeight*res + param.alpha_n*tmp +  param.FD1Weight*sum(FD1(:)) + param.FD2Weight*sum(FD2(:))+ param.WTWeight*sum(WT(:));


% compute gradient 


function grad = gradient(x,param)
% compute total gradient
grad(:,:,1) = param.FD1Weight*gFD1(x(:,:,1) + param.xn(:,:,1), param);
grad(:,:,2) = param.FD1Weight*gFD1(x(:,:,2) + param.xn(:,:,2), param);
grad(:,:,3) = param.FD2Weight*gFD2(x(:,:,3) + param.xn(:,:,3), param);

if param.WTWeight
  grad(:,:,1) = grad(:,:,1) + param.WTWeight*gWT(x(:,:,1) + param.xn(:,:,1), param);
  grad(:,:,2) = grad(:,:,2) + param.WTWeight*gWT(x(:,:,2) + param.xn(:,:,2), param);
end
grad = grad + param.dataWeight*gData(x,param) + param.alpha_n*2*x; 


function grad = gData(x,param)
grad = param.dataWeight*2*(param.DA'*(param.DA*x - param.data));



function grad = gFD1(x,param)
% gradient of the second order finite differences term
p    = param.pnorm_im;
fd   = param.FD1*x;
grad = p*fd.*(fd.*conj(fd) + param.mu).^(p/2-1);
grad = param.FD1'*grad;


function grad = gFD2(x,param)
% gradient of the second order finite differences term
p    = param.pnorm_off;
fd   = param.FD2*x;
grad = p*fd.*(fd.*conj(fd) + param.mu).^(p/2-1);
grad = param.FD2'*grad;


function grad = gWT(x,param)
% gradient of the second order finite differences term
p    = param.pnorm_im;
wt   = param.WT*x;
grad = p*wt.*(wt.*conj(wt) + param.mu).^(p/2-1);
grad = param.WT'*grad;

 

