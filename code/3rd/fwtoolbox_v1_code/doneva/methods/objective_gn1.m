function res = objective_gn1(x, dx, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t, param)
        
 res = (forward_op_mp1(x+t*dx, param.t, param.f_wf(2:end),param.rel_amp(2:end), param.mask) - param.k_data); 
 res = res(:)'*res(:);
 
 if param.FD1Weight  
     p   = param.pnorm_im;
     temp1 = FD1w + t*FD1dw;
     temp2 = FD1f + t*FD1df;
     FD1  = (temp1.*conj(temp1) + param.mu).^(p/2) + (temp2.*conj(temp2) + param.mu).^(p/2);
 else
     FD1=0;
 end
  

 if param.FD2Weight 
     p   = param.pnorm_off;
     temp = FD2o + t*FD2do;
     FD2  = (temp.*conj(temp) + param.mu).^(p/2);
 else
     FD2=0;
 end
 
 if param.WTWeight  
     p   = param.pnorm_im;
     temp1 = WTw + t*WTdw;
     temp2 = WTf + t*WTdf;
     WT    = (temp1.*conj(temp1) + param.mu).^(p/2) + (temp2.*conj(temp2) + param.mu).^(p/2);
 else
     WT=0;
 end
 
 tmp = t*dx;
 tmp = tmp(:)'*tmp(:);

 if param.debug
 disp('data consistency')
 disp(res)
 disp('fat water TV l1 term')
 disp(param.FD1Weight*sum(FD1(:)))
 disp('fat water wavelet l1 term')
 disp(param.WTWeight*sum(WT(:)))
 disp('offresonance term')
 disp(param.FD2Weight*sum(FD2(:)))
 end
 
 res    = param.dataWeight*res + param.alpha_n*tmp +  param.FD1Weight*sum(FD1(:)) + param.FD2Weight*sum(FD2(:))+ param.WTWeight*sum(WT(:));

