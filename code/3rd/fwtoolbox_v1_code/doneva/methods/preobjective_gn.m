function [FD1w, FD1dw, FD1f, FD1df, WTw, WTdw, WTf, WTdf, FD2o, FD2do] = preobjective_gn(x,dx,param)
 % finite differences or wavelet constraint for water and fat, for now with the same constant   
 % second order finite differences constraint for offresonance
  

 if param.FD1Weight 
      % water part
      FD1w  = param.FD1*(x(:,:,1));
      FD1dw = param.FD1*dx(:,:,1);
      % fat part
      FD1f  = param.FD1*(x(:,:,2));
      FD1df = param.FD1*dx(:,:,2);
 else
     FD1w  = 0;
     FD1dw = 0;
     FD1f  = 0;
     FD1df = 0;
 end  
 
 if param.FD2Weight 
     % offresonance part
      FD2o  = param.FD2*(x(:,:,3));
      FD2do = param.FD2*dx(:,:,3);
     
 else
     FD2o  = 0;
     FD2do = 0;
 end  
 
 
 
 if param.WTWeight 
      % water part
      WTw   = param.WT*(x(:,:,1));
      WTdw  = param.WT*(dx(:,:,1));
      
      % fat part
      WTf   = param.WT*(x(:,:,2));
      WTdf  = param.WT*(dx(:,:,2));      
      
 else
     WTw  = 0;
     WTdw = 0;
     WTf  = 0;
     WTdf = 0;
 end  
 
 
