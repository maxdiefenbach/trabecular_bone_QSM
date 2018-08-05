function [w,f,fm] = wfcs(kdata, x0, cg_param)


%% normalize data
n_factor = 100/sqrt(kdata(:)'*kdata(:));
k_datan  = kdata*n_factor;

%% set initial values
xn = x0;
imsize = size(x0);
sx = imsize(1);
sy = imsize(2);

z  = zeros(size(xn));
cg_param.k_data    = k_datan;


for i = 1:cg_param.num_outer_iter_wfcs
    
    
    if cg_param.debug
       disp('Iteration'); disp(i);
    end
    
    if(mod(i,1)==0)
        cg_param.FD1Weight = cg_param.FD1Weight*0.9;
        cg_param.WTWeight  = cg_param.WTWeight*0.9;
    end
    
    
    res = k_datan - forward_op_mp1(xn, cg_param.t, cg_param.f_wf(2:end),cg_param.rel_amp(2:end), cg_param.mask);
    
    if cg_param.debug
        disp('Residuum ');
        disp(sqrt(res(:)'*res(:)/(k_datan(:)'*k_datan(:))));
    end
    
    % Update residuum, current value of x in the cg algorithm
    cg_param.data      = res;
    cg_param.xn        = xn;
    cg_param.DA        = DAMP(xn, cg_param.t, cg_param.f_wf(2:end), cg_param.rel_amp(2:end), cg_param.mask,[sx,sy]);
    
    for p = 1:cg_param.num_inner_iter_wfcs
        z                = nonlin_cg_gn(z, cg_param);
    end
    
    
    % z(:,:,3) = real(z(:,:,3));
    
    
    [FD1w, FD1dw, FD1f, FD1df, WTw, WTdw, WTf, WTdf, FD2o, FD2do] = preobjective_gn(xn,z,cg_param);
    f0 = objective_gn1(xn, z, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, 0, cg_param);
    t0 = 1;
    f1 = objective_gn1(xn, z, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t0, cg_param);
    
    if cg_param.debug
        disp('f0=');disp(f0);
        disp('f1=');disp(f1);
    end
    
    cc = cg_param.DA'*res;
    cc = -abs(cc(:).'*z(:));
    
    
    k  = 1;
    
    while(f1 > f0 + cg_param.alpha*t0/cg_param.beta*cc)&&(k<50)
        
        k  = k+1;
        t0  = t0*cg_param.beta;
        f1 = objective_gn1(xn, z, FD1w, FD1dw, FD1f, FD1df, WTw, WTdw,  WTf, WTdf, FD2o, FD2do, t0, cg_param);
        
    end
    
    
    %    update image
    
    xn = xn + t0*z;
    
    if cg_param.display
        figure(1); imshow([abs(xn(:,:,1)) abs(xn(:,:,2))],[]); drawnow;
        figure(2); imshow(real(xn(:,:,3)),[]); drawnow;
        % figure(3); imshow(imag(xn(:,:,3)),[]); drawnow;
    end
    
    
end

w  = xn(:,:,1)/n_factor;
f  = xn(:,:,2)/n_factor;
fm = xn(:,:,3);