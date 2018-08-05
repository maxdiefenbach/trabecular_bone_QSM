
ii = sqrt(-1);
xx = ones(64,64,3);
mask = ones(64,64,3);
f_wf = [0.087 0.693 0.128 0.004 0.039 0.048];
t    = [ 1 2 3];
rel_amp = [0.087 0.693 0.128 0.004 0.039 0.048];




%x = rand(100,1)+ii*rand(100,1);
x = rand(64,64,3) + ii*rand(64,64,3);
y = rand(64,64,3) + ii*rand(64,64,3);


    
Ax  = dforward_op (x, xx, t, f_wf, rel_amp, mask);
Aty = dforward_opH(y, xx, t, f_wf, rel_amp, mask);


x(:)'*Aty(:)
Ax(:)'*y(:)


