function res = DAMP(x, t, f_wf, rel_amp,mask, imsize)
%
%
res.adjoint = 0;
res.x       = x;
%res.Trajectory = Trajectory;
res.t       = t;
res.f_wf    = f_wf;
res.rel_amp = rel_amp;
res.mask    = mask;
res.imsize  = imsize;
%res.no_sp = no_sp;
%res.no_samples = no_samples; 
res       = class(res,'DAMP');

