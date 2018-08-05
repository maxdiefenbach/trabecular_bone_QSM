function res = mtimes(a,b)

if a.adjoint
    % transposed operator
    res = dforward_opH(b, a.x, a.t, a.f_wf, a.rel_amp, a.mask);    
else    
    res = dforward_op(b, a.x, a.t, a.f_wf, a.rel_amp,a.mask);
end



