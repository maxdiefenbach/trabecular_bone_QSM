function res = mtimes(a,b)


if a.adjoint
	res = D_t(b);

else
	res = D(b);

end




    
