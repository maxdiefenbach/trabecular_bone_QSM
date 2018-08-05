function res = mtimes(a,b)


if a.adjoint
	res = D2_t(b);

else
	res = D2(b);

end




    
