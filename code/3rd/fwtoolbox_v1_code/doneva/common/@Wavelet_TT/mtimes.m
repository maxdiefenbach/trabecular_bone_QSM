function res = mtimes(a,b)

if isa(a,'Wavelet_TT') == 0
    error('In  A.*B only A can be Wavelet operator');
end

if a.adjoint
	if(length(size(b))==3)
		res = zeros(size(b,1),size(b,2),size(b,3));
		for k=1:size(b,3)
			res(:,:,k) = IWT2_PO_TT2(real(b(:,:,k)),a.Lx,a.Ly,a.qmf) + 1i*IWT2_PO_TT2(imag(b(:,:,k)),a.Lx,a.Ly,a.qmf);
		end
	else
		res = IWT2_PO_TT2(real(b),a.Lx,a.Ly,a.qmf) + 1i*IWT2_PO_TT2(imag(b),a.Lx,a.Ly,a.qmf);
	end
else
	if(length(size(b))==3)
		res = zeros(size(b,1),size(b,2),size(b,3));
		for k=1:size(b,3)
			res(:,:,k) = FWT2_PO_TT2(real(b(:,:,k)),a.Lx,a.Ly,a.qmf) + 1i* FWT2_PO_TT2(imag(b(:,:,k)),a.Lx,a.Ly,a.qmf);
		end
	else
		res = FWT2_PO_TT2(real(b),a.Lx,a.Ly,a.qmf) + 1i* FWT2_PO_TT2(imag(b),a.Lx,a.Ly,a.qmf);
	end
end