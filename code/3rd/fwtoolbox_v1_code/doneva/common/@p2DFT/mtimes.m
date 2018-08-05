function res = mtimes(a,b)
% res = mtimes(FT, x)
%

    % Conjugate FT operator
if a.adjoint
    bb = reshape(b,a.dataSize(1),a.dataSize(2));            % orders column vector b in a matrix
    res = zpad(bb.*a.mask,a.imSize(1),a.imSize(2));         % zero pads data matrix to image matrix
    res = ifft2c(res);                                      % ifft for zero centered 2d signals
    res = res.*conj(a.ph);                                  % multiplies by conjugate phase
    switch a.mode
    	case 0
		res = real(res);
        figure;imshow(res,[]);        
   	case 1
		res = real(res);
         figure;imshow(res,[]); 
    end



    % FT operator 
else
    bb = reshape(b,a.imSize(1),a.imSize(2));            
    
    switch a.mode
    	case 0
		bb = real(bb);
   	case 1
		bb = real(bb);
    end
    
    bb = bb.*a.ph; % phase correct
    res = fft2c(bb);
    res = crop(res,a.dataSize(1),a.dataSize(2));
    res = res.*a.mask;
end

if size(b,2) == 1
    res = res(:);
end



    
