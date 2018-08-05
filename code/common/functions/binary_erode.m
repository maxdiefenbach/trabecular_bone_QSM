function mask = binary_erode(mask, N)
    
    if nargin < 2
        N = 1;
    end

    mask = padarray(mask, [1, 1, 1]);

    SE = strel('sphere', 1);
    
    for n = 1:N
        mask = imerode(mask, SE);
    end
    
    mask = mask(2:end-1, 2:end-1, 2:end-1);

end