function [ip, wp, fp, pp, rp] = init_pyramid(s, D, ht, wd)
ip = cell(1, D); wp = cell(1, D); 
fp = cell(1, D); pp = cell(1, D); rp = cell(1, D);
for d=1:D,
    if d==1,
        % the finest level
        ip{d} = s;        
    else        
        ip{d} = imresize(ip{d-1}, .5, 'nea'); 
    end
    wp{d} = zeros(round(ht/(2^(d-1))), round(wd/(2^(d-1)))); 
    fp{d} = wp{d}; pp{d} = wp{d}; 
    rp{d}= wp{d};

end
