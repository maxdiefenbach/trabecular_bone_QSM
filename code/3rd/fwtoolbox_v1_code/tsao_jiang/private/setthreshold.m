%%
% Build histogram of image and calculate noise threshold
%

% Jeffrey Tsao, 2011
function level = setthreshold(img)
bin = (0.5:1:100-1)/100*double(max(abs(img(:))));
while 1,
    %-- Gamma variate fit to histogram
    [binfreq,bin] = hist(single(abs(img(:))),bin);
    [maxval,maxidx] = max(binfreq);
    tmp_mean = sum(bin.*binfreq)./sum(binfreq); % mean = k*theta
    tmp_mode = double(bin(maxidx));             % mode = k*theta-theta
    clear maxval maxidx;
    x = double([tmp_mean/(tmp_mean-tmp_mode),(tmp_mean-tmp_mode)]);   % k, theta
    x = fminsearch(@(x) OptFunc(x,bin,binfreq), x, optimset('Display','off'));
    y = real(exp( (x(1)-1)*log(bin) - bin/x(2) - x(1)*log(x(2)) - gammaln(x(1)) ));
    cum_y = cumsum(y);
    clear x y;
    if cum_y(end)>0,
        cum_y = cum_y/max(cum_y);
        [tmpu,tmpidx]=unique(cum_y);
        level = interp1([0,cum_y(tmpidx)],[0,bin(tmpidx)],0.99);
        clear tmpu tmpidx;
        if level>bin(4), break; end % Done
        bin = bin/3;  % Need finer division
    else
        % Histogram fitting did not work
        level=min(abs(img(:))); % Just return the smallest intensity
        break;
    end
    clear tmp_mean tmp_mode;
end
clear bin binfreq;
end

%%
function errval = OptFunc(x,bin,binfreq)
x(2) = max(x(2),0.00000001);
y = exp( (x(1)-1)*log(bin(:)) - bin(:)/x(2) - x(1)*log(x(2)) - gammaln(x(1)) );
assignin('base','x',x);
coef = pinv(y)*double(binfreq(:));

fit_y = y*coef;
errval = sum(abs( fit_y - double(binfreq(:)) ));
end
