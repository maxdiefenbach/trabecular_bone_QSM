function [subset] = getsubset(input, setmin, setmax)

subset = [];
if( isempty(input) )
  subset = setmin:setmax;  
else    
  if( input(1)<setmin )
    subset = setmin:setmax;
  else
    for ss=1:numel(input),
      if( input(ss)>=setmin && input(ss)<=setmax )
          subset = [subset input(ss)];
      end
    end
  end
end

return;