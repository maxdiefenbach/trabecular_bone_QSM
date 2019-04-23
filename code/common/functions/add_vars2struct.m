function Struct = add_vars2struct(Struct, varargin)
% adds variables to Struct

    for iArgin = 1:numel(varargin)
        
        varName = inputname(1 + iArgin);
        
        Struct.(varName) = varargin{iArgin};
    
    end

end