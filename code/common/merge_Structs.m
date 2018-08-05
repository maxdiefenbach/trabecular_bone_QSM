function StructMerged = merge_Structs(varargin)
    
    StructMerged = varargin{1};
    
    for iArg = 2:numel(varargin)
        StructMerged = merge_twoStructs(StructMerged, varargin{iArg});
    end

end


function StructMerged = merge_twoStructs(S1, S2, overwrite)
% add fields from S2 to S1
% do overwrite fields by default
    
    if nargin < 3
        overwrite = 1;
    end

    StructMerged = S1;

    if isempty(S2) | ~isstruct(S2)
        return
    else
        fieldNameCell = fieldnames(S2);
        for iField = 1:numel(fieldNameCell)
            fieldName = fieldNameCell{iField};
            if overwrite
                StructMerged.(fieldName) = S2.(fieldName);
            else
                if ~isfield(StructMerged, fieldName)
                    StructMerged.(fieldName) = S2.(fieldName);
                end
            end
        end
    end
end
