function unpack_struct2vars(Struct)
% unpack_struct2vars(Struct)
% 
% unpacks values from struct fields 
% to workspace variables

    fieldNameCell = fieldnames(Struct);
    for i = 1:numel(fieldNameCell)
        fieldName = fieldNameCell{i};
        assignin('caller', fieldName, Struct.(fieldName));
    end

end