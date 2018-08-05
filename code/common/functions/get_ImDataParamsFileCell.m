function ImDataParamsFileCell = get_ImDataParamsFileCell(filepath)
    FileStruct = dir(fullfile(filepath, '*_ImDataParams.mat'));
    filenameCell = {FileStruct.name};
    for i = 1:numel(filenameCell)
        filenameCell{i} = fullfile(filepath, filenameCell{i});
    end
    ImDataParamsFileCell = filenameCell;
end