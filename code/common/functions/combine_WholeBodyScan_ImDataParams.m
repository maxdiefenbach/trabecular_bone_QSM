function ImDataParams = combine_WholeBodyScan_ImDataParams(matFileCell, dimension)
    
    filename = matFileCell{1};
    
    SeriesDate = filename(1:8);
    SeriesTime = {filename(10:15)};
    SeriesNumber = {filename(17:20)};
    
    tmp = load(filename);
    ImDataParams = tmp.ImDataParams;
    signal = ImDataParams.signal;
    for i = 2:numel(matFileCell)
        filename = char(matFileCell(i));
        
        SeriesNumber{end+1} = filename(17:20);
        SeriesTime{end+1} = filename(10:15);
        
        tmp = load(filename);
        signal = cat(dimension, signal, tmp.ImDataParams.signal);
    end
    
    ImDataParams.SeriesDate = SeriesDate;
    ImDataParams.SeriesTime = SeriesTime;
    ImDataParams.SeriesNumber = SeriesNumber;
    ImDataParams.signal = signal;

end