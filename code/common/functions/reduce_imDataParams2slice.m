function imDataParams_1Slice = reduce_imDataParams2slice(imDataParams, iSlice)
    
    imDataParams_1Slice = imDataParams;
    imDataParams_1Slice.images = imDataParams.images(:, :, iSlice, :, :);

end
