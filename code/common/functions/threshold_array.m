function arrayThresholded = threshold_array(array, threshold_percent)

    maxVal = max(array(:));
    arrayThresholded = array >= threshold_percent/100 * maxVal .* ones(size(array));

end
