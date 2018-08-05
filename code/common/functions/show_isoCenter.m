function show_isoCenter(ImDataParams)

    isoCenter_REC = round(ImDataParams.isoCenter_REC);
    magnitude = abs(ImDataParams.signal);
    center = floor((size(magnitude) + 1) ./ 2);
    center = center(1:3);

    plot_slice(magnitude, isoCenter_REC(3));
    hold on;
    plot(isoCenter_REC(2), isoCenter_REC(1), 'x', 'Color', 'red', 'LineWidth', 2)
    if center(3) == isoCenter_REC(3)
        plot(center(2), center(1), 'x', 'Color', 'blue', 'LineWidth', 2)
    else
        fprintf('Image center coordinates: [%d, %d, %d]\n', center)
    end
    hold off;

end