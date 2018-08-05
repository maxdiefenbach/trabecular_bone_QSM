function h = plot_WFIparamsSlice(WFIparams, iSl)
    
    if ~isfield(WFIparams, 'voxelSize_mm')
        WFIparams.voxelSize_mm = [1, 1, 1];
    end
    if ~isfield(WFIparams, 'fatFraction_percent')
        WFIparams = get_fatFraction_percent(WFIparams);
    end
    if isfield(WFIparams, 'R2s_Hz')
        WFIparams = get_T2s_ms(WFIparams);
    else
        R2s_Hz = 0;
        T2s_ms = 0;
    end
    if ~isfield(WFIparams, 'fieldmap_Hz')
        fieldmap_Hz = 0
    end

    unpack_struct2vars(WFIparams);
    
    h = figure;
    if isfield(WFIparams, 'residual') & isfield(WFIparams, 'iterations')
            subplot(4, 2, 1)
            imagesc(remove_outliers(abs(water(:, :, iSl)), [0, 99.9]))
            title('water')
            axis off
            subplot(4, 2, 2)
            imagesc(remove_outliers(abs(fat(:, :, iSl)), [0, 99.9]))
            title('fat')
            axis off
            subplot(4, 2, 3)
            imagesc(fatFraction_percent(:, :, iSl))
            title('PDFF [%]')
            axis off
            caxis([-5, 105])
            subplot(4, 2, 4)
            imagesc(fieldmap_Hz(:, :, iSl))
            title('fieldmap [Hz]')
            axis off
            subplot(4, 2, 5)
            caxis([-500, 500])
            imagesc(R2s_Hz(:, :, iSl))
            title('R2s [Hz]')
            axis off
            caxis([0, 300])
            subplot(4, 2, 6)
            imagesc(T2s_ms(:, :, iSl))
            title('T2s [ms]')
            axis off
            caxis([0, 100])
            subplot(4, 2, 7)
            imagesc(residual(:, :, iSl))
            title('residual')
            axis off
            subplot(4, 2, 8)
            imagesc(iterations(:, :, iSl))
            title('iterations')
            axis off
    else
            if isfield(WFIparams, 'R2s_Hz') & isfield(WFIparams, 'T2s_ms')
                subplot(3, 2, 1)
                imagesc(abs(water(:, :, iSl)))
                title('water')
                axis off
                subplot(3, 2, 2)
                imagesc(abs(fat(:, :, iSl)))
                title('fat')
                axis off
                subplot(3, 2, 3)
                imagesc(fatFraction_percent(:, :, iSl))
                title('PDFF [%]')
                axis off
                caxis([-5, 105])
                subplot(3, 2, 4)
                imagesc(fieldmap_Hz(:, :, iSl))
                title('fieldmap [Hz]')
                axis off
                subplot(3, 2, 5)
                caxis([-500, 500])
                imagesc(R2s_Hz(:, :, iSl))
                title('R2s [Hz]')
                axis off
                caxis([0, 300])
                subplot(3, 2, 6)
                imagesc(T2s_ms(:, :, iSl))
                title('T2s [ms]')
                axis off
                caxis([0, 100])
            else
                subplot(2, 2, 1)
                imagesc(abs(water(:, :, iSl)))
                title('water')
                axis off
                subplot(2, 2, 2)
                imagesc(abs(fat(:, :, iSl)))
                title('fat')
                axis off
                subplot(2, 2, 3)
                imagesc(fatFraction_percent(:, :, iSl))
                title('PDFF [%]')
                axis off
                caxis([-5, 105])
                subplot(2, 2, 4)
                imagesc(fieldmap_Hz(:, :, iSl))
                title('fieldmap [Hz]')
                axis off
            end
            set(h, 'Name', WFIparams.method);

    end
