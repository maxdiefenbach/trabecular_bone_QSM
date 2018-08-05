function [wat, fat, psi, resd] = getSepResult(rawFname, snumSet, TE, alpha, biCorrFlag)
workPath ='.';
dat = recon3d(rawFname);

for ctr=1:length(snumSet)
    
snum = snumSet(ctr);
disp(['Process frame ' num2str(snum)]);
im1 = dat(:,:,1,snum); im2 = dat(:,:,2,snum); im3 = dat(:,:,3,snum);
[ht, wd] = size(im1);
lpMat = get_linPhaseMat(alpha, ht, wd);
im2 = im2.*lpMat;

if (exist('biCorrFlag'))
    if biCorrFlag
        im2c = im2;
        mask = imMagMask(dat(:,:,:,snum));
        m2 = mask;
        for k=2:ht,
            for l=1:wd,
                if (abs(im2(k, l))<.8*abs(im3(k, l)) || abs(im2(k, l))>1.2*abs(im1(k, l)))
                    im2c(k, l) = (abs(im1(k, l))*abs(im3(k, l)))^.5.*exp(i*angle(im2(k-1, l)));
                    m2(k, l) = 0;
                end
            end
        end
        im2 = im2c;
        mask = cat(3, mask, m2, m2); mask = mask.*1; figure; imshow(mask); drawnow;
    end
end

[wat, fat, psi, resd] = multiResSep(TE, cat(3, im1, im2, im3), 0);
wat = abs(wat); wat = wat./max(wat(:))*256;
fat = abs(fat); fat = fat./max(fat(:))*256;
psi = psi-min(psi(:)); psi = psi./max(psi(:))*256;
imgWatName = sprintf('%s%s%s%.4d-water.bmp', workPath, '/output/water/', 'water', snum);
imgFatName = sprintf('%s%s%s%.4d-fat.bmp', workPath, '/output/fat/', 'fat', snum);
imgFmapName = sprintf('%s%s%s%.4d-fmap.bmp', workPath, '/output/fmap/', 'fmap', snum);
imwrite(uint8(wat), imgWatName);
imwrite(uint8(fat), imgFatName);
imwrite(uint8(psi), imgFmapName);
        
end