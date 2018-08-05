clear;
TE = [4.3680    7.5440    5.9520];
workPath = 'E:/dicom/15T5';
imgNameFmt = 'IM';

% indices to the starting and end slices
startInd = 30;
endInd = 60;

% process type: 1 single coil (I, Q)
%               2 single coil (M, I, Q)
%               3 multi-coil  (I, Q);
%               4 multi-coil  (M, I, Q)
processType = 2; 
if processType>2,
    ncoils = 4;
end

%-------------------------------------------------------------------------------
% Data processing
%-------------------------------------------------------------------------------
for k = startInd:endInd,
    
    disp(['Process image ' num2str(k)]);
    
    switch processType,
    case 1, 
        im1 = []; im2 = []; im3 = [];
	for l=(k-1)*6+1:k*6,
	    imgName = sprintf('%s%s%s%.4d.dcm', workPath, '/data/', imgNameFmt, l);
	if mod(l, 6) == 1
            im1 = double(dicomread(imgName));
        elseif mod(l, 6)==2
            im1 = im1 + i*double(dicomread(imgName));
        elseif mod(l, 6)==3
            im2 = double(dicomread(imgName));
        elseif mod(l, 6)==4
            im2 = im2 + i*double(dicomread(imgName));
        elseif mod(l, 6)==5
            im3 = double(dicomread(imgName));
        elseif mod(l, 6)==0
            im3 = im3 + i*double(dicomread(imgName));
        end
        end  
        
	sigSet = cat(3, conj(-i.*im1), conj(-i.*im2), conj(-i.*im3));
	[w, f, psi] = multiResSep(TE, sigSet, 0);
       
    case 2,	
        im1 = []; im2 = []; im3 = [];
	for l=(k-1)*9+1:k*9,
	    imgName = sprintf('%s%s%s%.4d.dcm', workPath, '/data/', imgNameFmt, l);
	if mod(l, 9) == 2
            im1 = double(dicomread(imgName));
        elseif mod(l, 9)==3
            im1 = im1 + i*double(dicomread(imgName));
        elseif mod(l, 9)==5
            im2 = double(dicomread(imgName));
        elseif mod(l, 9)==6
            im2 = im2 + i*double(dicomread(imgName));
        elseif mod(l, 9)==8
            im3 = double(dicomread(imgName));
        elseif mod(l, 9)==0
            im3 = im3 + i*double(dicomread(imgName));
        end
        end  
        
	sigSet = cat(3, conj(-i.*im1), conj(-i.*im2), conj(-i.*im3));
	[w, f, psi] = multiResSep(TE, sigSet, 0);
    
    case 3,	
	im1 = []; im2 = []; im3 = [];
	for l=(k-1)*(3*2*(ncoils+1))+1:k*(3*2*(ncoils+1)), % 3: echoes; 2: I, Q; ncoil+1: sum-of-squares
       
	    imgName = sprintf('%s%s%s%.4d.dcm', workPath, '/data/', imgNameFmt, l);
        
	    if ceil((l-(k-1)*(3*2*(ncoils+1)))/(2*(ncoils+1)))<=1,
            ctr = l-(k-1)*(3*2*(ncoils+1));
                if mod(ctr, 2)==1,
                    tmp = double(dicomread(imgName));
                else
                    tmp = tmp+i.*double(dicomread(imgName));
                    im1 = cat(3, im1, tmp);
                    tmp = [];
                end
        elseif ceil((l-(k-1)*(3*2*(ncoils+1)))/(2*(ncoils+1)))<=2,
            ctr = l-k*(3*2*(ncoils+1));
                if mod(ctr, 2)==1,
                    tmp = double(dicomread(imgName));
                else
                    tmp = tmp+i.*double(dicomread(imgName));
                    im2 = cat(3, im2, tmp);
                    tmp = [];
                end
        elseif ceil((l-(k-1)*(3*2*(ncoils+1)))/(2*(ncoils+1)))<=3,
            ctr = l-2*k*(3*2*(ncoils+1));
                if mod(ctr, 2)==1,
                    tmp = double(dicomread(imgName));
                else
                    tmp = tmp+i.*double(dicomread(imgName));
                    im3 = cat(3, im3, tmp);
                    tmp = [];
                end
        end                     
    end
    % discard sum-of-squares images
    im1 = im1(:,:,1:end-1);
    im2 = im2(:,:,1:end-1);
    im3 = im3(:,:,1:end-1);
    sigSet = cat(4, conj(-i.*im1), conj(-i.*im2), conj(-i.*im3));
    sigSet = permute(sigSet, [1 2 4 3]);
    [w, f, psi] = multiCoilSep(TE, sigSet, 0, 0);
    
    case 4,
    im1 = []; im2 = []; im3 = [];
    % 3: echoes; 3: M, I, Q; ncoil+1: sum-of-squares
	for l=(k-1)*(3*3*(ncoils+1))+1:k*(3*3*(ncoils+1)), 
       
	    imgName = sprintf('%s%s%s%.4d.dcm', workPath, '/data/', imgNameFmt, l);
        
	    if ceil((l-(k-1)*(3*3*(ncoils+1)))/(3*(ncoils+1)))<=1,
            ctr = l-(k-1)*(3*3*(ncoils+1));
                if mod(ctr, 3)==2,
                    tmp = double(dicomread(imgName));
                elseif mod(ctr, 3)==0
                    tmp = tmp+i.*double(dicomread(imgName));
                    im1 = cat(3, im1, tmp);
                    tmp = [];
                end
        elseif ceil((l-(k-1)*(3*3*(ncoils+1)))/(3*(ncoils+1)))<=2,
            ctr = l-k*(3*3*(ncoils+1));
                if mod(ctr, 3)==2,
                    tmp = double(dicomread(imgName));
                elseif mod(ctr, 3)==0
                    tmp = tmp+i.*double(dicomread(imgName));
                    im2 = cat(3, im2, tmp);
                    tmp = [];
                end
        elseif ceil((l-(k-1)*(3*3*(ncoils+1)))/(3*(ncoils+1)))<=3,
            ctr = l-2*k*(3*3*(ncoils+1));
                if mod(ctr, 3)==2,
                    tmp = double(dicomread(imgName));
                elseif mod(ctr, 3)==0,
                    tmp = tmp+i.*double(dicomread(imgName));
                    im3 = cat(3, im3, tmp);
                    tmp = [];
                end
        end                     
    end
    % discard sum-of-squares images
    im1 = im1(:,:,1:end-1);
    im2 = im2(:,:,1:end-1);
    im3 = im3(:,:,1:end-1);
    sigSet = cat(4, im1, im2, im3);
    sigSet = permute(sigSet, [1 2 4 3]);
    % combine into single-coil data
    [sc, coilSense] = srcCoilCombineGEMR1(sigSet, TE/1000, -440, 20, ones(ncoils, 1));
%     [w, f, psi] = multiResSep(TE, sc, 0);
    
%     [w, f, psi] = multiCoilSep(TE, sigSet, 0, 0);
            
    end

%     w = w./max(w(:))*1024;
%     f = f./max(f(:))*1024;
%         imgWatName = sprintf('%s%s%swater-%.4d', workPath, '/output/water/', imgNameFmt, k);
%         imgFatName = sprintf('%s%s%sfat-%.4d', workPath, '/output/fat/', imgNameFmt, k);
%     fid = fopen(imgWatName, 'wb'); fwrite(fid, w); fclose(fid);
%     fid = fopen(imgFatName, 'wb'); fwrite(fid, f); fclose(fid);
% 
% 	w = abs(w); w = w./max(w(:))*256;
%         f = abs(f); f = f./max(f(:))*256;
%         psi = psi-min(psi(:)); psi = psi./max(psi(:))*256;
%         
%         imgWatName = sprintf('%s%s%s%.4d-water.bmp', workPath, '/output/water/', imgNameFmt, k);
%         imgFatName = sprintf('%s%s%s%.4d-fat.bmp', workPath, '/output/fat/', imgNameFmt, k);
%         imgFmapName = sprintf('%s%s%s%.4d-fmap.bmp', workPath, '/output/fmap/', imgNameFmt, k);
%         
%         imwrite(uint8(w), imgWatName);
%         imwrite(uint8(f), imgFatName);
%         imwrite(uint8(psi), imgFmapName);
    eval(['sc' num2str(k) ' =sc;'])

end
