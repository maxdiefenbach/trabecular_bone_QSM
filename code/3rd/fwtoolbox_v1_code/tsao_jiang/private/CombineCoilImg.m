% Phase correct based on maximum TE and combine coil images
%
% Jeffrey Tsao, 2011
function comboimg = CombineCoilImg(img) % img is in x,y,z,coil,TEs  
  % Find maximum echo
  [maxval,maxTE] = max(sum(sum(sum(sum(abs(img),1),2),3),4));
  clear maxval;

  % Calculate low-pass phase version of maximum echo image;
  CorrectImg = ifftn(img(:,:,:,:,maxTE));
  maxkspace = sum(abs(CorrectImg).^2,4);
  [maxval,maxidx] = max(maxkspace(:));
  [maxkx,maxky,maxkz] = ind2sub(size(maxkspace),maxidx); clear maxidx maxval maxkspace;
  filterwidth = max(bitshift([size(img,1),size(img,2),size(img,3)],-3),8);
  [kx,ky,kz] = ndgrid( mod(((1:size(img,1))-maxkx) + bitshift(size(img,1),-1), size(img,1))-bitshift(size(img,1),-1), ...
                       mod(((1:size(img,2))-maxky) + bitshift(size(img,2),-1), size(img,2))-bitshift(size(img,2),-1), ...
                       mod(((1:size(img,3))-maxkz) + bitshift(size(img,3),-1), size(img,3))-bitshift(size(img,3),-1) );
  tmpfilt = exp(-((kx./filterwidth(1)).^2 + (ky./filterwidth(2)).^2 + (kz./filterwidth(3)).^2));
  clear kx ky kz maxkx maxky maxkz filterwidth;
  CorrectImg = CorrectImg.*repmat(tmpfilt,[1,1,1,size(CorrectImg,4)]); clear tmpfilt;
  CorrectImg = exp(-1i*angle(fftn(CorrectImg)));
  clear maxTE;
  
  % Correct phase based on maximum echo
  img = img.*repmat(CorrectImg,[1,1,1,1,size(img,5)]); clear CorrectImg;
  
  if size(img,4)==1, 
    comboimg = img;
  else
    % Optimal combine
    img = img.*abs(img);   % Magnitude-squared weighted
    comboimg = exp(1i*angle(sum(img,4))) .* sqrt(mean(abs(img),4));  % Use RMS magnitude and SoS-combined phase
  end
end
