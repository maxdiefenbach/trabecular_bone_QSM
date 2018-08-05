function psiHat = wlsFit(psiHat, S)
mask = imMagMask(S);
[ht, wd] = size(mask);
magMask = mean(abs(S), 3).*mask;
% prepare basis images
xBasisImg = repmat((-wd/2:(wd/2-1)), [ht, 1]);
yBasisImg = repmat((-ht/2:(ht/2-1))', [1, wd]);
xyBasis   = yBasisImg*diag((-wd/2:(wd/2-1)));
xSquareBasis = xBasisImg.*xBasisImg;
ySquareBasis = yBasisImg.*yBasisImg;
dcBasis   = ones(ht, wd);
xb = xBasisImg(mask); xb = xb(:);
yb = yBasisImg(mask); yb = yb(:);
dc = dcBasis(mask); dc = dc(:);
xs = xSquareBasis(mask); xs = xs(:);
ys = ySquareBasis(mask); ys = ys(:);
xy = xyBasis(mask); xy = xy(:);
H = [dc xb yb xs ys xy];
wv = magMask(mask); wv = wv(:);
pdv = psiHat(mask); pdv = pdv(:);
HpWH = H'*((wv*ones(1,size(H, 2))).*H);
aw = HpWH\(H'*(wv.*pdv));
% figure; imshow(dcBasis*aw(1)+xBasisImg*aw(2)+yBasisImg*aw(3), []);
psiHat = dcBasis*aw(1)+xBasisImg*aw(2)+yBasisImg*aw(3)+xSquareBasis*aw(4)+...
    ySquareBasis*aw(5)+xyBasis*aw(6);