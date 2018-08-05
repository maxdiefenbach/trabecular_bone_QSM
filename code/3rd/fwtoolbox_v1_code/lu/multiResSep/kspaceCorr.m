load DAT;
TE = [4.1 7.3 10.5];

[ht, wd, nechoes] = size(dat);

%% before k-space misalignment correction
r1 = ift(dat(:,:,1)); r1 = r1./max(abs(r1(:)));
r2 = ift(dat(:,:,2)); r2 = r2./max(abs(r2(:)));
r3 = ift(dat(:,:,3)); r3 = r3./max(abs(r3(:)));
R = cat(2, r1(:), r2(:), r3(:));
cf = corrcoef(R);

%% after k-space misalignment correction
alpha = -.035; lpMat = get_linPhaseMat(alpha, ht, wd);
dat2 = dat;
dat2(:,:,2) = dat2(:,:,2).*lpMat;
c1 = ift(dat2(:,:,1)); c1 = c1./max(abs(c1(:)));
c2 = ift(dat2(:,:,2)); c2 = c2./max(abs(c2(:)));
c3 = ift(dat2(:,:,3)); c3 = c3./max(abs(c3(:)));
RC = cat(2, c1(:), c2(:), c3(:));
cfc = corrcoef(RC);