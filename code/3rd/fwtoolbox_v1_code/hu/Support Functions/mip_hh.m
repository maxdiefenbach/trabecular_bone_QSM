function output = mip_hh(input,dim)
%function mip_hh.m
%Simple projection calculation of maximum-intensity projection of 3D data.
%Syntax: output = mip_hh(input,dim);
%   Input: 3D data matrix
%   Dim: dimension 1-x, 2-y, 3-z

[xres yres zres]=size(input);

if dim == 1
    output = zeros(yres,zres);
    for ii = 1:yres
        for jj = 1:zres
            temp = input(:,ii,jj);
            output(ii,jj) = max(temp);
        end
    end
end

if dim == 2
    output = zeros(xres,zres);
    for ii = 1:xres
        for jj = 1:zres
            temp = input(ii,:,jj);
            output(ii,jj) = max(temp);
        end
    end
end

if dim == 3
    output = zeros(xres,yres);
    for ii = 1:xres
        for jj = 1:yres
            temp = input(ii,jj,:);
            output(ii,jj) = max(temp);
        end
    end
end