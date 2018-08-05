function [x,y] = cofmass (data)

%   Function finds the Center of Mass (moment 1) of an image
%   and returns the coordinates [x,y].
%   [x,y] = com(data);
%   data: 2D image

[n1,n2]=size(data);

M0 = sum(sum(abs(data))); % zero moment

temp_sumx = 0; temp_sumy = 0;
for r = 1:n1
    for c = 1:n2
        temp_sumx = temp_sumx + abs(data(r,c)) * r;
        temp_sumy = temp_sumy + abs(data(r,c)) * c;
    end
end

x = floor(temp_sumx/M0);
y = floor(temp_sumy/M0);