function d = euclidean(p1,p2)
% Finds the Euclidean distance between 2 pts.
% d = euclidean(p1,p2)
% p1 and p2 are in the form of [x1,y1] and [x2,y2]
d = sqrt ( power(abs(p1(1)-p2(1)),2) + power(abs(p1(2)-p2(2)),2) );