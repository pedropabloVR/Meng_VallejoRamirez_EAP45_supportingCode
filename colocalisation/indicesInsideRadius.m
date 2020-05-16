% This function takes in two sets of points and finds the indices of points
% in the second set within a given distance from the first set.
% 
% INPUT:
%   X1,Y1..........: nx1 doubles that contain xy-coordinates of the first
%                    set of points
%   X2,Y2..........: mx1 doubles that contain xy-coordinates of the second
%                    set of points
%   R_search.......: maximum distance that two points can be removed from
%                    each other to be considered associated
% 
% OUTPUT:
%   idx2 ..............: indices of all the x2,y2 inside search radius
%   x2_close,y2_close..: ordered coordinates of points from the second set of 
%                    points within R_search from a point from first set
% 
% Pedro Vallejo Ramirez, Laser Analytics Group

% Last updated on 18 Jan 2020


function [idx2,x2_close,y2_close] = ...
    indicesInsideRadius(X1,Y1,X2,Y2,R_search)

% Get the number of localizations in each channel
N_loc1 = size(X1,1);
N_loc2 = size(X2,1);


% Calculate distance between every combination of (x1,y1) and (x2,y2)
R = zeros(N_loc2,1);
for j = 1:N_loc2
    R(j) = sqrt((X1-X2(j))^2 + (Y1-Y2(j))^2);
end

% find elements in (x2,y2) which are inside the search radius

idx2 = (R < R_search ) & (R>0); % indices of all the x2,y2 inside search radius, except x1,y1 (itself)
x2_close = X2(idx2);
y2_close = Y2(idx2);


end