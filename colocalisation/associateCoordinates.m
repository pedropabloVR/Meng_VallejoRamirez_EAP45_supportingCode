% This function takes in two sets of points and matches points that are
% closer together than some specified search radius.
% 
% The functions loops over all points in set 1 and checks whether there are
% points in set 2 that are closer than the search radius removed from the
% current point from set 1. If so, the point in set 2 that is the closest
% to the current point from set 1 is considered to be 'associated' with
% that point. The ordered coordinates of the matched points are returned.
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
%   X1,Y1..........: ordered coordinates of points from the first set of
%                    points that were matched with a point from second set
%   X2_min,Y2_min..: ordered coordinates of points from the second set of 
%                    points that were matched with a point from first set
% 
% Pedro Vallejo Ramirez, Laser Analytics Group
% Modified by Ezra Bruggeman on 09/09/2018 to also return number of
% neighbours within R_search + a bug fix.
% 
% Last updated on 6 Oct 2018


function [X1,Y1,X2_min,Y2_min,N_local] = ...
    associateCoordinates(X1,Y1,X2,Y2,R_search)

% Get the number of localizations in each channel
N_loc1 = size(X1,1);
N_loc2 = size(X2,1);

% Initialize vectors
N_local = zeros(N_loc1,1);
R_min   = zeros(N_loc1,1);
X2_min  = zeros(N_loc1,1);
Y2_min  = zeros(N_loc1,1);

% Loop over all localizations (x1,y1)
for i = 1:N_loc1

    % Calculate distance between every combination of (x1,y1) and (x2,y2)
    R = zeros(N_loc2,1);
    for j = 1:N_loc2
        R(j) = sqrt((X1(i)-X2(j))^2 + (Y1(i)-Y2(j))^2);
    end
    
    % get indices of R that are within R_search (that is, the indices of
    % the particles in ch2 which are found in close proximity to ch1).
    
    
    % Count number of localizations in (X2,Y2) with r < R_search to the i-th localization in (X1,Y1)
    N_local(i) = sum(R < R_search);
    
    % Get distance between the i-th localization in (X1,Y1) to the closest localization in (X2,Y2)
    R_min(i) = min(R);
    
    % Get index of this value
    idx = find(R == min(R));
    idx = idx(1); % take the first one if there are multiple

    % Get the coordinates of the nearest localization in (X2,Y2)
    X2_min(i) = X2(idx);
    Y2_min(i) = Y2(idx);
end

% Filter out all localizations (x1,y1) that were not matched with another
% localization (x2,y2)
X1(N_local == 0) = [];
 Y1(N_local == 0) = [];
X2_min(N_local == 0) = [];
Y2_min(N_local == 0) = [];
% X1(N_local ~= 1) = [];
% Y1(N_local ~= 1) = [];
% X2_min(N_local ~= 1) = [];
% Y2_min(N_local ~= 1) = [];
% currently this removes values of N_local that are zero, but it also
% removes values larger than 1 (particles with more than 1 neighbor at
% R_search). Changing the code to account for this 26/05/2019

% % Display
% figure('Color','white','name','Histogram of R_min','Units',...
%     'normalized','OuterPosition',[0.2 0.2 0.6 0.5]);
% 
% subplot(1,2,1)
% plot(X1,Y1,'+')
% hold on
% plot(X2_min,Y2_min,'r+')
% axis equal
% legend('Red channel','Green channel');
% 
% subplot(1,2,2)
% hist(R_min,0:10:500)
% xlim([0 500])
% xlabel('R_{offset} (nm)')
% title('Histogram of chromatic offset')
end
