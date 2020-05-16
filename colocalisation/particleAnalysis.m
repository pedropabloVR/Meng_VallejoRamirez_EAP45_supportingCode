% particleAnalysis.m
% 
% This function performs a particle analysis on blobs in a binary mask.
% 
% INPUT:
%   mask..........: binary mask
%   magnification.: magnification factor if generateImage.m was used while
%                   making the mask (otherwise magnification = 1)
%   locs..........: 
% 
% OUTPUT:
%   measurements..: table containing different measurements for individual
%                   cluster
%   mask_labeled..: the mask, with every cluster labeled with different
%                   integers (cluster one will be filled with 1, the second
%                   cluster with 2, ...)
%   mask_colored..: pseudo-coloured version of mask_labeled (for
%                   visualization purposes)
% 
% Author: Ezra Bruggeman, Laser Analytics Group
% 
% Last updated on 6 Sept 2018

% Modified by PVR to be applied for widefield images - no need to feed in
% localisations and it no longer returns a density measurement. 

function [measurements, mask_labeled, mask_colored] = ...
    particleAnalysis(mask, magnification)

% Identify individual clusters and label them with different integer values
mask_labeled = bwlabel(mask);

% Get a pseudo-coloured RGB version of the labeled mask (only for display
% purposes)
mask_colored = label2rgb(mask_labeled, 'hsv', 'k', 'shuffle');

% Get cluster properties using the built-in MATLAB function regionprops
clusterMeasurements = regionprops(mask_labeled, mask, 'all');
numberOfClusters = size(clusterMeasurements, 1);

% Get coordinates of centroids of clusters
allClusterCentroids = [clusterMeasurements.Centroid];
x_centroid = allClusterCentroids(1:2:end);
y_centroid = allClusterCentroids(2:2:end);

% Get area of clusters
allClusterArea = [clusterMeasurements.Area];
allClusterArea = allClusterArea*(magnification^2);

% Get equivalent diameter of clusters
allClusterEquivDiameter = [clusterMeasurements.EquivDiameter];
allClusterEquivDiameter = allClusterEquivDiameter*magnification;

% Get eccentricity of clusters
allClusterEccentricity = [clusterMeasurements.Eccentricity];

% % Get number of localizations per cluster
% %x = locs.x/magnification; y = locs.y/magnification;
% x = locs(:,1)/magnification; y = locs(:,2)/magnification;
% NumberOfLocsInCluster = zeros(numberOfClusters,1);
% for i=1:numberOfClusters
%     % Get mask for current cluster only
%     maskCurrentCluster = mask_labeled;
%     maskCurrentCluster(maskCurrentCluster~=i) = 0;
%     maskCurrentCluster(maskCurrentCluster==i) = 1;
%     maskCurrentCluster = logical(maskCurrentCluster);
%     % Get number of localizations inside mask
%     counter = 0;
%     for j=1:size(x,1)
%         if maskCurrentCluster(round(x(j)),round(y(j))) == 1
%             counter = counter + 1;
%         end
%     end
%     NumberOfLocsInCluster(i) = counter;
% end
% 
% % Get density of localizations per cluster (localizations/area)
% DensityLocsInCluster = NumberOfLocsInCluster'./allClusterArea;
% 
if numberOfClusters > 1
    % Get nearest neighbour distance
    NearestNeighbourDist = zeros(numberOfClusters,1);
    for i=1:numberOfClusters
        xi = x_centroid(i);
        yi = y_centroid(i);
        ptCloud = [x_centroid' y_centroid'];
        ptCloud(i,:) = [];
        k = dsearchn(ptCloud,[xi yi]);
        x_closest = x_centroid(k);
        y_closest = y_centroid(k);
        dist = sqrt((xi - x_closest)^2 + (yi - y_closest)^2);
        NearestNeighbourDist(i) = dist;
    end
    NearestNeighbourDist = NearestNeighbourDist*magnification;
else
    NearestNeighbourDist = NaN;
end
% Convert centroid coordinates to right scale
x_centroid = x_centroid*magnification;
y_centroid = y_centroid*magnification;

% Save measurements in a table
id = 1:numberOfClusters;
measurements = [id' x_centroid' y_centroid' allClusterArea' allClusterEquivDiameter' allClusterEccentricity' NearestNeighbourDist];
measurements = array2table(measurements,'VariableNames',{'ClusterID','XCentroid','YCentroid','Area','EquivDiameter','Eccentricity','NearestNeighbourDist'});

end