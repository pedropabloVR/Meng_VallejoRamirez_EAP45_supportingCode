% wf_colocalisation_fromMask.m imports two binary masks with puncta-like 
% particles segmented from two-colour microscopy images. The script determines 
% the number of neighbouring particles from mask 2 within a given search
% search radius from mask 1, and viceversa. 


% The script:
% - Imports binary masks for two channels of an HIV-ESCRTII acquisition.
% - Classifies the thresholded regions using region props
% - Removes large clusters using an area threshold
% - Finds the number of clusters associated between the two channels
%   within a given search radius. 

% Depends on functions:
% - particleAnalysis.m
% - associateCoordinates.m
% - indicesInsideRadius.m
% - natsortfiles.m
% - natsort.m

% Output:
% - .csv file with the results from the association analysis
% - overlaid images of the segmented binary mask used to compute the
% centroids of the particles and the raw image data, to check for the
% faithfulness of the segmentation compared to the raw data. 

% Pedro Vallejo Ramirez
% Created: 12/05/2019
%   natsort.m and natsortfiles.m kindly provided by Stephen Cobeldick through
%   Matlab file exchange. Stephen Cobeldick (2020). Natural-Order Filename Sort (https://www.mathworks.com/matlabcentral/fileexchange/47434-natural-order-filename-sort), MATLAB Central File Exchange. Retrieved May 15, 2020. 

% Updates:

%       13/01/2020 to extract the cluster sizes of all particles (EAP45
%       and Gag), as well as the specific cluster sizes of the EAP45 particles
%       associated with Gag and those not associated with Gag for comparison.

%       18/01/2020 added lines to extract sizes of EAP45 particles in the
%       neighborhood of the associated EAP45 particles for size comparison. 

 
clear all
close all

% Load image data

% ch1 is green (Gag) and ch2 is magenta (EAP45)

%% HeLa Cells 

% Sample directories for [experiment repeat]
% path_ch1    = 'D:\Experiments\hiv escrt\2019_04_02_Pedro_Bo_HIVESCRT_exp6\slide1well3\AvgInt\weka'; % Reference channel - 647
% path_ch2    = 'D:\Experiments\hiv escrt\2019_04_02_Pedro_Bo_HIVESCRT_exp6\slide1well3\AvgInt\warped\corrected\weka'; % registered channel

%% Set parameters

ch1         = dir([path_ch1,filesep,'*','_weka.tif']); % files imported end in _weka.tif
ch2         = dir([path_ch2,filesep,'*','_weka.tif']);

ch1_token     = '_647';
ch2_token     = '_488';

% p_ch1         = 5; % minimum number of pixels for a blob for background reduction.
% p_ch2         = 5; 
save_flag     = 1; % check here to save all results
crop_flag     = 0; % check here to crop the sides of the image to remove background
show          = 0; % check here to show intermediate results 
NumBins       = 100;   % number of bins
max_dist      = 3; % 10 pixels is equivalent to ~1 um for these diffraction-limited images. 
r_search      = 0:1:1; %r_search limits from half a pixel to 10 pixels away
area_max      = 100; % maximum area in pixels for a blob in regionprops
count         = 1;
save_qualityControl = 1;

% radii of disk structural elements for erosion and dilation
se_erode_radius  = 1;
se_dilate_radius = 10;
se_erode  = strel('disk',se_erode_radius);
se_dilate = strel('disk',se_dilate_radius);

% Make output folder
if save_flag
output_dir_path = fullfile(path_ch1,'result');
if exist(output_dir_path, 'dir')
    opts.Interpreter = 'tex';
    opts.Default = 'Continue';
    quest = '\fontsize{12}An output folder ''tracking\_results'' already exists. If you continue, data in this folder might be overwritten.';
    answer = questdlg(quest,'Message','Cancel','Continue',opts);
else
    mkdir(output_dir_path);
end
end 

% Extract filenames and sort the numbers in ascending order.
ch1         = extractfield(ch1,'name');
ch2         = extractfield(ch2,'name');
ch1_sort    = natsortfiles(ch1);
ch2_sort    = natsortfiles(ch2);

% initalize results array to hold results
%results = cell(size(ch1,2),6);
results = cell(size(ch1,2)*size(r_search,2),8);

% save parameters used for analysis
save([output_dir_path filesep 'parameters']);
size_results_ch1        = [];
size_results_ch2        = [];
ch1_min_sizes_assoc     = [];
ch1_min_sizes_nonAssoc  = [];
ch1_min_sizes_NonAssoc_neighbours = [];
count_empty = 0; % variable to keep count on when an array is first filled
count_empty_assoc = 0;
size_dist = 15;

%% Loop over images in folder 

for i = 1:size(ch1,2)
    
    %%---Import masks red and green channel---%%
    
    filename_ch1 = ch1_sort{i};
    filename_ch2 = ch2_sort{i};
    
    % Get file names
    area_id_split = strsplit(filename_ch1,'_');
    area_id = [area_id_split{1}];
    
    % print filename info
    disp('## ESCRTII Mask analysis ##')
    disp(' ') 
    disp(['488 channel: ' filename_ch1])
    disp(['647 channel: ' filename_ch2])
    
    % read in files
    ch1_bw = imread(fullfile(path_ch1,filename_ch1));
    ch2_bw = imread(fullfile(path_ch2,filename_ch2));
    
    % check here if image is already inverted or requires inversion
    mode_ch1 = mode(mode(ch1_bw));
    mode_ch2 = mode(mode(ch2_bw));
   
    if mode_ch1 && mode_ch2 == 255
        ch1_bw = imcomplement(ch1_bw);
        ch2_bw = imcomplement(ch2_bw);
    elseif mode_ch1 == 0 && mode_ch2 ~= 0
         ch2_bw = imcomplement(ch2_bw);
    elseif mode_ch2 == 0 && mode_ch1 ~= 0
         ch1_bw = imcomplement(ch1_bw);
    end

%     
    %% Perform particle analysis
    % use region props to extract regions within the cell area and discard
    % background

     magnification = 1;
    [ch1_measurements, ch1_mask_labeled, ch1_mask_colored] = particleAnalysis(ch1_bw, magnification);
    [ch2_measurements, ch2_mask_labeled, ch2_mask_colored] = particleAnalysis(ch2_bw, magnification);

%    % Filter out regions which are too large
%    ch1_measurements = ch1_measurements(ch1_measurements.Area      < area_max ,:);
%    ch2_measurements = ch2_measurements(ch2_measurements.Area      < area_max ,:);
%    
   % get mask after filtering out the large regions
   ch1_mask_filtered = ismember(ch1_mask_labeled,ch1_measurements.ClusterID); % ok so this magic worked! 
   ch2_mask_filtered = ismember(ch2_mask_labeled,ch2_measurements.ClusterID);  
     
   % show mask with only small regions chosen
   % Here the user can check whether the thresholded image corresponds well to the original image data.
   % In other words, if the segmented blobs correspond to bright spots in
   % the image.
%    if show 
%         figure('units','normalized','outerposition',[0 0 1 1]);
%         subplot(2,2,1), imshow(ch1_mask_filtered);
%         title(['Binarized Ch1 after removing large elements']);
%         subplot(2,2,2), imshow(ch1_im,[]);
%         title(['Ch1 raw: ' area_id_split{1} ' ' area_id_split{3}]);
%         subplot(2,2,3), imshow(ch2_mask_filtered);
%         title(['Binarized Ch2 Ch1 after removing large elements']);
%         subplot(2,2,4), imshow(ch2_im,[]);
%         title(['Ch2 raw: ' area_id_split{1} ' ' area_id_split{3}]);
%         pause % wait for user to click a button to continue script 
%    end 
   
    if save_qualityControl
         figure;imshowpair(ch1_bw,ch2_bw); % save image of the gfp and its nanobooster together 
         im = imfuse(ch1_bw,ch2_bw);
%        fig1 = figure;imshowpair(bwlabel(ch1_mask_filtered),ch1_im)
         imwrite((im),fullfile(output_dir_path,[area_id ch1_token 'weka_maskOverlay.tif']))

        %saveas(fig1,fullfile(output_dir_path,[area_id ch1_token 'weka_maskOverlay.tif']));
%        fig2 = figure;imshowpair(bwlabel(ch2_mask_filtered),ch2_im)
%        saveas(fig2,fullfile(output_dir_path,[area_id ch2_token '_maskOverlay.tif']));
%        
    end 
     
   % Check again if only acceptable regions are found
   % to-do 
   % - potentially erode the filtered mask, since the regionprops generally
   % chooses areas much larger than the original area of the blob.  
   
   % You can use the function vislabels to visualize the number
   % corresponding to each cluster in the labelled mask, to check the
   % corresponding features output by region props. The user can then use
   % this by visual inspection to filter out large regions (which are
   % likely many blobs together).     
   
   % Save imagepair of the raw data with the mask on top, for later
   % inspection 
   
   % code addition 2020-01-10
   % Save the sizes of the EAP45 clusters for each region analyzed
   % Bonus points if you can indicate which ones colocalise with Gag
   % repeat for Gag
   
   sizes_ch1 = ch1_measurements.Area;
   sizes_ch2 = ch2_measurements.Area;
   % create cell to store area results
   
   size_results_ch1 = vertcat(size_results_ch1,sizes_ch1);   
   size_results_ch2 = vertcat(size_results_ch2,sizes_ch2);   
   
   %size_results_ch2(i,1) = sizes_ch2;   
     
  % Get centroid coordinates from both channels
    X1 = ch1_measurements.XCentroid;
    Y1 = ch1_measurements.YCentroid;
    X2 = ch2_measurements.XCentroid;
    Y2 = ch2_measurements.YCentroid;
   
  % Loop over array of search distances
   
  for j = 1:size(r_search,2)
    
      max_dist = r_search(j);
   [X1_min,Y1_min,X2_min,Y2_min,N12_local] = associateCoordinates(X1,Y1,X2,Y2,max_dist);
   N12_num = sum(N12_local);
   if isempty(X1_min)
        disp(['! Data from ROI ' area_id ' were skipped because there were no clusters in channel 2 that were ']);
        disp(['less than max_dist = ' num2str(max_dist) ' removed from a cluster in channel 1.']);
        disp(' ');
        N12 = 0;
        %continue
    else
        dist12 = sqrt((X1_min-X2_min).^2 + (Y1_min-Y2_min).^2);
        dist12 = array2table(dist12,'VariableNames',{'dist'});
        
        % Get sampleID column
        clear var current_sampleID
        [current_sampleID{1:size(dist12,1)}] = deal(area_id);
        dist12.sampleID = current_sampleID';
        dist12 = dist12(:,[end 1:end-1]);

        N12 = size(dist12,1);
   end
    
    % Find neighbours within some radius, with second channel as reference channel
    [X2_min,Y2_min,X1_min,Y1_min,N21_local] = associateCoordinates(X2,Y2,X1,Y1,max_dist);
    N21_num = sum(N21_local);
    if isempty(X1_min)
        disp(['! Data from ROI ' area_id ' were skipped because there were no clusters in channel 1 that were ']);
        disp(['less than max_dist = ' num2str(max_dist) ' removed from a cluster in channel 2.']);
        disp(' ');
        N21 = 0;
        %continue
    else
        dist21 = sqrt((X1_min-X2_min).^2 + (Y1_min-Y2_min).^2);
        dist21 = array2table(dist21,'VariableNames',{'dist'});
        % Get sampleID column
        clear var current_sampleID
        [current_sampleID{1:size(dist21,1)}] = deal(area_id);
        dist21.sampleID = current_sampleID';
        dist21 = dist21(:,[end 1:end-1]);
   
        N21 = size(dist21,1);
               
        % get sizes of the elements from ch1 which are within R_search from
        % ch2 (matching X1_min and Y1 min to its index), i.e. the
        % associated elements. 
        [tf,ch1_idx] = ismember(X1_min,X1);
        ch1_min_sizes_assoc = vertcat(ch1_min_sizes_assoc,ch1_measurements.Area(ch1_idx)); %these are the sizes of the EAP45 blobs which were found in proximity of a Gag blob
        
        if ~isempty(ch1_min_sizes_assoc)
            count_empty_assoc = count_empty_assoc +1;
        end
        
        if i == 1
          sizes_labels = cell(size(ch1_min_sizes_assoc,1),1); 
          [sizes_labels{1:end}] = deal(area_id);
        elseif count_empty_assoc == 1
          sizes_labels = cell(size(ch1_min_sizes_assoc,1),1); 
          [sizes_labels{1:end}] = deal(area_id);
        else 
         [sizes_labels{end:end+size(ch1_measurements.Area(ch1_idx),1)}] = deal(area_id);   
        end 
          
        %% Find the neighbours of EAP45 particles within R pixels, and extract get their sizes
        % Iterate through the indices of the ch1 particles
        nonAssoc_size_change = size(ch1_min_sizes_NonAssoc_neighbours,1);
        for k = 1:size(ch1_idx,1)
            % get centroid of first colocalised particle
            x_k = ch1_measurements.XCentroid(ch1_idx(k));
            y_k = ch1_measurements.YCentroid(ch1_idx(k));
            
            % get indices of neighbouring EAP45 particles to an associated
            % particle
            [idx2,x2_close,y2_close] = indicesInsideRadius(x_k,y_k,X1,Y1,size_dist);
            % extract the sizes of these particles
            
            ch1_min_sizes_NonAssoc_neighbours = vertcat(ch1_min_sizes_NonAssoc_neighbours,ch1_measurements.Area(idx2));
            % set a counter for first instance in which variable above
            % isn't empty
 
        end
        
        if ~isempty(ch1_min_sizes_NonAssoc_neighbours)
            count_empty = count_empty +1;
        end
        
        if i == 1
         
          sizes_labels_nonAssoc = cell(size(ch1_min_sizes_NonAssoc_neighbours,1),1); 
          [sizes_labels_nonAssoc{1:end}] = deal(area_id);
        elseif count_empty == 1
           sizes_labels_nonAssoc = cell(size(ch1_min_sizes_NonAssoc_neighbours,1),1); 
          [sizes_labels_nonAssoc{1:end}] = deal(area_id);
        else 
         
         size_after = size(ch1_min_sizes_NonAssoc_neighbours,1);   
         size_change =  size_after - nonAssoc_size_change;
         if ~isempty(ch1_min_sizes_NonAssoc_neighbours)
            [sizes_labels_nonAssoc{end:end+size_change,1}] = deal(area_id);
         end 
        end 
         
        % retrieve all other indices of the non colocalising elements
        % get their sizes
        [C, ia, ib] = setxor(X1_min,X1);
        ch1_min_sizes_nonAssoc = vertcat(ch1_min_sizes_nonAssoc,ch1_measurements.Area(ib)); 
        
        % record for 1 pix distance and accumulate for all areas. 
    end
      
    % Get percentages of colocalizations
    percWithNeighbour1 = round((N12/size(N12_local,1))*100, 1);
    percWithNeighbour2 = round((N21/size(N21_local,1))*100, 1);
    str1 = sprintf('%d out of %d (%.1f%%) clusters in channel 1 have at least one neighbour in channel 2 within search radius %.2f',N12, size(N12_local,1), percWithNeighbour1, max_dist);
    str2 = sprintf('%d out of %d (%.1f%%) clusters in channel 2 have at least one neighbour in channel 1 within search radius %.2f',N21, size(N21_local,1), percWithNeighbour2, max_dist);
    disp(str1); disp(str2);
    
    
        % Write away results
        
        results{count,1} = area_id;
        results{count,2} = size(N12_local,1); % total number of red blobs
        results{count,3} = size(N21_local,1); % total number of green blobs
        results{count,4} = max_dist;
        results{count,5} = N12; % number of red blobs with at least 1 green blob within max_dist
        results{count,6} = N21; % number of green blobs with at least 1 red blob within max_dist
        results{count,7} = percWithNeighbour1; % percentage of red blobs associated
        results{count,8} = percWithNeighbour2; % percentage of red blobs associated
 
        
        count = count+1 ; 
  end 
  
  
end 


% concatenate rough size results into table as well 
        results_tab_sizech1 = array2table(size_results_ch1,'VariableNames',{'area'}); % sizes of 647 spots
        results_tab_sizech2 = array2table(size_results_ch2,'VariableNames',{'area'}); % sizes of 647 spots
         
%unused
%results_tab_ch1_min_nonAssoc = array2table(ch1_min_sizes_nonAssoc,'VariableNames',{'area'});

% concatenate sizes of particles and labels for their area ids
if size(sizes_labels,2) ~= 1
  results_size_label_ch1_assoc = table(sizes_labels',ch1_min_sizes_assoc,'VariableNames',{'id','area'});
else
  results_size_label_ch1_assoc = table(sizes_labels,ch1_min_sizes_assoc,'VariableNames',{'id','area'});
end

% this is the array of sizes of EAP45 particles which are neighbours to a
% colocalising particle
results_size_label_ch1_Nonassoc = table(sizes_labels_nonAssoc,ch1_min_sizes_NonAssoc_neighbours,'VariableNames',{'id','area'});

if save_flag;writetable(results_size_label_ch1_assoc,[output_dir_path filesep 'results_ch1_associated_size.csv']);end 
if save_flag;writetable(results_size_label_ch1_Nonassoc,[output_dir_path filesep 'results_ch1_nonAssociated_neighbours_size.csv']);end 

% Concatenate colocalisation results into a table
        results_tab = cell2table(results,'VariableNames',{'area_id','tot_red_blobs','tot_green_blobs','r_search','assoc_red_blobs','assoc_green_blobs','perc_assoc_red','perc_assoc_green'});
    
% save result in csv file for plotting later
 if save_flag;writetable(results_tab,[output_dir_path filesep 'results.csv']);end 
 if save_flag;writetable(results_tab_sizech1,[output_dir_path filesep 'results_ch1_size.csv']);end 
 if save_flag;writetable(results_tab_sizech2,[output_dir_path filesep 'results_ch2_size.csv']);end 


