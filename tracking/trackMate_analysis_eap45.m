% Script to analyze two-colour live cell movies taken on the SIM in
% widefield TIRF mode from fluorescently labelled HIV virus and the EAP45
% protein. The script takes as input:

%       - 2 tiff stacks representing a ROI in the two channels (488 and 640) 
%       - *_Tracks.xml file exported from Trackmate for each of the two tiff stacks 
 
% The output is a stabilized image of the cyan channel (a rectangular ROI is cropped
% the center x,y coordinates of the cyan track), with an outline of the
% magenta tracks over it. 
% The script also now outputs a plot of the distance between the two
% particle centroids vs time, along with a histogram of the number of
% frames for all distances between particles. 

% For tracking HIV-ESCRTII data sets
% Created: 15/10/2019

% now accounting for tracks in which the frame numbers are not equivalent
% and only overlap for a portion of the tracks. 

clear all
close all


%% Set up directories

% Images must be in the same directory
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live8_slide1well2_1_ROI2';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live8_slide1well2_1_ROI1';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live8_slide1well2_4_ROI1';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live8_slide1well2_4_ROI3';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live7_slide1well2_8_ROI3';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live7_slide1well2_8_ROI4';
PathInput         = 'D:\OneDriveNew\OneDrive - University of Cambridge\lag\microscopy work\hiv-eap45\tracking data\live7_slide1well2_8_ROI2';

% on the work pc
PathInput         = 'D:\OneDriveNew\OneDrive - University of Cambridge\lag\microscopy work\hiv-eap45\tracking data\live8_slide1well2_1_ROI2';


% % Inputs from Live6_ROI12 
PathInput         = 'Y:\Users\ppv23\2019_08_15_Pedro_Bo_Katharina_LiveExperiment6\tracking';
PathInput         = 'Y:\Users\ppv23\2019_08_15_Pedro_Bo_Katharina_LiveExperiment6\tracking\live6_slide1well3_12_ROI2'; % not useful, too close to edge
PathInput         = 'Y:\Users\ppv23\2019_08_15_Pedro_Bo_Katharina_LiveExperiment6\tracking\live6_slide1well3_12_ROI3';
PathInput         = 'Y:\Users\ppv23\2019_08_15_Pedro_Bo_Katharina_LiveExperiment6\tracking\live6_slide1well3_12_ROI5';
  
% Inputs from Live7 
PathInput         = 'Y:\Users\ppv23\2019_09_05_Pedro_Bo_Katharina_liveExperiment7\ROIs_Katharina\live7_slide1well2_4_ROI1';
PathInput         = 'Y:\Users\ppv23\2019_09_05_Pedro_Bo_Katharina_liveExperiment7\ROIs_Katharina\live7_slide1well2_4_ROI2'; %no video yet 
%PathInput         = 'Y:\Users\ppv23\2019_09_05_Pedro_Bo_Katharina_liveExperiment7\ROIs_Katharina\live7_slide1well2_5_ROI1';
%PathInput         = 'Y:\Users\ppv23\2019_09_05_Pedro_Bo_Katharina_liveExperiment7\ROIs_Katharina\live7_slide1well2_7_ROI1';

% for the paper figure
PathInput = 'D:\OneDriveNew\OneDrive - University of Cambridge\lag\microscopy work\hiv-eap45\tracking data\live7_slide1well2_8_ROI3';
PathInput         = 'D:\OneDriveNew\OneDrive - University of Cambridge\lag\microscopy work\hiv-eap45\tracking data\live8_slide1well2_1_ROI1';
PathInput         = '/Users/pedrovallejo/OneDrive - University Of Cambridge/lag/microscopy work/hiv-eap45/tracking data/live7_slide1well2_8_ROI3';

[filename,Path]   = uigetfile([PathInput filesep '*.tif']); % pick the cyan file
ch1               = '_488';
ch2               = '_640';
filename_magenta  = strsplit(filename,'_');
filename_magenta  = [filename_magenta{1},'_',filename_magenta{2},'_',filename_magenta{3},ch2,'_',filename_magenta{5},'_',filename_magenta{6}];
info              = imfinfo([Path filename]);
N_frames          = length(info);

% Change time frame and pixel value depending on what the original image
% values are - if it's microns and time or pixels and frame
dT                  = 5.22; % time interval between frames 5.22
pix                 = 1; % pixel size in microns 0.05
edge                = 14; % pick 10 pixels on each side of the centroid to crop the image

% turn on or off intermediate output steps
flag                = 1;
show                = flag;
video_flag          = flag;
show_original       = flag;
video_flag_original = flag;

% Variables for making nice figures in matlab
width               = 10;  % width in centimeters
height              = 10;  % height in centimeters
font_size           = 12;
lw                  = 1.0; %linewidth
msz                 = 4;  %markersize


%% Import the tracks from TrackMate and account for any skipped frames

% import tracks automatically using the name of the 488 channel
if ismac
    addpath('/Applications/Fiji.app/scripts'); % to obtain function ImportTrackMateTracks
else
    addpath('D:\Fiji.app\scripts');
end

track_filename_cyan     = strsplit(filename,'.');
track_filename_magenta  = strsplit(filename_magenta,'.');
path_tracks_cyan        = [Path track_filename_cyan{1} '_Tracks.xml'];
path_tracks_magenta     = [Path track_filename_magenta{1} '_Tracks.xml'];
clipZ                   = true; % Remove Z coordinates, if you know you can.
scaleT                  = true; % Use physical time for T.

[tracks_cyan,md]        = importTrackMateTracks(path_tracks_cyan,clipZ, scaleT);
[tracks_magenta,md]     = importTrackMateTracks(path_tracks_magenta,clipZ, scaleT);

%tracks_cyan{1}(5, :); % this is how you index into the tracks

n_tracks_cyan         = numel( tracks_cyan );
n_tracks_magenta      = numel( tracks_magenta );

fprintf('Found %d tracks in the cyan channel.\n', n_tracks_cyan)
fprintf('Found %d tracks in the magenta channel.\n', n_tracks_magenta)

% choose track and turn into table
% Add 1 to the frames to avoid indexing from zero
track1_cyan     = tracks_cyan{1};
tr_cyan         = [(track1_cyan(:,1)/dT)+1,track1_cyan(:,2),track1_cyan(:,3)];
tr_cyan         = array2table(tr_cyan,'VariableNames',{'frame','x','y'});

track1_magenta  = tracks_magenta{1};
tr_magenta      = [(track1_magenta(:,1)/dT)+1,track1_magenta(:,2),track1_magenta(:,3)];
tr_magenta      = array2table(tr_magenta,'VariableNames',{'frame','x','y'});

% Check if there are missing frames
% if the frame numbers skip locations, duplicate the last frame before
% movement or the next frame to advance the movement, depending on which
% one is more representative to what you observe in the videos. 

% PENDING WORK: INTERPOLATE BETWEEN POSITIONS OF FRAME JUMPS TO CREATE
% SMOOTHER VISUALIZATIONS, INSTEAD OF FILLING IN WITH PREVIOUS/POST
% POSITIONS.

cyan_arr            = table2array(tr_cyan);
magenta_arr         = table2array(tr_magenta);
cyan_arr_new        = [];
magenta_arr_new     = [];

tr_cyan             = consecutiveCheck(cyan_arr);
tr_magenta          = consecutiveCheck(magenta_arr);


%% Make folders to store videos

if video_flag
    file = strsplit(filename,'.');
    video_dir_path = fullfile(Path,sprintf('video_%s',file{1}));
    if exist(video_dir_path, 'dir')
        opts.Interpreter = 'tex';
        opts.Default = 'Continue';
        quest = '\fontsize{12}An output folder ''video'' already exists. If you continue, data in this folder might be overwritten.';
        answer = questdlg(quest,'Message','Cancel','Continue',opts);
    else
        mkdir(video_dir_path);
    end

end 

if video_flag_original
    file = strsplit(filename,'.');
    video_dir_path_orig = fullfile(Path,sprintf('video_cropRect_%s',file{1}));
    if exist(video_dir_path_orig, 'dir')
        opts.Interpreter = 'tex';
        opts.Default = 'Continue';
        quest = '\fontsize{12}An output folder ''video'' already exists. If you continue, data in this folder might be overwritten.';
        answer = questdlg(quest,'Message','Cancel','Continue',opts);
    else
        mkdir(video_dir_path_orig);
    end

end 

%% Read images, create visualizations of the tracks, and compute distance between particles in every frame

nframes_cyan    = numel(tr_cyan.frame);
nframes_magenta = numel(tr_magenta.frame);

% array to hold cropped images
im_stable       = cell(nframes_cyan,1);
x_magenta       = zeros(nframes_magenta,1);
y_magenta       = zeros(nframes_magenta,1);
d               = zeros(nframes_magenta,1);
frame_array     = zeros(nframes_magenta,1); 
common_counter = 1; % counter to pick the first commmon frame between cyan and magenta

% Loop over the frames in the EAP45 (moving) channel to generate pictures
% and quantify the distance between the two particles. 

for i = 1:nframes_magenta
    
    % start from the first matching track between the magenta and the cyan,
    % and stop at the last frame of the magenta tracks.
    
    frame_id_magenta  = uint8(tr_magenta.frame(i));    
    frame_id_cyan     = find(uint8(tr_cyan.frame) == frame_id_magenta); %corresponding cyan frame number to the first point in the magenta tracks
    
    % if matching track not found on the first magenta track id, continue
    % iterating until found 
    
    if isempty(frame_id_cyan)
        continue
    elseif common_counter
        first_common_frame_magenta = i; % first frame in the magenta channel which matches a cyan channel
        first_common_frame_cyan = frame_id_cyan;
        common_counter = 0;
    end 
    
    max_frame_magenta = uint8(max(tr_magenta.frame));
    
    % continue the loop until the last frame of the magenta tracks
    if frame_id_cyan > max_frame_magenta
        return
    end 
    
    % read in the two mages
    im_cyan      = imread([Path filename], frame_id_magenta);
    im_magenta   = imread([Path filename_magenta], frame_id_magenta);
    
    x_cyan = tr_cyan.x(frame_id_cyan)/pix;
    y_cyan = tr_cyan.y(frame_id_cyan)/pix;
    
    % if x or y coordinates from cyan are too close to edge, skip the frame
    % and proceed to the next iteration of the loop
    if x_cyan < edge
        continue
    elseif y_cyan < edge
        continue 
    end 
        
    % recenter frame based on centroid of the cyan channel - crop using
    % imcrop
    crop_rect      = [ceil(x_cyan) - edge, ceil(y_cyan) - edge, 2*edge+1,2*edge+1];
    im_stable{i,1} = imcrop(im_cyan,crop_rect);
   
    % transform magenta tracks into new coordinate system
    
    %matching_id = find(uint8(tr_magenta.frame) == frame_id);
    
    x_magenta(i) = (tr_magenta.x(i)/pix) - (ceil(x_cyan) - edge);
    y_magenta(i) = (tr_magenta.y(i)/pix) - (ceil(y_cyan) - edge);
    
    % show images to check green and red spots on cyan picture 
    if show_original 
        figure;
        imshow(im_cyan,[],'InitialMagnification',1500);
        hold on; 
        rectangle('Position',crop_rect,'EdgeColor', 'c', 'LineWidth', 4, 'LineStyle','--');
        plot(x_cyan,y_cyan,'+','color','c','MarkerSize',20,'LineWidth',4);
        plot(tr_magenta.x(i)/pix,(tr_magenta.y(i)/pix),'o','color','m','MarkerSize',20,'LineWidth',8);
        frame(i) = getframe(gcf);
        fig2 = gcf;        
    end 
    
     if video_flag_original
        %print(fig2,[video_dir_path_orig filesep 'tracking_' num2str(i)],'-dpng','-r300') % print at 300 dpi resolution
        f=getframe;
        imwrite(f.cdata,[video_dir_path_orig filesep 'tracking_' num2str(i) '.png']) % print at 300 dpi resolution
     end
    
    
    if show
        figure
        imshow(im_stable{i},[],'InitialMagnification',1500); 
        hold on
        plot(edge+1,edge+1,'o','color','c','MarkerSize',20,'LineWidth',4);
        plot(ceil(size(im_stable{i},1)/2),ceil(size(im_stable{i},1)/2),'+','color','c','MarkerSize',20,'LineWidth',4);
        plot(x_magenta(i),y_magenta(i),'o','color','m','MarkerSize',25,'LineWidth',4);
        frame(i) = getframe(gcf);
        fig2 = gcf;     
    end 
    
    if video_flag
        %print(fig2,[video_dir_path filesep 'tracking_' num2str(i)],'-dpng','-r300') % print at 300 dpi resolution
        f=getframe;
        imwrite(f.cdata,[video_dir_path filesep 'tracking_' num2str(i) '.png']) % print at 300 dpi resolution
    end
    
    close(gcf)
    
    % now use the raw tracks for quantification of mean distance between
    % centroids
    
    d(i) = sqrt( (x_cyan - tr_magenta.x(i))^2 + (y_cyan - tr_magenta.y(i))^2);
    frame_array(i) = frame_id_magenta;

end

% outputs from this for loop are:
% - Shifted x and y coordinates for the magenta channel
% - distance between x and y as a function of the frame number
% - stabilized cropped image of the cyan channel

%% Make nice visualizations of the data

% in the case in which the magenta tracks leave the cropped rectangle
% check when the magenta coordinates become negative in the x or y domain

mask_inside_rectangle_y = y_magenta > 0;
mask_inside_rectangle_x = x_magenta > 0;

sum_y                   = sum(mask_inside_rectangle_y);
sum_x                   = sum(mask_inside_rectangle_x);

if sum_y < sum_x
    mask_inside_rectangle = mask_inside_rectangle_y;
    fprintf('Magenta particle leaves the rectangle in the Y coordinate \n')

else 
    mask_inside_rectangle = mask_inside_rectangle_x;
    fprintf('Magenta particle leaves the rectangle in the X coordinate \n')
end 

% get a logical mask with these indices
% use the logical mask to select the tracks to plot over the average Gag
% image

% plot lines of the magenta track over a representative image of the Gag
% remove tracks with zeros in them for the magenta x and y coordinates

x_magenta_cropped = x_magenta(mask_inside_rectangle);
y_magenta_cropped = y_magenta(mask_inside_rectangle);


% If any frames were removed due to the xy coordinates of the cyan channel
% being close to the edges, remove these frames as well from the xy magenta
% coordinates, the distance between x and y, and the stable image. 

% indicex of xy coordinates of cyan channel inside box
inside_box = find(~cellfun(@isempty,im_stable));

% obtain average Gag image 
im_stable_old = im_stable;
% remove empty cells in im_stable
im_stable = im_stable(inside_box);
% (only consider frames in common between cyan and magenta channel)
%im_stable = im_stable(first_common_frame_magenta:end,:);

for i = 1: size(im_stable,1)
    if i ==1
        avg_im = im_stable{i};
    else 
        avg_im = im_stable{i} + avg_im;
    end 
    
end

avg_im = avg_im/size(im_stable,1);

% this doesn't perfectly need to correspond with the magenta channel, since
% it's only to get the average gag image

% to plot the one-colored track with jet-colormap points over it in a
% stationary image

figure
c = jet(size(x_magenta_cropped,1));
imshow(avg_im,[],'InitialMagnification',3000);
hold on
%plot(x_magenta_cropped,y_magenta_cropped,'-','Color', [0.8500 0.3250 0.0980], 'LineWidth',2)
lwd = 1.5;
for j = 1:size(x_magenta_cropped,1)
    if j == 1 
        plot(x_magenta_cropped(1:j),y_magenta_cropped(1:j),'-','Color', c(j, :), 'LineWidth',lwd)
    elseif j ~= size(x_magenta_cropped,1)
        plot(x_magenta_cropped(j-1:j+1),y_magenta_cropped(j-1:j+1),'-','Color', c(j, :), 'LineWidth',lwd)
    else
        plot(x_magenta_cropped(j-1:j),y_magenta_cropped(j-1:j),'-','Color', c(j, :), 'LineWidth',lwd)
    end 
end 

hold off
fig3 = gcf;
%print(fig3,[Path filesep track_filename_cyan{1} '_Co-moving_frame_wTracks'],'-dpng','-r300') % print at 300 dpi resolution
f=getframe;
imwrite(f.cdata,[Path filesep track_filename_cyan{1} '_Co-moving_frame_wTracks.png']) % print at 300 dpi resolution

% % print colorbar for jet colormap
% figure
% colormap(c)
% colorbar
% fig_colorbar = gcf;
% f=getframe;
% imwrite(f.cdata,[Path filesep track_filename_cyan{1} '_Co-moving_frame_wTracks_colorbar.pdf']) % print at 300 dpi resolution

%% Quantification

% get the mean distance between the raw cyan and raw magent tracks, and
% determine the duration the magenta particle spends within 2 pixels of
% the Gag particle.

% get distances only from the common frames between magenta and cyan
d_old           = d;
frame_array_old = frame_array;

%d           = d(first_common_frame_magenta:end);
%frame_array = frame_array(first_common_frame_magenta:end);

% remove elements in which the xy coords of the cyan channel were at the
% edges
d           = d(inside_box);
frame_array = frame_array(inside_box);

associated  = d(d<2);
% open a summary file for writing results
summary_file = fopen(fullfile(Path,strcat(track_filename_cyan{1},'_summary.txt')),'wt');
fprintf(summary_file,"Total number of frames the magenta particle stays within 2 pixels of the Gag particle is %d \n",size(associated,1));

% Check how many consecutive frames the EAP45 stays within 2 pixels of the
% Gag
% Average number of consecutive frames it stays within 2 pixels'

assoc_mask = d<2; % mask for indices of all instances distance between regions is smaller than 2
count      = 0;
avg        = zeros(size(assoc_mask ,1),1); % initialize array to hold some stuff

for i = 1: size(assoc_mask ,1)
    consecutive = assoc_mask(i);
    if consecutive == 1
        count = count+1;
    else
        count = 0;
    end 
    avg(i) = count;
end

max_fr = max(max(avg)); % Max number of consecutive frames the EAP45 is within 2 pixels of Gag
fr_sum = sum(avg ~= 0); % Number of instances in which the EAP45 is within 2 pixels of the Gag 
avg_fr = ceil(sum(avg)/fr_sum); % Average number of consecutive frames the EAP45 stays within 2 pixs of the Gag

fprintf(summary_file,"The magenta particle stays within 2 pixels of the Gag particle a maximum of %d consecutive frames \n",max_fr);
fprintf(summary_file,"The magenta particle stays within 2 pixels of the Gag particle on average for %d consecutive frames \n",avg_fr);


% Plot the # of consecutive frames the EAP45 stays close to the Gag 

% histogram of consecutive frames 
figure;
histogram(avg); 

% plot of distance from Gag as a function of time
time_per_frame = 5.22;
figure;
plot(frame_array*time_per_frame,d)
hold on
xlabel('Time (s)');
ylabel('Distance between EAP45 and Gag centroid (pixels)');
yticks([0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
plot(frame_array*time_per_frame,repmat(2,size(frame_array,1)),'--')
set(gca,...
                'FontSize',font_size,...
                'Units','Normalized',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontName','Arial');

fig4 = gcf;          
print(fig4,[Path filesep track_filename_cyan{1} '_distance_vs_time'],'-dpng','-r300') % print at 300 dpi resolution

figure; histogram(d);
xlabel('Distance between EAP45 and Gag centroid (pixels)');
ylabel('Number of candidates');
xticks([0 2 4 6 8 10 12 14 16 18 20 22 24 26 28 30])
set(gca,...
                'FontSize',font_size,...
                'Units','Normalized',...
                'FontUnits','points',...
                'FontWeight','normal',...
                'FontName','Arial');

fig5 = gcf;          
print(fig5,[Path filesep track_filename_cyan{1} '_distance_vs_time_histogram'],'-dpng','-r300') % print at 300 dpi resolution

%% Print out distance vs time data to plot in R
% Concatenate results into a table
results = [frame_array,d];
results_tab = array2table(results,'VariableNames',{'frame','distance'});
%results = table(X.pearsonCoeff,X.spearmanCoeff,X.mandersCoeff,X.M1,X.M2,'VariableNames',{'PC','SPC','MoC','M1RG','M2GR'});

% save result in csv file for plotting later
writetable(results_tab,[Path filesep track_filename_cyan{1} 'distance_vs_time_results.csv']);

%% Deprecated code

% Take a Fourier transform of the distance plot to see if there are any
% signicant frequencies of oscillation of the EAP45
% 
% Fs = 1/time_per_frame;            % Sampling frequency                    
% T = 1/Fs;             % Sampling period       
% L = 211;             % Length of signal
% t = (0:L-1)*T;        % Time vector
% 
% Y = fft(d);
% L = max(max(frame_array));
% 
% P2 = abs(Y/L);
% P1 = P2(1:L/2+1);
% P1(2:end-1) = 2*P1(2:end-1);
% 
% figure;
% f = Fs*(0:(L/2))/L;
% plot(f,P1) 
% title('Single-Sided Amplitude Spectrum of X(t)')
% xlabel('f ')
% ylabel('|P1(f)|')

% Plotting the track dynamically over the picture of the co-moving frames
% You can use this to visualize dynamically how the magenta particle moves
% over the cyan channel. 

% figure
% imshow(avg_im,[],'InitialMagnification',2000);
% hold on
% for j = 1:size(x_magenta_cropped,1)
%  plot(x_magenta_cropped(1:j),y_magenta_cropped(1:j),'>','Color', 'y', 'LineWidth', 2)
%  pause
% end 
% hold off


%% old plot tracks 
% figure
% hold on
% c = jet(n_tracks_cyan);
% for s = 1 : n_tracks_cyan
% x_cyan = tracks_cyan{s}(:, 2);
% y_cyan = tracks_cyan{s}(:, 3);
% plot(x_cyan, y_cyan, '.-', 'Color', c(s, :))
% end
% axis equal
% xlabel( [ 'X (' md.spaceUnits ')' ] )
% ylabel( [ 'Y (' md.spaceUnits ')' ] )

%% old code to check if consecutive frames are being skipped

% counter             = 2;
% consecutive_counter = 1;
% 
% % check if any frames are being skipped or jumped in the cyan channel
% for i = 2:size(cyan_arr,1)
%     
%     % don't do anything for the first frame
%     % check if next frame is immediately following the previous
%     prev_fr     = cyan_arr(i-1,1,1);
%     fr          = cyan_arr(i,1,1);
%     consecutive = fr-prev_fr;
%     
%     % keep track of the frame numbers added with the counter 
%     if uint8(consecutive) ~= 1
%        % insert a copy of the previous row to fill in the gap 
%        %copy_row = repmat(cyan_arr(i-1,:,:),[consecutive-1 1]); % copy of the previous row
%        copy_row = repmat(cyan_arr(i,:,:),[consecutive-1 1]); 
%        % add frame numbers to fill consecutive gap
%        %for j = 1:size(copy_row,1)
%        k = 1;
%        for j = size(copy_row,1):-1:1        
%            %copy_row(j,1,1) = copy_row(j,1,1)+j;
%            copy_row(k,1,1) = copy_row(k,1,1)-j;
%            k = k+1;
%        end 
%        % append new rows to the previous 
%        cyan_arr_new         = cat(1,cyan_arr_new(1:counter-1,:,:),copy_row);
%        counter              = counter+consecutive-1;
%        consecutive_counter  = 1;
%     else 
%         % add the two consecutive rows after a gap
%         % if there is no gap, add a single row       
%         copy_row = cyan_arr(i-1:i,:,:);
%         if counter == 2 %for the first iteration
%             cyan_arr_new = copy_row;
%             counter      = counter+1;
%         else 
%             if consecutive_counter > 1 % after two consecutive runs, only add one row
%                 copy_row     = cyan_arr(i,:,:);
%                 cyan_arr_new = cat(1,cyan_arr_new(1:counter-1,:,:),copy_row);
%                 counter      = counter+1;
%             else 
%                 cyan_arr_new = cat(1,cyan_arr_new(1:counter-1,:,:),copy_row);
%                 counter      = counter+2;
%             end 
%         end 
%         consecutive_counter = 1+consecutive_counter;
%     end 
%        
% end 
% 
% % compare to original array - if different, make new tables
% if size(cyan_arr,1) ~= size(cyan_arr_new,1)
%     tr_cyan     = array2table(cyan_arr_new,'VariableNames',{'frame','x','y'});
%     
% end 
% 
% % check if any frames are being skipped or jumped in the magenta channel
% counter             = 2;
% consecutive_counter = 1;
% for i = 2:size(magenta_arr,1)
%     
%     % don't do anything for the first frame
%     % check if next frame is immediately following the previous
%     prev_fr     = magenta_arr(i-1,1,1);
%     fr          = magenta_arr(i,1,1);
%     consecutive = fr-prev_fr;
%     
%     % keep track of the frame numbers added with the counter 
%     if uint8(consecutive) ~= 1
%        % insert a copy of the previous row to fill in the gap 
%        %copy_row = repmat(cyan_arr(i-1,:,:),[consecutive-1 1]); % copy of the previous row
%        copy_row = repmat(magenta_arr(i,:,:),[consecutive-1 1]); 
%        % add frame numbers to fill consecutive gap
%        %for j = 1:size(copy_row,1)
%        k = 1;
%        for j = size(copy_row,1):-1:1        
%            %copy_row(j,1,1) = copy_row(j,1,1)+j;
%            copy_row(k,1,1) = copy_row(k,1,1)-j;
%            k = k+1;
%        end 
%        % append new rows to the previous 
%        magenta_arr_new         = cat(1,magenta_arr_new(1:counter-1,:,:),copy_row);
%        counter              = counter+consecutive-1;
%        consecutive_counter  = 1;
%     else 
%         % add the two consecutive rows after a gap
%         % if there is no gap, add a single row       
%         copy_row = magenta_arr(i-1:i,:,:);
%         if counter == 2 %for the first iteration
%             magenta_arr_new = copy_row;
%             counter      = counter+1;
%         else 
%             if consecutive_counter > 1 % after two consecutive runs, only add one row
%                 copy_row     = magenta_arr(i,:,:);
%                 magenta_arr_new = cat(1,magenta_arr_new(1:counter-1,:,:),copy_row);
%                 counter      = counter+1;
%             else 
%                 magenta_arr_new = cat(1,magenta_arr_new(1:counter-1,:,:),copy_row);
%                 counter      = counter+2;
%             end 
%         end 
%         consecutive_counter = 1+consecutive_counter;
%     end 
%        
% end 
% 
% % compare to original array - if different, make new tables
% if size(magenta_arr,1) ~= size(magenta_arr_new,1)
%     tr_magenta     = array2table(magenta_arr_new,'VariableNames',{'frame','x','y'});
%     
% end 
