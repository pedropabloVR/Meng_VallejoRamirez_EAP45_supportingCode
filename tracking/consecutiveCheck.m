function [table_arr] = consecutiveCheck(arr1)
%consecutiveCheck takes in an array containing the tracks (frame, x, and y)
%for single particle tracking data, and it checks for non-consecutive
%frames in the tracks. If it finds any jumps between frames, it fills in
%the gap with the previous or following location of the particle. The
%output of the code is a table containing the frame, x, and y coordinates
%of the filled-in track. 
%   Pedro Vallejo Ramirez
%   14/10/2019

counter             = 2;
consecutive_counter = 1;

% check if any frames are being skipped or jumped in the cyan channel
for i = 2:size(arr1,1)
    
    % don't do anything for the first frame
    % check if next frame is immediately following the previous
    prev_fr     = arr1(i-1,1,1);
    fr          = arr1(i,1,1);
    consecutive = fr-prev_fr;
    
    % keep track of the frame numbers added with the counter 
    if uint8(consecutive) ~= 1
        consecutive = uint8(consecutive);
       % insert a copy of the previous row to fill in the gap 
       %copy_row = repmat(arr1(i-1,:,:),[consecutive-1 1]); % copy of the previous row
       copy_row = repmat(arr1(i,:,:),[consecutive-1 1]); 
       % add frame numbers to fill consecutive gap
       %for j = 1:size(copy_row,1)
       k = 1;
       for j = size(copy_row,1):-1:1        
           %copy_row(j,1,1) = copy_row(j,1,1)+j;
           copy_row(k,1,1) = copy_row(k,1,1)-j;
           k = k+1;
       end 
       % append new rows to the previous 
       arr1_new         = cat(1,arr1_new(1:counter-1,:,:),copy_row);
       counter              = counter+consecutive-1;
       consecutive_counter  = 1;
    else 
        % add the two consecutive rows after a gap
        % if there is no gap, add a single row       
        copy_row = arr1(i-1:i,:,:);
        if counter == 2 %for the first iteration
            arr1_new = copy_row;
            counter      = counter+1;
        else 
            if consecutive_counter > 1 % after two consecutive runs, only add one row
                copy_row     = arr1(i,:,:);
                arr1_new = cat(1,arr1_new(1:counter-1,:,:),copy_row);
                counter      = counter+1;
            else 
                arr1_new = cat(1,arr1_new(1:counter-1,:,:),copy_row);
                counter      = counter+2;
            end 
        end 
        consecutive_counter = 1+consecutive_counter;
    end 
       
end 

% compare to original array - if different, make new tables
if size(arr1,1) ~= size(arr1_new,1)
    table_arr     = array2table(arr1_new,'VariableNames',{'frame','x','y'});
else 
    table_arr = array2table(arr1,'VariableNames',{'frame','x','y'});
end 

end

