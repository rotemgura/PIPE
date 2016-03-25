function [t_start, y0, x0] = PIPE_find_pulse_center(tiffs,t_start,apply_mask)
    %%% This function receives a 3D movie matrix and finds the frame where
    %%% the photo-conversion pulse starts and the center of the
    %%% photo-converted protein ensemble.
    %%% Arguments:
    %%% tiffs -- a 3D movie matrix
    %%%
    %%% t_start -- the frame in which the photo-conversion pulse starts.
    %%% Should be 0 if not known 
    %%%
    %%% apply_mask -- Boolean. If the movie is imaged within a cell with visible
    %%% boundaries then this flag ignores the signal outside the cell
    %%%
    %%% Output
    %%% t_start -- the frame in which the photo-conversion pulse starts.
    %%%
    %%% y0 -- y coordinate of the pulse center at frame t_start
    %%%
    %%% x0 -- x coordinate of the pulse center at frame t_start
    
    %%% set constants
    MIN_LAG_BEFORE_PULSE=10;
    GRAD_RATIO_THRESH=5;
    MAX_AREA_CUTOFF=pi*100^2;
    MIN_AREA_CUTOFF=pi*15^2;
    MAX_ITERS=110;
    MIN_MOVIE_LENGTH=30;
    
    %%% segment image
    if apply_mask 
        mask=PIPE_identify_cell(mean(tiffs(:,:,end-MIN_MOVIE_LENGTH:end),3)); %identify the cell from data at the end of movie
    else
        mask=ones(size(tiffs,1),size(tiffs,2));
    end
    
    %%% find the time point when the photo-conversion pulse starts (denoted t_start) 
    %%% by looking for a suddent increase in intensity gradient
    grad_length=1;
    while t_start<1 && grad_length<5
        avg_int_array=arrayfun(@(x)mean(mean(immultiply(mask,tiffs(:,:,x)))),1:size(tiffs,3)); %calculate mean intensity in the image
        grad_array=avg_int_array(1+grad_length:end)-avg_int_array(1:end-grad_length); %calculate intensity gradient
        for i=MIN_LAG_BEFORE_PULSE:size(tiffs,3)-grad_length
            if grad_array(i)/std(grad_array(1:i-1))>GRAD_RATIO_THRESH %gradient is much steeper at point i than before i
                t_start=i+grad_length;
                break
            end
        end
        grad_length=grad_length+1; %steep gradient was not detected, so calculate the gradient with a longer delay to ignore short fluctuations
    end
    if t_start==0 %failed to detect steep gradient
        y0=0;
        x0=0;
        return
    end
    
    %%% find the pulse center in image t_start
    pulse_image = immultiply(mask,tiffs(:,:,t_start))-mean(tiffs(:,:,1:5),3); %A pre-photo-conversion image is subtracted from the examined image to minimize effects of existing red signal
    found_center=0;
    thresh=1.0;
    iter_num=1;
    
    %%% detect pulse center by thresholding
    while found_center==0 && iter_num<MAX_ITERS
        
        %find pixels above threshold
        above_thresh_image = zeros(size(pulse_image));
        above_thresh_image(pulse_image>thresh*max(pulse_image(:)))=1;
        
        %find islands above threshold
        dilated = imdilate(above_thresh_image,strel('disk',2));
        filled = imfill(dilated,'holes');
        eroded = imerode(filled,strel('disk',3));
        labeled = bwlabel(eroded);
        
        %choose largest island
        stats = regionprops(labeled,'Area','Centroid');
        areas = [stats.Area];
        
        %check if island area is within boundaries
        if isempty(areas) % decrease threshold
            thresh=thresh-0.1;
        else
            [max_area, ind] = max(areas);
            if max_area<MIN_AREA_CUTOFF % decrease threshold
                thresh=thresh-0.1;
            elseif max_area>MAX_AREA_CUTOFF % increase threshold
                thresh=thresh+0.001;
            else % found center of pulse
                x0 = round(stats(ind).Centroid(1));
                y0 = round(stats(ind).Centroid(2));
                found_center=1;
            end
        end
        iter_num=iter_num+1;
    end
    if iter_num==MAX_ITERS %failed to find center of pulse
        x0=0;
        y0=0;
    end
end