function [t0,t_start,y0,x0]=PIPE_detect_pulse_end(tiffs,t_start,y0_start,x0_start)
    %%% This function receives the x,y,t coordinates of the center of the PhoC
    %%% pulse and determines the time of the end of the pulse by tracking the
    %%% total fluorescence at the area of the pulse. It returns the time of the
    %%% end of the pulse and the center of the pulse in each frame during the
    %%% time when the pulse is on.
    %%%
    %%% Arguments:
    %%% tiffs -- 3D movie matrix
    %%%
    %%% t_start -- frame in which the photo-conversion pulse starts
    %%%
    %%% y0_start -- y coordinate of the pulse center at frame t_start
    %%%
    %%% x0_start -- x coordinate of the pulse center at frame t_start
    %%%
    %%% Output:
    %%% t0 -- frame in which the photo-conversion pulse ends
    %%%
    %%% t_start -- frame in which the photo-conversion pulse starts
    %%%
    %%% y0 -- an array of y coordinate of pulse center for each time point
    %%% between t_start and t0
    %%%
    %%% x0 -- an array of x coordinate of pulse center for each time point
    %%% between t_start and t0
    
    THRESH=150;
    SEGMENT_LENGTH=15;
    ROI_RADIUS=30;
    %%% check data analyzability
    if size(tiffs,3)<t_start+30 || t_start<1 || x0_start<ROI_RADIUS || y0_start<ROI_RADIUS || size(tiffs,1)-y0_start<ROI_RADIUS || size(tiffs,2)-x0_start<ROI_RADIUS
        t0=0;
        x0=0;
        y0=0;
        return
    end
    
    %%% find duration of photo-conversion pulse
    %define an ROI around pulse center
    roi_mask = zeros(size(tiffs(:,:,1)));
    disk = strel('disk',ROI_RADIUS);
    roi_mask(y0_start-ROI_RADIUS+1:y0_start+ROI_RADIUS-1,x0_start-ROI_RADIUS+1:x0_start+ROI_RADIUS-1)=disk.getnhood();
    
    %calculate mean intensity in ROI vs time
    data=arrayfun(@(x)mean(mean(immultiply(tiffs(:,:,x),roi_mask))),t_start:size(tiffs,3));
    
    %calculate duration of photo-conversion pulse
    duration=PIPE_find_phoc_peak(data);
    if duration==-1 %failed to calculate duration (data might be too noisy) 
        [t0,y0,x0]=deal(0);
        return
    end
    t0=t_start-1+duration;

    %find the center of the pulse in all frames between the pulse's beginning & end
    x0=zeros(1,t0-t_start+1);
    y0=zeros(1,t0-t_start+1);
    for i=t_start:t0
        [~,y,x]=PIPE_find_pulse_center(tiffs,i,0);
        x0(i-t_start+1) = x;
        y0(i-t_start+1) = y;
    end
end