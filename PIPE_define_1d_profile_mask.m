function [mask,perp_direction]=PIPE_define_1d_profile_mask(y0,x0,direction,thickness,mask_size)
    %%% This function creates a 2D mask of a 1D profile pre-specified angle 
    %%% and thickness, passing through a pre-specified point (y0,x0).
    %%% Angle 0 is defined as the positive x axis.
    %%% Arguments:
    %%% y0,x0 -- coordinates of a point through which the 1D profile passes
    %%%
    %%% direction -- an angle between -90 and 90, where 0 is positive x
    %%%
    %%% thickness -- number of lines to average over perpendicular to 'direction'
    %%%
    %%% mask_size -- size of 2D image in which the 1D profile would be defined
    %%%
    %%% Output:
    %%% mask -- a 2D boolean mask where ones define the 1D profile
    %%%
    %%% perp_direction -- the direction in which to average over lines. 
    %%% The value is either 1 or 2, corresponding to y or x.
    
    mask=zeros(mask_size);
    mask1=zeros(mask_size);
    t=-tan(direction*pi/180); % the minus sign compensates for Matlab's upside-down representation of images.
    
    % The function takes the positive x axis to be the forward direction, 
    % so if the angle is > 45 or <-45 then
    % flip x and y, construct the mask and then transpose the mask before returning
    if abs(t)<=1
        perp_direction=1;
        transpose=0;
    else
        perp_direction=2;
        transpose=1;
        tmp=x0;
        x0=y0;
        y0=tmp;
        mask_size=fliplr(mask_size);
        t=1/t;
    end

    %%% create the 1D profile in mask
    %The 1D profile is defined as y=t*x+B, where x0,y0 are given and t was already calculated above.
    B=y0-t*x0;
    
    if abs(1/t)>mask_size(2) % angle is sharp enough to define a vertical profile
        mask(y0-round(thickness/2):y0+thickness-round(thickness/2)-1,:)=1;
        if transpose
            mask=mask';
        end
        return
    end
    
    %determine the range of the x coordinate within the matrix boundaries
    if t>0 %the slope is positive
        if t+B>=1
            xmin=1;
        else
            xmin=round((1-B)/t);
        end
        if t*mask_size(2)+B>mask_size(1)
            xmax=round((mask_size(1)-B)/t);
        else
            xmax=mask_size(2);
        end
    else %the slope is negative
        if t+B<mask_size(1)
            xmin=1;
        else
            xmin=round((mask_size(1)-B)/t);
        end
        if t*mask_size(2)+B<1
            xmax=round((1-B)/t);
        else
            xmax=mask_size(2);
        end
    end
    
    % fill mask according to the functional definition of the line y(x)
    for x=xmin:xmax
        mask1(round(t*x+B),x)=1;
    end
    
    %%% broaden the line in the perpendicular direction
    %determine the thickness of the line in both directions
    go_back=-round(thickness/2)+1;
    go_forward=thickness+go_back-1;
    
    %consider only the perpendicular coordinate of the line and broaden
    %each point by moving back and forward along this coordinate.
    %first check constrains.
    line_locations=find(mask1);
    if line_locations(1)+go_back<1 %ignore the first point because broadening it would go outside the matrix
        line_locations(1)=line_locations(2);
    end
    if line_locations(end)+go_forward>numel(mask1) %ignore the last point because broadening it would go outside the matrix
        line_locations(end)=line_locations(end-1);
    end
    
    %now broaden the line
    for i=go_back:go_forward
        mask1(line_locations+i)=1;
    end
    
    %%% There are small islands of ones that were formed due to the cyclic
    %%% shift when broadedning the line. Copy the main body of ones into
    %%% the final mask
    bwlabeled=bwlabel(mask1);
    ind=bwlabeled(y0,x0); %this is the index of the main body of ones
    mask(bwlabeled==ind) = 1;
    
    % if the angle was >45 or <-45 we flipped x and y, so now transpose the mask
    if transpose
        mask=mask';
    end
end