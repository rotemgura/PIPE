function min_distance_to_boundary=PIPE_find_min_distance_to_boundary(boundary_mask,coords)
    %%% This function receives the cell mask and calculates the minimal 
    %%% distance between the center of the photo-conversion pulse and the
    %%% the cell boundary along the pre-defined direction of the intensity profile
    
    %%% unpack coordinates
    temp=num2cell(coords);
    [~,y0,x0,direction,thickness]=temp{:};
    
    %%% define area of intensity profile
    [profile_mask,perp_direction]=PIPE_define_1d_profile_mask(y0,x0,direction,thickness,size(boundary_mask)); 
    mask=profile_mask.*boundary_mask;
    
    %%% find boundaries
    beginning_boundary=find(sum(mask,perp_direction)>thickness/2,1);
    end_boundary=find(sum(mask,perp_direction)>thickness/2,1,'last');
    
    %%% find distance from boundaries
    distance1=coords(3-perp_direction)-beginning_boundary; % if perp_direction is y then distance is along x and vice versa
    distance2=end_boundary-coords(3-perp_direction);
    min_distance_to_boundary=min(distance1,distance2);
end
    
    