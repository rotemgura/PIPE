function profile=PIPE_apply_profile_mask(image, mask, perp_direction, thickness)
    %%% This function receives an image and a mask that represents a 1D 
    %%% profile passing through the image. The function returns a vector of
    %%% values from the original image along the profile represented by the 
    %%% mask, after some averaging along the perpendicular dimension.
    %%%
    %%% Arguments:
    %%% image - a 2D matrix of intensity values
    %%%
    %%% mask - a 2D boolean matrix of the same size as image
    %%%
    %%% perp_direction - the effective direction of averaging perpendicular  
    %%% to the main direction of the 1D profile (2 or 1 for x or y)
    %%%
    %%% thickness - # of pixels to average over perpendicular to the line
    
    
    %%% Prepare data
    % multiply image and profile mask
    filtered_image = immultiply(image, mask);
    
    % determine beginning and end points of profile (profile should be thick enough and not extend beyond cell boundaries)
    beginning_from_mask=find(sum(mask,perp_direction)>thickness/2,1);
    end_from_mask=find(sum(mask,perp_direction)>thickness/2,1,'last');
    
    beginning_from_data=find(sum(logical(filtered_image),perp_direction)>thickness/2,1);
    end_from_data=find(sum(logical(filtered_image),perp_direction)>thickness/2,1,'last');
    
    from=max(beginning_from_mask,beginning_from_data);
    to=min(end_from_mask,end_from_data);
    
    %%% average over thickness of profile
    collapse_perp_direction=sum(filtered_image,perp_direction); 
    normalization=sum(logical(filtered_image),perp_direction);
    profile = collapse_perp_direction(from:to)./normalization(from:to);
    
    if sum(isnan(profile))>0
        ['Detected a problem in continuity of intensity profile. Try running PIPE with ' ...
            'exclude_extra_cellular=0, or increase the threshold for cell segmentation']
        profile='';
    end
    
    %%% make sure profile is a column vector
    profile=profile(:); 
end