function reduced=PIPE_reduce_n_time_points(time_points,max_n_time_points)
    %%% This function receives an array of time points and reduces the number
    %%% of points below the threshold max_n_time_points in a way that keeps the
    %%% density of time points as close as possible to uniform on a log scale
    
    %%% find threshold on the number of points in log scale intervals
    %bin the time points uniformly in log scale
    [h,centers]=hist(log(time_points(time_points>0)),min(10,max_n_time_points));
    
    % find threshold for the number of points to sample from each bin
    for i=1:max(h)
        if sum(min(h(:),ones(size(h(:)))*i))>max_n_time_points
            threshold=i-1;
            break
        end
    end
    
    %%% initialize screened index list. if 0 is included in original time 
    %%% points then add it now, because it does notnexist in the log scale
    %%% histogram above.
    if ismember(0,time_points)
        post_screen_indices=find(time_points==0);
    else
        post_screen_indices=[];
    end
    
    %%% dilute bins with many time points
    lower=exp([-1 (centers(1:end-1)+centers(2:end))/2]);
    upper=exp([(centers(1:end-1)+centers(2:end))/2 max(time_points)]);
    for i=1:numel(h)
        range_indices=find(time_points>lower(i) & time_points<=upper(i));
        if numel(range_indices)<=threshold % no need to dilute -- all time points in current bin pass.
            post_screen_indices=[post_screen_indices, range_indices];
        else
            subset=range_indices(round(numel(range_indices)*(1:threshold)/threshold)); % This way of choosing the indices ensures a uniform distribution of the time points within the range
            post_screen_indices=[post_screen_indices, subset];
        end
    end
    
    %%% done with diluting -- return screened time points
    reduced=time_points(post_screen_indices);
end