function corrected=PIPE_correct_immobile_fraction(data,coords,time_points)
    %%% This function receives data that has been diagnosed to contain an
    %%% immobile fraction of fluorophores. The function averages and smooths 
    %%% the data between the second-to-last and the last time points, and
    %%% subtracts that from the entire data.
    
    % calculate average of images between last and second-to-last time points
    average_image=mean(data(:,:,coords(1)+time_points(end-1):coords(1)+time_points(end)),3);
    
    % smooth the average image to obtain an average shape of the immobile fraction
    immobile_fraction_snapshot=PIPE_smooth_with_gaussian_kernel(average_image,3);
    
    % subtract the immbile fraction from the data
    corrected=zeros(size(data));
    for i=1:size(data,3)
        corrected(:,:,i)=data(:,:,i)-immobile_fraction_snapshot;
    end
end