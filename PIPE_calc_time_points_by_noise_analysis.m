function time_points=PIPE_calc_time_points_by_noise_analysis(data,coords,t_,varargin)
    %%% This function analyzes the noise in the data to find a set of time
    %%% points that allow sufficient averaging of the data to potentially
    %%% achieve good Gaussian fits. The function calculates how noisy the
    %%% intensity profile is at each time point and uses a predetermined
    %%% function to estimate how much averaging is needed to push the noise
    %%% below a certain threshold.
    %%%
    %%% Arguments:
    %%% data -- a 3D matrix
    %%%
    %%% coords -- a vector of center-of-pulse coordinates
    %%% (t0,y0,x0,direction,thickness)
    %%%
    %%% t_ -- time of beginning of photo-conversion pulse
    %%%
    %%% Optional arguments:
    %%%
    %%% max_time_point -- maximal time point to analyze, counting from t0
    %%%
    %%% r2_threshold -- threshold for acceptance of Gaussian fit
    
    MAX_NOISE=300; % corresponds to a low SNR of 1/300
    
    pars=inputParser;
    addParameter(pars,'max_time_point',0)
    addParameter(pars,'r2_threshold',0.65)
    parse(pars,varargin{:})

    %%% set parameters for analysis
    % unpack coordinates
    temp=num2cell(coords);
    [t0,y0,x0,direction,thickness]=temp{:};
    
    % determine max time point to analyze
    if pars.Results.max_time_point
        max_time=min(t0+pars.Results.max_time_point,size(data,3));
    else
        max_time=size(data,3);
    end
    
    %%% estimate the background from the pre-PhoC data
    prephoc_data=data(:,:,1:t_-1);
    background=sum(prephoc_data(:))/sum(logical(prephoc_data(:)));
    
    %%% get the intensity profiles
    % define profile mask
    [profile_mask,perp_direction]=PIPE_define_1d_profile_mask(y0,x0,direction,thickness,size(data(:,:,1)));

    % get first profile to determine its length
    profile1=PIPE_apply_profile_mask(data(:,:,t0),profile_mask, perp_direction, thickness);
    if isempty(profile1)
        'Could not obtain an intensity profile. Check image segmentation.'
        time_points=0;
        return
    end
    profiles=zeros(max_time-t0+1,numel(profile1));

    % get the profiles from the rest of the images
    for i=t0:max_time
        profiles(i-t0+1,:)=PIPE_apply_profile_mask(data(:,:,i),profile_mask, perp_direction, thickness);
    end
    
    
    %%% calculate noise (inverse SNR) for all time points
    % find s.d. at the first 30 data points of the profile, where the shape should be mostly flat
    stddev=std(profiles(:,1:30),0,2); % this is a vector of fluctuations vs time
    
    % define the signal as the Gaussian amplitude
    peak=max(PIPE_smooth_with_gaussian_kernel(profiles,15,'smoothing_direction',2),[],2);
    bottom=mean(profiles(:,1:30),2);
    signal=PIPE_smooth_with_gaussian_kernel(peak-bottom,3); % this is a smoothed vector of amplitude vs time
    
    % calculate the noise in percentage
    noise_vs_time=PIPE_smooth_with_gaussian_kernel(stddev./signal,3)*100;

    %%% construct a mapping from noise to number of time points to average over 
    % The following formula was empirically obtained from analysis of simulated data:
    % 1) random walk was simulated with 1e5 particles in 2D .
    % 2) the data was binned to pixels and Poissonian noise was added to each pixel.
    % 3) consecutive images were averaged together pixel by pixel, which decreased the noise.
    % 4) the following formula describes how much averaging is needed to decrease the noise from its initial value to 15%.
    n_points=@(noise)noise^1.9/exp(4.99); 

    %%% test the estimate from formula above and increase averaging if needed 
    time_points=[0];
    cur_time=ceil(n_points(noise_vs_time(1)));
    while cur_time<size(profiles,1) && noise_vs_time(cur_time)<MAX_NOISE % stop adding time points at a maximal cutoff on noise
        
        % average over profiles from consecutive time points
        profile=mean(profiles(time_points(end)+1:cur_time,:),1);
        
        % estimate starting values for fitting parameters
        [a0,b0]=max(profile);
        c0=numel(profile)/4;
        
        % fit averaged profile to a Gaussian
        [~,params]=fit((1:numel(profile))',profile',@(a,b,c,x)a*exp(-((x-b)/c).^2)+background,'Startpoint',[a0 b0 c0],'Lower',[0 1 0]);
        
        % check whether R^2 is sufficient
        if params.rsquare>pars.Results.r2_threshold 
            
            % if R^2 was sufficient, make sure it is not sensitive to fluctuations in background 
            [~,params2]=fit([1:numel(profile)]',profile',@(a,b,c,x)a*exp(-((x-b)/c).^2)+background*0.9,'Startpoint',[a0 b0 c0],'Lower',[0 1 0]);
            [~,params3]=fit([1:numel(profile)]',profile',@(a,b,c,x)a*exp(-((x-b)/c).^2)+background*1.1,'Startpoint',[a0 b0 c0],'Lower',[0 1 0]);
            if params2.rsquare>pars.Results.r2_threshold && params3.rsquare>pars.Results.r2_threshold
                
                % the average profile fits well to a Gaussian. Record the
                % number of profiles averaged together by adding to the list
                % of time points the index of the last profile that
                % participated in the current averaging.
                time_points(end+1)=cur_time;
                
                % move on to the next averaging procedure
                cur_time=cur_time+ceil(n_points(noise_vs_time(cur_time)));
            else
                cur_time=cur_time+1;
            end
        else
            cur_time=cur_time+1;
        end        
    end
end