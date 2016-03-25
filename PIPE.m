function results=PIPE(varargin)
    %%% This analysis operates on microscopy movies that depict
    %%% photo-conversion experiments. The analysis receives a stack of .tif
    %%% images, tracks the expansion of fluorescence intensity profiles and
    %%% returns a diffusion coefficent describing the observed expansion.
    %%%
    %%% Usage: all arguents are optional, but make sure to supply the raw
    %%% data either from a 3D matrix (argument 'movie') or a folder
    %%% containing .tif images ('path').
    %%%
    %%% When using the 'path' argument, check that 'channel_number' states
    %%% the correct channel of the data. The channel number should appear
    %%% in the .tif file names after the letter C. See further requirements
    %%% on file names in documentation of function PIPE_read_2d_movie.m.
    %%%
    %%% This program outputs both graphical information and a data
    %%% structure that contains the results of the analysis. Information 
    %%% about the diffusion coefficient is found in the fields results.D,
    %%% results.err_D and results.units_of_D. See the user manual for a
    %%% complete documentation of the other output fields.
    %%%
    %%% To present the diffusion coefficient in units of micrometers^2/sec,
    %%% provide the values of the argument 'pixel_size' and 'imaging_rate'
    %%% from the imaging system that was used to record the data. Without
    %%% these values, the default units used are pixels^2/frame.
    %%% 
    %%% Other arguments allow optimization of the PIPE analysis by
    %%% constraining the time points to be analyzed, specifying smoothing
    %%% lengths and determining how to estimate background intensity. See
    %%% the user manual for a complete documentation of the arguments.
    %%%
    %%% Examples:
    %%% read data from an image folder from channel 0:
    %%% results=PIPE('path','myfolder/movie1/','channel_number',0)
    %%%
    %%% analyze an existing 3D data matrix (y,x,t) without plotting results
    %%% results=PIPE('movie','movie1','to_plot',0)
    %%%
    %%% Copyright: Created by Rotem Gura Sadovsky as part of work at MIT, 
    %%% March 2016. Read license information for public use of this 
    %%% software and citation guidelines.
    

    pars=inputParser;
    addParameter(pars,'movie','')
    addParameter(pars,'path','')
    addParameter(pars,'channel_number',1)
    addParameter(pars,'to_plot',1)
    addParameter(pars,'to_plot_extended',0)
    addParameter(pars,'name','')
    addParameter(pars,'direction',0)
    addParameter(pars,'thickness_of_profile',30)
    addParameter(pars,'time_points','')
    addParameter(pars,'correct_baseline',1)
    addParameter(pars,'exclude_extra_cellular',1)
    addParameter(pars,'max_time_point',0)
    addParameter(pars,'pulse_coords','')
    addParameter(pars,'max_n_data_points',30)
    addParameter(pars,'estimate_baseline_from','beginning')
    addParameter(pars,'baseline_estimation_window',100)
    addParameter(pars,'smoothing_window',25)
    addParameter(pars,'pixel_size',0) %in micrometers
    addParameter(pars,'imaging_rate',0) %in frames per second
    addParameter(pars,'correct_immobile_fraction',1)
    parse(pars,varargin{:})
    
    %%% read data
    name=pars.Results.name;
    if ~isempty(pars.Results.movie) % data was already read
        data=pars.Results.movie;
    else
        if ~isempty(pars.Results.path) % read data from given path
            data=PIPE_read_2d_movie(pars.Results.path,'beginning','end',pars.Results.channel_number);
            if data==0
                ['Could not read movie in ' pars.Results.path]
                results='';
                return
            end
            if isempty(name) % give the movie a name from the file path
                if ismember('/',pars.Results.path)
                    temp=strsplit(pars.Results.path,'/');
                elseif ismember('\',pars.Results.path)
                    temp=strsplit(pars.Results.path,'\');
                else
                    temp={pars.Results.path};
                end
                if isempty(temp{end})
                    name=temp{end-1};
                else
                    name=temp{end};
                end
            end
        else
            error('No data source found. Please input either raw data or .tif file folder path')
        end
    end
    
    % transform zero-intensity pixels to a small finite value to avoid 
    % confusion with excluded data, which will contain zeros
    data(~logical(data))=max(data(:))/1e5;
    
    %%% find coordinates of photo-conversion pulse
    if pars.Results.pulse_coords
        t0=pars.Results.pulse_coords(1);
        t_=pars.Results.pulse_coords(2);
        y0=pars.Results.pulse_coords(3);
        x0=pars.Results.pulse_coords(4);
    else
        [temp1, temp2, temp3]=PIPE_find_pulse_center(data,0,pars.Results.exclude_extra_cellular);
        [t0,t_,y0,x0]=PIPE_detect_pulse_end(data,temp1,temp2,temp3);
        if t0==0
            ['Could not detect the center of pulse in' pars.Results.path]
            results='';
            return
        end
    end
    coords=[t0,y0(end),x0(end),pars.Results.direction,pars.Results.thickness_of_profile];
    
    
    %%% pre-process the data
    % segment the images and record distance from edge of field of view to cell boundary
    min_distance_to_boundary=''; 
    if pars.Results.exclude_extra_cellular
        [mask,is_cell]=PIPE_identify_cell(mean(data(:,:,end-30:end),3));
        if is_cell %if cell boundary is visible then calculate distance to boundary
            min_distance_to_boundary=PIPE_find_min_distance_to_boundary(mask,coords);
        end
    end
    
    % correct baseline
    if pars.Results.correct_baseline
        data=PIPE_correct_baseline(data,t_,'estimate_from',pars.Results.estimate_baseline_from,'estimation_window',pars.Results.baseline_estimation_window,'smoothing_window',pars.Results.smoothing_window);
    end
    
    % exclude data outside cell membrane or other boundaries
    if pars.Results.exclude_extra_cellular
        for i=1:size(data,3)
            data(:,:,i)=data(:,:,i).*mask;
        end
    end
    
    % check that the cell is wide enough compared with the initial width of intensity profile
    width0=PIPE_estimate_width0(data,coords);
    if width0==0
        ['Could not estimate width of initial intensity profile']
    elseif ~isempty(min_distance_to_boundary) && width0/double(min_distance_to_boundary)>1/3
        ['Warning: PIPE detected a boundary close to the point of photo-conversion. The signal might ' ...
            'not have enough room to expand, which might result in unreliable output.']
    end
    
    %%% construct averaging scheme: find the time points that optimally extract information from data
    if ~isempty(pars.Results.time_points) %averaging scheme was specified by user
        time_points=pars.Results.time_points;
        if numel(time_points)>pars.Results.max_n_data_points
            'Too many time points have been entered manually. Enter fewer time points or change the threshold'
            results='';
            return
        end
    else %call averaging scheme constructor
        time_points=PIPE_calc_time_points_by_noise_analysis(data,coords,t_,'max_time_point',pars.Results.max_time_point);
        if time_points==0
            ['Could not obtain high enough SNR from movie ' name '. Try changing the threshold or choose the time points manually.']
            results='';
            return
        end
    end
    
    % reduce the number of time points if total number is capped
    if numel(time_points)>pars.Results.max_n_data_points+1
        time_points=PIPE_reduce_n_time_points(time_points,pars.Results.max_n_data_points);
    end
    
    %%% calculate diffusion coefficients
    results=PIPE_calc_D_of_expanding_gaussians(data,coords,time_points,'to_plot',pars.Results.to_plot, ...
        'plot_params',pars.Results.to_plot_extended,'name_of_movie',name,'pixel_size',pars.Results.pixel_size, ...
        'imaging_rate',pars.Results.imaging_rate,'correct_immobile_fraction',pars.Results.correct_immobile_fraction);
end