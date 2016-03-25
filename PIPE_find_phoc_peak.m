function T=PIPE_find_phoc_peak(data,varargin)
    %%% This function looks for a peak in a noisy time series.
    %%% The function clusters the points that pass a threshold and looks for
    %%% a big enough stretch (default size spans 10 consecutive time points)
    %%% The function returns the index of the peak
    
    pars=inputParser;
    addParameter(pars,'start_index',1)
    addParameter(pars,'end_index',numel(data))
    parse(pars,varargin{:});
    
    %%% subset data and make sure it's a column vector
    data=data(pars.Results.start_index:pars.Results.end_index);
    data=data(:);
    
    %%% use thresholding to find stretches of high fluorescence signal
    range=1:-0.01:0;
    for thresh=range
        vec=data>max(data)*thresh;
        labeled=bwlabel(imerode(imdilate(vec',ones(1,3)),ones(1,3)));
        stats=regionprops(labeled,'Area');
        [max_area,max_label]=max([stats.Area]); %find longest stretch
        if max_area>10 %check that stretch is long enough
            indices=find(labeled==max_label); %choose longest stretch
            break
        end
    end
    if thresh==0 %failed to find peak
        T=-1;
        'Could not detect the end of pulse using ''cluster peaks'' algorithm'
        return
    end
    
    %%% find the peak location
    if indices(end)==max_area % this means that the wanted peak is in the beginning of data so it is probably steep.
        [~,T]=max(data'.*(labeled==max_label)); %time of peak is just the time of maximal signal
    else % this means that the wanted peak is in the middle of data and might be noisy.
        T=round(mean(indices)); %assume peak is symmetric and thus time of peak is in the middle of the high intensity stretch
    end
end
