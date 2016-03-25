function corrected=PIPE_correct_baseline(tiffs,t_,varargin)
    %%% This function calculates the baseline fluorescence from a movie and
    %%% corrects for that baseline assuming it's static and independent on
    %%% laser power (can be subtracted out).
    pars=inputParser;
    addParameter(pars,'smoothing_window',25)
    addParameter(pars,'estimate_from','beginning')
    addParameter(pars,'estimation_window',100)
    parse(pars,varargin{:})
    
    corrected=zeros(size(tiffs));
    
    t_average=pars.Results.estimation_window;
     
    %%% estimate baseline
    if isequal(pars.Results.estimate_from,'beginning')
        baseline=PIPE_smooth_with_gaussian_kernel(mean(tiffs(:,:,max(1,t_-t_average):t_-1),3),pars.Results.smoothing_window);
    elseif isequal(pars.Results.estimate_from,'end')
        baseline=PIPE_smooth_with_gaussian_kernel(mean(tiffs(:,:,end-t_average+1:end),3),pars.Results.smoothing_window);
    end
        
    %%% subtract out baseline
    for i=1:size(tiffs,3)
        corrected(:,:,i)=tiffs(:,:,i)-baseline;
    end
 end