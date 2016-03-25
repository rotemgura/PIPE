function results=PIPE_calc_D_of_expanding_gaussians(data,coords,time_points,varargin)
    %%% This function calculates the diffusion coefficient of expanding
    %%% spatial profiles by fitting them to Gaussians and fitting their 
    %%% square widths vs. time to a straight line.
    %%%
    %%% Arguments:
    %%% data -- a 3D matrix of photo-conversion movie
    %%%
    %%% coords - a vector of center of pulse coordinates
    %%% (t0,y0,x0,direction,thickness)
    %%%
    %%% time_points -- a vector of indices of images to analyze, counting
    %%% from t0. Images between time points are averaged over to increase
    %%% SNR.
    %%% 
    %%% Optional arguments:
    %%% to_plot - a boolean denoting whether graphical results of fitting
    %%% procedure would be presented to the user
    %%%
    %%% plot_params - a blolean denoting whether other fitting parameters
    %%% (amplitude, background and others) should also be plotted
    %%%
    %%% name_of_movie -- a string that will be included in the plots
    %%%
    %%% pixel_size -- a length in micrometers that corresponds to the
    %%% side-length of one pixel. By default this field is empty and the
    %%% output is presented in units of pixels^2/frame
    %%%
    %%% imaging_rate -- a rate in frames per second that states the imaging
    %%% rate of the analyzed images
    %%%
    %%% correct_immobile_fraction -- a boolean determines whether PIPE
    %%% would try to correct for an immobile fraction if it finds one
    
    pars=inputParser;
    addParameter(pars,'to_plot',0)
    addParameter(pars,'plot_params',0)
    addParameter(pars,'name_of_movie','')
    addParameter(pars,'pixel_size',0)
    addParameter(pars,'imaging_rate',0)
    addParameter(pars,'correct_immobile_fraction',1)
    parse(pars,varargin{:})
    
    N_POINTS_THRESHOLD=5; %This threshold sets the minimum number of points required for linear regression
    
    %%% calculate the widths by fitting the intensity profiles to Gaussians
    results=PIPE_calc_gaussian_widths(data,coords,time_points,'to_plot',pars.Results.to_plot, ...
        'plot_params',pars.Results.plot_params,'name_of_movie',pars.Results.name_of_movie); %for simplicity the r2 threshold is defined within that function
    if isempty(results)
        return
    end
    if results.immobile_fraction & pars.Results.correct_immobile_fraction
        data=PIPE_correct_immobile_fraction(data,coords,time_points);
        results=PIPE_calc_gaussian_widths(data,coords,time_points,'to_plot',pars.Results.to_plot, ...
            'plot_params',pars.Results.plot_params,'name_of_movie',pars.Results.name_of_movie);
        results.immobile_fraction_correction_attempted=true;
    else
        results.immobile_fraction_correction_attempted=false;
    end
    
    %%% make sure there are enough widths that passed the r2 threshold
    if numel(results.widths)<N_POINTS_THRESHOLD
        'There aren''t enough data points that fit well enough to Gaussians; try to increase the gap between time points or lower the goodness of fit threshold'
        results.D=0;
        results.err_D=0;
        return
    end
    
    %%% fit the widths to a straight line
    x=results.time_points;
    y=results.widths.^2;
    err=2*results.err_widths.*results.widths;
    linfit=fit(x(:),y(:),'poly1','Weights',err(:).^(-2));
    ci=confint(linfit,0.67);
    
    % obtain D and err_D from line slope and confidence intervals
    results.D=linfit.p1/4;
    results.err_D=(ci(2,1)-ci(1,1))/2;
    if pars.Results.pixel_size & pars.Results.imaging_rate
        results.D=results.D*pars.Results.pixel_size^2*pars.Results.imaging_rate;
        results.err_D=results.err_D*pars.Results.pixel_size^2*pars.Results.imaging_rate;
        results.units_of_D='micrometers^2/sec';
    else
        results.units_of_D='pixels^2/frame';
    end
    
    %%% plot the widths and the linear fit
    if pars.Results.to_plot
        figure('units','inches','position',[5 5 6 6])
        errorbar(x,y,err,'o','markersize',16,'linewidth',3)
        hold on
        y2=plot(linfit);
        legend('off')
        set(y2,'linewidth',3,'color',[0.8500 0.3250 0.0980])
        set(gca,'fontsize',20,'xlim',[0, get(gca,'xlim')*[0;1]],'ylim',[0, get(gca,'ylim')*[0;1]])
        if strcmp(results.units_of_D,'pixels^2/frame')
            text_ins=['D=' num2str(round(results.D*100)/100.0) '_{px^2/frame}']; % can use round(results.D,2) instead, but it only works in new Matlab releases
        else
            text_ins=['D=' num2str(round(results.D*100)/100.0) '_{\mum^2/sec}'];
        end
        text(max(x)*.5,max(y)/4,text_ins,'fontsize',20,'fontweight','bold')
        xlabel('Time (# frames)')
        ylabel('Width^2 (# pixels)^2')
        set(gca,'fontsize',30,'linewidth',2,'ticklength',[0.02,0.02])
        set(gcf,'color','white')
        htitle=title(pars.Results.name_of_movie,'fontsize',25);
        set(htitle,'interpreter','none')
    end
end