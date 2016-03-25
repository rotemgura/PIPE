function results=PIPE_calc_gaussian_widths(data,coords,time_points,varargin)
%%% This function calculates the widths of spatial profiles in given time
%%% points by fitting the profiles to Gaussians.
    pars=inputParser;
    addParameter(pars,'to_plot',0)
    addParameter(pars,'same_plot',0)
    addParameter(pars,'plot_params',0)
    addParameter(pars,'plot_residuals',0)
    addParameter(pars,'name_of_movie','')
    parse(pars,varargin{:})
    
    %%% initialize variables
    R2_THRESHOLD=0.65;
    MODEL_RESOLUTION=0.1;
    [amp,err_amp,mean_loc,err_mean_loc,width,err_width,r2,sum_r,background,err_background]=deal(zeros(1,numel(time_points)));
    temp=num2cell(coords);
    [t0,y0,x0,direction,thickness]=temp{:};
    movie_duration=size(data,3);
    
    %%% calculate how many frames to average over based on the given time_points vector
    average_over=diff(time_points)-1;
    time_to_end_of_data=movie_duration-time_points(end)-t0;
    average_over(end+1)=min(average_over(end),time_to_end_of_data);
    post_averaging_time_points=time_points+average_over/2;
    
    %%% define the spatial profiles
    [profile_mask,perp_direction]=PIPE_define_1d_profile_mask(y0,x0,direction,thickness,size(data(:,:,1)));
      
    %%% handle pre-plotting logistics
    if pars.Results.to_plot
        h1=figure();
        if pars.Results.plot_residuals
            h2=figure();
        end
        figure(h1)
        num_plot=1;
        total_num_plots = length(time_points);
        subplot_cols = 4;
        subplot_rows = ceil(total_num_plots/4.0);
    end
    
    %%% start operating on each time point
    for i=1:numel(time_points)
        
        %%% get the intensity profile for the current time point
        %average over the images between the current and the next time point
        avg=mean(data(:,:,t0+time_points(i):t0+time_points(i)+average_over(i)),3);
        
        % get the intensity profile from the average image
        profile=PIPE_apply_profile_mask(avg,profile_mask, perp_direction, thickness);
        if isempty(profile)
            results='';
            return
        end
        
        %%% fit the intensity profile to a Gaussian
        fit_from=1;
        fit_to=numel(profile);
        
        if i==1 % create plot matrix
            plot_mat_x=repmat(fit_from:fit_to,[numel(time_points),1]);
            plot_mat_y=repmat(post_averaging_time_points(:),[1,fit_to-fit_from+1]);
            plot_mat_z=zeros(numel(time_points),fit_to-fit_from+1);
            plot_mat_model=zeros(numel(time_points),fit_to-fit_from+1);
        end
        
        % estimate starting points for the Gaussian parameters
        [a0,b0]=max(profile);
        c0=numel(profile)/4;
        d0=0;
        
        % define the fitting function
        % notice that the width c is defined without a factor of 2, and
        % that a background parameter d is also fitted
        gauss=@(a,b,c,d,x)a*exp(-((x-b)/c).^2)+d;
        
        % fit
        [model,params]=fit((fit_from:fit_to)',profile(fit_from:fit_to),gauss,'Startpoint',[a0 b0 c0 d0],'Lower',[0 fit_from 0 -abs(min(profile))],'Upper',[100*(max(profile)-min(profile)) 100*numel(profile) 100*numel(profile) abs(max(profile))]);
        
        % get error from confidence intervals
        ci=confint(model,.67);
        
        
        amp(i)=model.a;
        err_amp(i)=(ci(2,1)-ci(1,1))/2;
        mean_loc(i)=model.b;
        err_mean_loc(i)=(ci(2,2)-ci(1,2))/2;
        width(i)=model.c;
        err_width(i)=(ci(2,3)-ci(1,3))/2;
        r2(i)=params.rsquare;
        background(i)=model.d;
        err_background(i)=(ci(2,4)-ci(1,4))/2;
        residuals=profile(fit_from:fit_to)-model(fit_from:fit_to);
        sum_r(i)=sum(residuals);
        
        %%% plot data and model
        smooth=background(i)+model.a*exp(-(([fit_from:MODEL_RESOLUTION:fit_to]-model.b)./model.c).^2);
        if pars.Results.to_plot
            if pars.Results.same_plot
                plot_mat_z(i,:)=profile;
                plot_mat_model(i,:)=smooth(1:1/MODEL_RESOLUTION:end);
            else
                subplot(subplot_rows,subplot_cols,num_plot)
                plot(profile)
                axis([1,length(profile),min(profile),max(profile)+0.0001]);
                hold on
                plot(fit_from:fit_to,smooth(1:1/MODEL_RESOLUTION:end))
                title(['t=' num2str(post_averaging_time_points(i))],'fontsize',10)
            end
            %plot residuals
            if pars.Results.plot_residuals
                %set(groot,'CurrentFigure',h2)
                figure(h2)
                subplot(subplot_rows,subplot_cols,num_plot)
                plot(fit_from:fit_to,smooth_with_gaussian_kernel(residuals,15),'LineWidth',1.5)
                xlim([fit_from,fit_to])
                %set(groot,'CurrentFigure',h1)
                figure(h1)
            end
            num_plot=num_plot+1;
        end
    end
    if pars.Results.to_plot & pars.Results.same_plot
        fh1=figure(h1);
        y1=waterfall(plot_mat_x,plot_mat_y,plot_mat_z);
        set(y1,'linewidth',2)
        hold on
        y1=waterfall(plot_mat_x,plot_mat_y,plot_mat_model);
        set(y1,'linewidth',2)
        xlabel('Location (px)')
%         hax=get(gca,'xlabel');
%         get(hax,'position')
%         set(hax,'rotation',18,'position',[550 23 0])
        ylabel('Time (frames)')
%         hax=get(gca,'ylabel');
%         set(hax,'rotation',-20,'position',[750 16 0])
%         set(gca,'xtick',[100 300 500])
        zlabel('Intensity (A.U.)')
        set(gca,'fontsize',30,'fontname','helvetica','linewidth',3,'ticklength',[0.02,0.02],'xlim',[fit_from,fit_to],'zlim',[min(plot_mat_z(:)) get(gca,'zlim')*[0;1]])
        set(gcf,'color','white')
        set(fh1,'units','inches','Position',[1 1 12 9])
        view(135,15)
        htitle=title(pars.Results.name_of_movie);
        set(htitle,'interpreter','none')
    end
        
    if pars.Results.plot_params
        figure()
        subplot(2,2,1)
        if post_averaging_time_points(1)==0
            errorbar(post_averaging_time_points(2:end),amp(2:end),err_amp(2:end));
        else
            errorbar(post_averaging_time_points,amp,err_amp);
        end
        set(gca,'yscale','log','xscale','linear','fontsize',20,'xlim',[0,get(gca,'xlim')*[0;1]])
        title('Amplitude(t) (A.U.)','fontsize',20)
        subplot(2,2,2)
        errorbar(post_averaging_time_points,mean_loc,err_mean_loc);
        set(gca,'fontsize',20,'xlim',[0,get(gca,'xlim')*[0;1]])
        title('Mean(t) (# pixels)','fontsize',20)
        subplot(2,2,3)
        errorbar(post_averaging_time_points,background,err_background);
        set(gca,'fontsize',20,'xlim',[0,get(gca,'xlim')*[0;1]])
        title('Background(t) (A.U.)','fontsize',20)
        subplot(2,2,4)
        plot(post_averaging_time_points,r2);
        set(gca,'fontsize',20,'xlim',[0,get(gca,'xlim')*[0;1]])
        title('R^2(t)','fontsize',20)
        h1=suptitle('Fitting parameters vs. time');
        set(h1,'fontsize',25)
        sup=suptitle(pars.Results.name_of_movie);
        set(sup,'fontsize',25,'interpreter','none')
        set(gcf,'color','white')
    end
    %%% return fits that pass the r2 threshold
    indices=r2>R2_THRESHOLD;
    good_widths=width(indices);
    err_good_widths=err_width(indices);
    post_averaging_good_time_points=post_averaging_time_points(indices);
    if is_immobile_fraction(post_averaging_time_points,width,err_width)
        immobile_fraction=true;
    else
        immobile_fraction=false;
    end
    
    %%% pack info into data structure
    results=struct('widths',good_widths,'err_widths',err_good_widths, ...
    'time_points',post_averaging_good_time_points, 'mean_loc', mean_loc, ...
    'sum_r',sum_r,'background',background,'err_background',err_background, ...
    'immobile_fraction',immobile_fraction);
end

function is_immobile=is_immobile_fraction(post_averaging_time_points,width,err_width)
%%% This function receives the parameters of the Gaussian fits and
%%% diagnoses whether there is an immobile fraction of fluorophores in the
%%% data.
    if post_averaging_time_points(end)-post_averaging_time_points(1)<20 || numel(post_averaging_time_points)<7
        is_immobile=false;
        ['cannot diagnose immobile fraction due to small number or range of time points']
        return
    end
    
    
    %%% check condition 1: a significant decrease in Gaussian width toward the end of the movie
    if mean(width(end-2:end))-width(1)<median(width-width(1))/10 && mean(err_width(end-2:end))<mean(width(end-2:end))
        is_immobile=true;
    %%% check condition 2: strong saturation of the width^2
    else
        temp=polyfit(post_averaging_time_points(1:5), width(1:5).^2,1);
        beginning_slope=temp(1);
        temp=polyfit(post_averaging_time_points(end-4:end), width(end-4:end).^2,1);
        end_slope=temp(1);
        if end_slope<beginning_slope/10 && mean(err_width(end-4:end))<mean(width(end-4:end))
            is_immobile=true;
        else
            is_immobile=false;
        end
    end
end
