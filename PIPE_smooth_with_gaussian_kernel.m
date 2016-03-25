function smoothed=PIPE_smooth_with_gaussian_kernel(data,sigma,varargin)
    %%% This function receives a matrix 'data' and a length scale 'sigma' and
    %%% smooths data by filtering it with a Gaussian kernel of width sigma
    
    pars=inputParser;
    addParameter(pars,'smoothing_direction','both')
    parse(pars,varargin{:})
    
    if size(data,1)==1 || size(data,2)==1 
        
        %%%data is 1D
        g=@(x0,sigma,x)exp(-(x-x0).^2/(2*sigma^2)); %define 1D gaussian
        range=-(numel(data)-1):numel(data)-1;
        if size(data,2)==1
            range=range'; %make range a column vector like data
        end
        gauss=g(0,sigma,range); %define the 1D smoothing kernel
        
    else
        
        %%% data is 2D
        if isequal(pars.Results.smoothing_direction,'both') %smooth data in both directions
            g=@(x0,y0,sigma,x,y)exp(-((x-x0).^2+(y-y0).^2)/(2*sigma^2)); %define 2D gaussian
            range_x=-(size(data,2)-1):size(data,2)-1;
            range_y=(-(size(data,1)-1):size(data,1)-1)';
            gauss=g(0,0,sigma,repmat(range_x,[2*size(data,1)-1,1]),repmat(range_y,[1,2*size(data,2)-1])); %define the 2D smoothing kernel
        
        else
            g=@(x0,sigma,x)exp(-(x-x0).^2/(2*sigma^2));
            if pars.Results.smoothing_direction==1 % smooth 2D data in the y direction
                range=(-(size(data,1)-1):size(data,1)-1)';
            else % smooth 2D data in the x direction
                range=-(size(data,2)-1):size(data,2)-1;
            end
            gauss=g(0,sigma,range); %define the 1D smoothing kernel for the 2D data
        end
    end
    
    %%%convolute the data with the Gaussian and normalize accordingly
    smoothed=conv2(data,gauss,'same')./conv2(ones(size(data)),gauss,'same');
end