function movie=PIPE_read_2d_movie(folder_path, start, finish, channel_ind, varargin)
    %%% This function reads a folder of .tif images and holds the data in a
    %%% 3D matrix. The .tif file names must be written in a certain
    %%% structure, usually <string>T#C#.tif. The T index denotes the
    %%% time of the image in the sequence; the files must be identified 
    %%% with a unique T index, without gaps in the index sequence. The C 
    %%% index denotes the channel number in case the movie is separated 
    %%% into distinct channels (e.g. green and red fluorescence). If the 
    %%% movie contains only one channel, the file names can lack the 'C#', 
    %%% and the argument channel_ind has to be empty ('').
    %%% 
    %%% Arguments:
    %%% folder_path -- the folder containing the .tif images
    %%%
    %%% start -- the T index from which to start reading. It can be either
    %%% a number or the word 'beginning', which will start reading from the
    %%% earliest T index.
    %%% 
    %%% finish -- the latest T index to read. It can be either a number or
    %%% the word 'end', which will end reading at the latest T index.
    %%%
    %%% channel_ind -- the C index appearing in the file names to be read.
    %%%
    %%% Optional arguments:
    %%% min_movie_length -- the minimum number of images in a movie. This
    %%% argument is 30 by default.
    %%%
    %%% Example: to read the whole red channel from a green & red channel
    %%% movie, run
    %%% mymov=read_2d_movie('/myfolder/movie/','beginning','end',1);
    
    pars=inputParser;
    addParameter(pars,'min_movie_length',30)
    parse(pars,varargin{:});
    
    %%% check arguments
    % add separator to end of path if missing
    if ~strcmp(folder_path(end),'/') && ~strcmp(folder_path(end),'\')
        if ismember('\',folder_path)
            folder_path(end+1)='\';
        else
            folder_path(end+1)='/';
        end
    end
    
    % check for sufficient files in given folder
    d = dir(folder_path);
    if length(d)<3
        'Could not find image files in folder'
        movie=0;
        return
    end
    if length(d)<15
        'Could not find enough image files for a movie in folder'
        movie=0;
        return
    end
    all_filenames = {d.name};
    
    %%% filter files and prepare for reading
    % filter files by channel
    if ~isempty(channel_ind)
        hits=cellfun(@(x)~isempty(regexp(x,['[Cc]' num2str(channel_ind) '.tif'])),all_filenames);
    else
        hits=cellfun(@(x)~isempty(regexp(x,['.tif'])),all_filenames);
    end
    filenames=all_filenames(hits);
    if length(filenames)==0
        'Could not find image files of the specified channel number'
        movie=0;
        return
    end
    
    % find the first time point in list
    first_time_point=Inf;
    for i=1:numel(filenames)
        time_point=regexp(filenames{i},'[Tt]\d+','match'); %look for the T index in file name
        index=str2num(time_point{end}(2:end)); %take last match, in case file name contains the pattern
        if index<first_time_point
            first_time_point=index;
        end
    end 
    
    % sort list according to time points
    sorted_filenames=cell(numel(filenames),1);
    for i=1:numel(filenames)
        time_point=regexp(filenames{i},'[Tt]\d+','match'); %look for the T index in file name
        index=str2num(time_point{end}(2:end)); %take last match, in case file name contains the pattern
        sorted_filenames(index-first_time_point+1) = filenames(i);
    end
    
    % read the whole image if arguments instruct so
    if strcmp(start,'beginning') %put the first frame number into 'start'
        start = first_time_point;
    end
    if strcmp(finish,'end') %put the last frame number into 'finish'
        finish = numel(sorted_filenames)+first_time_point-1;
    end
    
    % check that movie is sufficiently long
    if finish-start<pars.Results.min_movie_length
        ['Movie is shorter than the minimal defined length (' num2str(pars.Results.min_movie_length) ' frames)']
        movie=0;
        return
    end
    
    %%% read images
    for i=start-first_time_point+1:finish-first_time_point+1
        im=imread([folder_path sorted_filenames{i}]);
        if i==start-first_time_point+1 %create data matrix
            [ymax,xmax]=size(im); % y is the first index in matlab matrices
            movie=zeros(ymax, xmax, finish - start + 1);
        end
        movie(:,:,i - (start - first_time_point)) = im;
    end
end