function [cell_mask,is_cell] = PIPE_identify_cell(image)
    %%% This function receives an image of a cell and detect the pixels
    %%% that lie inside the cell. 
    %%% Output:
    %%% cell_mask -- a 2D mask of the cell
    %%% is_cell -- a boolean that notes whether there is an actual cell in
    %%% the image, by testing whether the area outside the cell is smaller
    %%% than 10% of the entire image area.
    
    CELL_EXTERIOR_PIXEL_NUMBER_THRESHOLD=numel(image(:))/10;
    
    %%% detect continuous islands of bright fluorescence
    bwedges = edge(image, 'log');
    bwdilated = imdilate(bwedges,strel('disk',2));
    bwfilled = imfill(bwdilated,'holes');
    bweroded = imerode(bwfilled,strel('diamond',1));
    bwlabeled = bwlabel(bweroded);
    
    %%% choose the largest island
    stats = regionprops(bwlabeled,'Area');
    areas = [stats.Area];
    [~, ind] = max(areas);
    
    %%% create cell mask
    cell_mask = zeros(size(image));
    cell_mask(bwlabeled==ind) = 1;
    if sum(sum(~cell_mask))<CELL_EXTERIOR_PIXEL_NUMBER_THRESHOLD
         is_cell=0;
     else
         is_cell=1;
    end
end