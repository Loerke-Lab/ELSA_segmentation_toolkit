function [] = ELSA_3segmentCellsV5(data,checkInterval,tvec)
% This function segments the images specified by the ImageFileList (a field
% contained in 'data') and the time vector 'tvec' using a seeded watershed 
% technique. Before running this function the seeds must be initialized for 
% the first frame. 
%
% Input:    data: structure that contains lists of image files (as created by the
%                   function ELSA_1loadImageList) and the source (or path) to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%           checkInterval: Number of frames to pass before checking
%                           segmentation and allowing user to modify. If
%                           the movie is prone to frequent errors then use
%                           a smaller checkInterval for more frequent
%                           checking.
%           tvec: time vector specifying the frames to be segmented. 
%                   (e.g. [1:20])
%
% Output:   no function output - results are written into the specified
%           SegmentationData/frame directory
%
% 03/05/2011 Dinah Loerke
% 11/30/2021 Tim Vanderleest (most recent update)


% Record original directory (to return to it at the end)
od = cd;

% 2D Gaussian filter size used for segmentation. Note: I commonly use the 
% value 2 but you can change if it improves segmentation quality.
filtersize = 2;

% Get image file list for the movie from the data structure. Note: if this is a dual
% channel movie then the other list (ImageFileListCh2) may be the
% channel used for segmentation (so change this line if needed).
list = data.ImageFileList;

% Get the total number of frames available in the Movie
Nframes = length(list);


% The time vector can either increase (e.g. 1:1200) or decrease (e.g.
% 1200:-1:1) if segmentation is easier in reverse. So based on the time
% vector we will determine whether the time step is positive or negative
timeStepSign = sign(tvec(2)-tvec(1));


% loop over all timepoints given by tvec
for tidx = 1:length(tvec)
    
    % specify t for time point in the loop
    t = tvec(tidx);
    
    
    % exitflag when true will break from the loop, the user can choose to exit from the
    % modification menu.
    exitflag = false;
    

    % display the frame number being processed
    fprintf('segmentation @ timepoint %04d',t);
    
    % load current image from list and perform 2D Gaussian image filtering.   
    Image = mat2gray(imread(list{t}));
    ImageFiltered = filterImage3DpaddedEdges(Image,'Gauss',filtersize);
    
    % get folder names for both the current and the previous frame. The
    % current frame is the one we are segmenting and the previous one is
    % the one we get seeds from.
    currentfolder = sprintf('frame%04d',t);
        
    % This If block gets the proper path to get seeds
    if t==1 && timeStepSign>0 % Advancing forward at the first frame there are no previous seeds, so use current seeds if they have been initialized.
        previousfolder = currentfolder;
    elseif t==Nframes && timeStepSign<0 % Advancing backwards at the last frame there are no previous seeds, so use current seeds if starting from the end.
        previousfolder = currentfolder;
    else % in all other cases take seeds from the previous time point which depends on which direction we are stepping in time 
        previousfolder = sprintf('frame%04d',t-timeStepSign);
    end
    
    % Move to the SegmentationData directory where all segmentation results are saved
    cd(data.Source);
    cd('SegmentationData');
    
    % Load seeds and mask from previous folder
    cd(previousfolder);
    previousSeeds = load('seeds.mat').seeds;
    previousMask = load('mask.mat').mask;
    
    
    % At frames determined by the checkInterval make any necessary changes to
    % the seeds and mask if there are errors or new cells need to be added
    if ~logical(mod(t,checkInterval)) 
        
        title_string = ['Frame ',num2str(t)];
        [previousSeeds,previousMask,exitflag] = modifySeedsAndMaskV5(ImageFiltered,previousSeeds,previousMask,title_string);
  
    end % Check Interval seed modification loop
    
    % perform watershed segmentation to get segmentation lines
    % (ImageSegment), a label matrix image (ImageBWlabel), updated seeds
    % and updated mask. Seed Type can be set to 'Point' or 'Extended'. 
    [ImageSegment,ImageBWlabel,seeds,mask] = wsSegmentSingleImageV5(ImageFiltered,previousSeeds,previousMask);
    
    
    % Now switch to the folder for the current time frame to save new results to
    cd ..
    cd(currentfolder);
    
    % save results into the current results folder
    save('ImageSegment', 'ImageSegment');
    save('ImageBWlabel', 'ImageBWlabel');
    save('seeds','seeds');
    save('mask','mask');
    
    if exitflag % if exit flag has been switched to true break from function
        break
    end 

    % delete the frame number being processed before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for t

% enter a new line
fprintf('\n');

% return to the original directory
cd(od);
    
end % of function


%% =======================================================================
%
%                               subfunctions
%
%  =======================================================================

function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageV5(image,seeds_in,mask_in)
% Perform seeded watershed segmentation on image given the input seeds and 
% mask logical arrays. This function will output two different types of
% seeds to be saved for segmentation of the next frame
% Inputs:       
%           image: an image that has already been filtered for segmentation
%           seeds_in: logical matrix of seed regions for cells.
%           mask_in: logical matrix of background mask regions.
%
%
% Outputs:       
%           imageSegmented: logical matrix where segmented regions are true and
%           boundary lines (interfaces) are false.
%           imageBWlabel: an integer labeled image where cells have unique
%           integer labels, interfaces are zeros, and background is -1.
%           seeds_out: logical matrix of new seed regions for cells
%           mask_in: logical matrix of new background mask regions.
%


% impose minima for both the watershed seeds and the background mask
imageMasked = imimposemin(image, seeds_in | mask_in);

% perform watershed, set all labeled regions to 1, convert to logical (BW)
imageSegmented = watershed(imageMasked,8);
imageSegmented(imageSegmented>0) = 1;
imageSegmented = logical(imageSegmented);

% convert BW to label image
imageBWlabel = bwlabel(imageSegmented,8);


% Note that any areas that are close to the edge of the image or are part 
% of the background mask are set to a value of -1 to mark them as 
% background - this will ensure in the subsequent analysis that these cells 
% are not used for direct characterization of nodes, since they do not 
% yield full neighborhood information
edgepixels1 = imageBWlabel(:,1:4);
edgepixels2 = imageBWlabel(1:4,:);
edgepixels3 = imageBWlabel(:,(end-3):end);
edgepixels4 = imageBWlabel((end-3):end,:);
maskpixels  = imageBWlabel(mask_in);

%contiguous areas that touch the edge    
edgeOrMaskAreas = unique( [edgepixels1(:);edgepixels2(:);edgepixels3(:);edgepixels4(:);maskpixels(:)]);

% set these areas to value -1
for ie=1:length(edgeOrMaskAreas)
    cvalue = edgeOrMaskAreas(ie);
    if cvalue>0
        imageBWlabel(imageBWlabel==cvalue) = -1;
    end
end

% we get seeds by a "shrink" operation of the regions that are positive
seeds_out = bwmorph(imageBWlabel > 0,'shrink',7);

% we get the mask by negating the dilation of the positive regions
mask_out = ~imdilate(imageBWlabel>0,strel('disk',7));


end % of subfunction



function [seeds_out,mask_out,exitflag] = modifySeedsAndMaskV5(I,seeds_in,mask_in,titlestr)
%MODIFYSEEDSANDMASK allows a user to easily make changes to watershed seeds
%and background mask.
%
% INPUTS:       I:      The image for verifying seeds and masks.
%               seeds_in:  Watershed seeds (2D binary).
%               mask_in:   Mask of background region or regions (2D binary).
%               titlestr: title to be displayed on the figure (optional).
%
% OUTPUTS:      seeds_out: Modified seeds.
%               mask_out:  Modified mask.
%               exitflag: true if user selects the 'Quit' option. exitflag is used to
%               break from the parent function.
%
% 12/18/12 Timothy Vanderleest
%  6/28/22 Added an Undo option which undoes the last change

% set magnification value, exitflag default, and structuring element
magVal = 240;
exitflag = false;
SEdisk10   = strel('disk',10);

% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);

% initialize the modified seeds and mask to the input seeds and mask
seeds_out = seeds_in;
mask_out = mask_in;


% generate segmentation to see segmentation overlay along with seeds and mask
[Ifilt] = filterImage3DpaddedEdges(I, 'Gauss', 2);
imageMasked = imimposemin(Ifilt, seeds_out|mask_out);
imageSegmented = watershed(imageMasked,8);
segboundarylines = imageSegmented<1;

%set mask regions to zero
maskLabels = unique(imageSegmented(mask_out));
for ii=1:length(maskLabels)
    imageSegmented(imageSegmented==maskLabels(ii)) = 0;
end

% initialize RGB image where 2nd and 3rd layers will have seed and mask information
imgrgb = zeros([size(I),3]);

% set unsatisfactory to true
unsatisfactory = true;

% continue to loop while unsatisfactory is true
while unsatisfactory

    % update translucent RGB image with modseeds and modmask and display
    imgrgb(:,:,1) = 0.4*im2double(segboundarylines)+0.6*I;
    imgrgb(:,:,2) = 0.4*im2double(seeds_out)+0.6*I;
    imgrgb(:,:,3) = 0.4*im2double(mask_out)+0.6*I;
    imshow(imgrgb,'InitialMagnification',magVal)
    %figure('units','normalized','outerposition',[0 0 1 1])
%     imshow(imgrgb)
%     set(gcf,'Position',[141 29 1300 776])
    
    % add title string if it was an input
    if nargin > 3
        title(titlestr,'FontSize',16)
    end
    

    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Mask to Seed','Seed to Mask',...
        'Add Line Seed','Add Polygon Seed','Add Polygon Mask','Click to Remove',...
        'Remove Polygon','Remove Spot','re-shrink seeds','Undo','Quit');
    
    switch choice

        % Save and Quit
        case 1
            unsatisfactory = false;
            

        % Mask to Seed
        case 2
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
             
            % this case is designed to make adding seeds super easy by
            % just a click of the mouse
            imshow(imgrgb,'InitialMagnification',magVal); 
            
            % generate fresh seeds and perform watershed
            [i2seeds] = im2seeds(I,9);
            imageMasked = imimposemin(I, i2seeds);
            bwl = watershed(imageMasked,8);
            
            title(gca,'Click on Mask Region to Convert to Seed','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)
                if mask_out(x(k),y(k)) % if mask region is clicked dilate it only
                    newseedlabel = bwl(x(k),y(k));
                    newseedmask  = bwl == newseedlabel;
                    BWdilate = imdilate(newseedmask,strel('disk',7));
                    %BWerode = imerode(newsegmask,strel('disk',4));
                    BWnewseed= bwmorph(newseedmask,'shrink',7);
                    mask_out(BWdilate) = false;
                    seeds_out(BWnewseed) = true;
                end
            end
            

        % convert Seed to Mask    
        case 3
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            % this case is designed to make converting seeds to become
            % part of the background mask
            imshow(imgrgb,'InitialMagnification',magVal); 
            title(gca,'Click on Seed to Convert to Mask','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            
            BW = false(size(I));
            for k=1:length(x)
                if seeds_out(x(k),y(k)) 
                    L = bwlabel(seeds_out);
                    label = L(x(k),y(k));
                    BW(L == label) = true;
                    seeds_out(L==label) = false;
                end
            end
            BW = imclose(BW | mask_out,strel('disk',12));
            mask_out = BW;
            
     
        % Add line seed    
        case 4
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            h = imline;
            BW = h.createMask();
            seeds_out = seeds_out | BW;
            

        % Add Polygon Seed
        case 5
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            BW = roipoly;
            title(gca,'Add Polygon Seed','FontSize',16);                
            seeds_out = seeds_out | BW;
            

        % Add Polygon Mask
        case 6
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            mask_out = mask_out | BW;


        % Remove Region
        case 7 
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            title(gca,'Click on Regions to Remove','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            
            for k=1:length(x)
                if mask_out(x(k),y(k))
                    L = bwlabel(mask_out);
                    label = L(x(k),y(k));
                    idx = L == label;
                    mask_out(idx) = false;
                elseif seeds_out(x(k),y(k))
                    L = bwlabel(seeds_out);
                    label = L(x(k),y(k));
                    idx = L == label;
                    seeds_out(idx) = false;
                end
            end
            

        % Remove Polygon (removes polygon from both seeds and mask)
        case 8
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            seeds_out(BW) = false;
            mask_out(BW) = false;
            

        % Remove Spot Region    
        case 9
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
            
            imshow(imgrgb,'InitialMagnification',magVal); 
            title(gca,'Click Spot to Remove','FontSize',16);
            [y, x] = ginput(1);
            x = round(x); y = round(y);
            BW = false(size(seeds_in));
            BW(x,y) = true;
            BW = imdilate(BW,strel('disk',8));
            seeds_out(BW) = false;
            mask_out(BW) = false;
            

        % Re-shrink seeds
         case 10 
            % last seeds and mask for the UNDO option 
            last_seeds = seeds_in;
            last_mask = mask_in;
             
            seeds_out = bwmorph(imageSegmented>0,'shrink',7);
            

        % Undo last change
        case 11 
            seeds_out = last_seeds;
            mask_out  = last_mask;
            

        % Exit Function
        case 12 
            exitflag = true;
            break         

    end
    
    % after modifications have been made use the new seeds and mask to
    % generate a new segmentation and new "shrunken" seeds
    imageMasked = imimposemin(Ifilt, seeds_out|mask_out);
    imageSegmented = watershed(imageMasked,8);
    segboundarylines = imageSegmented<1;
    maskLabels = unique(imageSegmented(mask_out));
    for ii=1:length(maskLabels)
        imageSegmented(imageSegmented==maskLabels(ii)) = 0;
    end
    mask_out = ~imdilate(imageSegmented>0,SEdisk10);


end

% close everything
close all
end





function  [seeds] = im2seeds(I,filtersize)
%IM2SEEDS produces point-like watershed seeds for a 2D image based on
%regional minima of the filtered image.
%
% Input:    I:              2D image.
%           filtersize:     Gaussian filter size. Must have the function
%                           "filterImage3DpaddedEdges.m"
%
% Output:   seeds = mask of points of minima in the image for seeded
%                   watershed.
%
% 12/18/12 Timothy Vanderleest



% filter image
[Igf] = filterImage3DpaddedEdges(I, 'Gauss', filtersize);

% find the regional minimum in the filtered image (the output are points).
IRM = imregionalmin(Igf);


% slightly dilate these points
SE2 = strel('disk',5);
seeds = imdilate(IRM,SE2);

end