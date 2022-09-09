function [ ] = ELSA_segmentCellsInTime_LS(data,MovieNum,checkInterval,tvec,zlayer)
% ELSA_segmentCellsInTime_LS segments the images over time using a seeded 
% watershed technique. It also tracks them during the time vector inputted.
%
% INPUT:    data: structure that contains lists of image files (as created by the
%                   function ELSA_1loadImageList) and the source (or path) to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%           MovieNum: the movie number to be analyzed in the case that
%                     multiple movies are stored in the data structure
%           checkInterval: Number of frames to pass before checking
%                           segmentation and allowing user to modify. If
%                           the movie is prone to frequent errors then use
%                           a smaller checkInterval for more frequent
%                           checking.
%           tvec: time vector specifying the frames to be segmented. 
%                   (e.g. [1:20])
%           zlayer: the z layer (depth) of the movie to be segmented. 
%
% OUTPUT:   no function output - results are written into the specified
%           SegmentationData/frame directory


% record original directory (to return to it at the end)
od = cd;

%filtering parameter for watershed segmentation
filtersize = [2 2 1];

% Structuring element for the erosion of seeds and mask, 12 worked well for
% faster movies (up to 7 fps), switched to 18 for Zipper
SEerodeMask = strel('disk',10);
% BWmorph shrink magnitude, 10 worked for faster movies, switched to 15 for
% zipper movies
shrinkVal = 10;

% get list of images
list = data(MovieNum).ImageFileList;

% number of frames in movie and number of layers
Nframes = numel(list);

% DEFAULT time and section numbers are set to min and max available,
% unless different input is specified
tstart = 1;
tend = Nframes;
tlvec = tstart:1:tend;
if nargin>3
    if ~isempty(tvec)
        tlvec = tvec;
        tstart = tvec(1);
    end
    if tlvec(2)>tlvec(1)
        tstep = 1;
    else
        tstep =-1;
    end
end

% set exit flag to false
eflag2 = false;

% loop over all timepoints
for ti = 1:length(tlvec)
    
    % specify t for time point in the loop
    t = tlvec(ti);
    
    % display current progress of processing
    fprintf('segmentation @ timepoint %04d',t);
    
    
    % change directory to frame specific segmentation folder
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    if (t==tstart) && (tstep ==1)
        lframefoldername = sprintf('frame%04d',t);
    elseif (t==Nframes) && (tstep==-1)
        lframefoldername = sprintf('frame%04d',t);
    else
        lframefoldername = sprintf('frame%04d',t-tstep);
    end
    cd(lframefoldername);

    % load last seeds (or current if starting point)
    seeds = load('seeds.mat').seeds;
    lseeds_z = seeds(:,:,zlayer);
    mask = load('mask.mat').mask;
    lmask_z = mask(:,:,zlayer);

    % load ImageBWlabel
    %ImageSegment = load('ImageSegment.mat').ImageSegment;
    ImageBWlabel = load('ImageBWlabel.mat').ImageBWlabel;
    ImageBWlabelOld = ImageBWlabel(:,:,zlayer);
    
    % throw error if there are nans in ImageBWlabel for the zlayer
    if any(isnan(ImageBWlabelOld))
        disp('error')
    end

    % change directory up two levels
    cd ../..
    % get current image
%     I = tiffreadVolume(list{t});
    for z = 1:size(list,2)
        I(:,:,z) = imread(list{t,z});
    end
    Background = I(:,:,zlayer)==0;
    Background = bwareaopen(Background,10);

    % filter image
    Ifilt = imgaussfilt3(I,filtersize);
    
    % fix border seeds and make modifications if necessary.
    if ~logical(mod(t,checkInterval))
        titlestr = strcat('Layer ',num2str(zlayer),', Frame ',num2str(t));
        [lseeds_z,lmask_z,eflag,eflag2] = modifySeedsAndMask_LS(I(:,:,zlayer),Ifilt(:,:,zlayer),lseeds_z,lmask_z,titlestr,SEerodeMask,shrinkVal);
        
        if eflag
            break
        end
  
    end 
    

    % move to current timepoint subfolder...
    cd('SegmentationData');
    cframefoldername = sprintf('frame%04d',t);   
    cd(cframefoldername);
    

    % perform watershed segmentation, which also automatically yields a
    % label matrix and updated seeds and mask for current time point
    [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),lseeds_z,lmask_z,Background,SEerodeMask,shrinkVal);
    
    % load ImageSegment, ImageBWlabel, seeds, and mask
    ImageSegment = load('ImageSegment.mat').ImageSegment;
    ImageBWlabel = load('ImageBWlabel.mat').ImageBWlabel;
    seeds = load('seeds.mat').seeds;
    mask = load('mask.mat').mask;
    
    % perform tracking for ImageBWlabel
    [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabelOld);
    
    % select the specified z layer for ImageSegment, ImageBWlabel, seeds and mask
    ImageSegment(:,:,zlayer) = ImageSegment_z;
    ImageBWlabel(:,:,zlayer) = ImageBWlabel_z;
    seeds(:,:,zlayer) = seeds_z;
    mask(:,:,zlayer) = mask_z;
    
    % save results into the current results folder
    save('ImageSegment', 'ImageSegment');
    save('ImageBWlabel', 'ImageBWlabel');
    save('seeds','seeds');
    save('mask','mask');
  
    % delete the frame number being processed before the next loop   
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
    % if exit flag has been switched to true break from function
    if eflag2
        break
    end
    
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

function [ImageBWlabelNew_Track] = labelTracker(ImageBWlabelNew,ImageBWlabelOld)
% labelTracker tracks the cell labels over time or space in ImageBWlabel
%
% INPUT:    ImageBWlabelNew: an integer labeled image where cells have unique
%                            integer labels, interfaces are zeros, and background is -1. 
%                            The new ImageBWlabel comes after segmentation
%                            is performed
%           ImageBWlabelOld: an integer labeled image where cells have unique
%                            integer labels, interfaces are zeros, and background is -1.
%                            The original ImageBWlabel before segmentation
%                            is performed
%
% OUTPUT:   ImageBWlabelNew_Track: an integer labeled image where cells have unique
%                                  integer labels, interfaces are zeros, and background is -1.
%                                  ImageBWlabel tracked in time or space

% find the maximum value of the old image label
maxOldLabel = max(unique(ImageBWlabelOld(:)));

% make output equal to the new image label input
ImageBWlabelNew_Track = ImageBWlabelNew;

% find unique numbers in ImageBWlabelNew
ulabelsNew = unique(ImageBWlabelNew(:));

% set any numbers less than 1 to empty
ulabelsNew(ulabelsNew<1) = [];

% find the number of unique numbers
NumNew = numel(ulabelsNew);

% loop over the unique numbers
for ii=1:NumNew

    % create a cell mask where ImageBWlabelNew is equal to the unique
    % number in the loop
    maskCell = ImageBWlabelNew == ulabelsNew(ii);
    
    % find the overlap of the unique number in ImageBWlabelOld
    overlapOld = ImageBWlabelOld(maskCell);

    % find the most frequent number in the overlap for the old label
    OldLabel = mode(overlapOld);
    
    % if there is no cell there then a new cell is being tracked
    if OldLabel < 1
        ImageBWlabelNew_Track(maskCell) = maxOldLabel + 1;
        maxOldLabel = maxOldLabel + 1;

    % otherwise use the old label in the tracked label maatrix
    else
        ImageBWlabelNew_Track(maskCell) = OldLabel;
        
    end

end
end % labelTracker

function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageV5(image,seeds_in,mask_in,background,SEerodeMask,shrinkVal)
% Perform seeded watershed segmentation on image given the input seeds and 
% mask logical arrays. This function will output two different types of
% seeds to be saved for segmentation of the next frame.
%      
% INPUT:    image: an image that has already been filtered in the parent
%                  function.
%           seeds_in: logical matrix of seed regions (same size as image).
%           mask_in: logical matrix of mask regions (same size as image).
%           background: the background of the image to be used in the case
%                       that anything is touching too close to the edge
%           SEerodeMask: structuring element for the erosion of seeds and
%                        mask.
%           shrinkVal: BWmorph shrink magnitude.
%
% OUTPUT:   imageSegmented: logical matrix where segmented regions are true and
%                           boundary lines (interfaces) are false.
%           imageBWlabel: an integer labeled image where cells have unique
%                         integer labels, interfaces are zeros, and background is -1.
%           seeds_out: logical matrix of new seed regions for cells.
%           mask_out: logical matrix of new background mask regions.

% find the number of connected objects in seeds and the minimum pixel number
CCin = bwconncomp(seeds_in);
numPixelsIn = cellfun(@numel,CCin.PixelIdxList);
minnumPin = min(numPixelsIn);

% throw an error if the minimum pixel is less than 10
if minnumPin < 10
    disp('Problem!!!')
end

% impose minima for both the watershed seeds and the mask
imageMasked = imimposemin(image, seeds_in | mask_in);

% perform watershed, set all labeled regions to 1, convert to logical (BW)
imageSegmented = watershed(imageMasked,8);
imageSegmented(imageSegmented>0) = 1;
imageSegmented = logical(imageSegmented);

% convert BW to label image
imageBWlabel = bwlabel(imageSegmented,8);

% any regions that touch the background or are within a certain distance to it
background(:,1:5) = true;
background(1:5,:) = true;
background(:,(end-4):end) = true;
background((end-4):end,:) = true;

% make the edge areas equal to 0
edgeAreas = unique(imageBWlabel(background));
for ie=1:length(edgeAreas)
    cvalue = edgeAreas(ie);
    if cvalue>0
        imageBWlabel(imageBWlabel==cvalue) = 0;
    end
end

% shrink seeds for propagation
seeds_out = bwmorph(imageBWlabel > 0,'shrink',shrinkVal);

% find the number of connected objects in seeds and the minimum pixel number
CCout = bwconncomp(seeds_out);
numPixelsOut = cellfun(@numel,CCout.PixelIdxList);
minnumPout = min(numPixelsOut);

% throw an error if the minimum pixel is less than 0
if  minnumPout < 0
    disp('Problem!!!')
end
    

% % shrunken seeds sometimes get very narrow for small or elongated cells and
% % become poor seeds. Identify these and dilate them.
seedsLabelmat = bwlabel(seeds_out);
seed_labelvec = unique(seedsLabelmat(:));
seeds_eroded = imerode(seeds_out,strel('disk',1));
eroded_labelvec = unique(seedsLabelmat(seeds_eroded));

smallseed_labels = setdiff(seed_labelvec,eroded_labelvec);
smallseed_labels(smallseed_labels<1) = [];
seeds_small = false(size(seeds_out));
for ii=1:numel(smallseed_labels)
    seeds_small(seedsLabelmat==smallseed_labels(ii)) = true;
end
seeds_small = imdilate(seeds_small,strel('disk',2));

seeds_out = seeds_out | seeds_small;

mask_out = ~imdilate(imageBWlabel>0,SEerodeMask);


end % of subfunction



function [seeds,mask,eflag,eflag2] = modifySeedsAndMask_LS(I,Ifilt,seeds,mask,titlestr,SEerodeMask,shrinkVal)
% this function allows the user to make easy changes to watershed seeds and
% background mask
%
% INPUT:    I:  The image to display while modifications are made with RGB
%               values
%           Ifilt:  The filtered image used for segmentation
%           seeds:  Logical matrix of seed regions for cells
%           mask:   Logical matrix of background mask regions
%           titlestr: The title to diplayed on the figure
%           SEerodeMask: structuring element for the erosion of seeds and
%           mask.
%           shrinkVal: BWmorph shrink magnitude.
%
% OUTPUT:   seeds: Logical matrix of new seed regions for cells
%           mask: Logical matrix of new background mask regions
%           eflag: exit flag for quitting the function
%           eflag2: exit flag for saving and quitting the function



% set the figure position
figPosition = [266 52 859 745];

% set the structuring element
SEseg    = strel('disk',1);

% set the seedfilter (gaussian filter)
seedfilter = 12; % was 10

% set exit flags to false
eflag = false;
eflag2 = false;
 
% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);

% set unsatisfactory to true
unsatisfactory = true;

% continue to loop while unsatisfactory is true
while unsatisfactory

    %seeds = bwareaopen(seeds,5);
    mask  = bwareaopen(mask,50);
    
    imageMasked = imimposemin(Ifilt, seeds|mask);
    L = watershed(imageMasked,8);
    seg = imdilate(L==0,SEseg);
    
    % update translucent RGB image with modseeds and modmask and display
    imgR = I;
    imgG = I;
    imgB = I;
    imgR(seg) = 1;
    imgG(seg) = 0;
    imgB(seg) = 0;
    imgG(seeds) = I(seeds)+ 0.5;
    imgR(mask) = 0.75*I(mask);
    imgG(mask) = 0.75*I(mask);
    imgB(mask) = I(mask)+ 0.5;
    imgrgb = cat(3,imgR,imgG,imgB);

    imshow(imgrgb)
    if nargin > 3
        title(titlestr,'FontSize',18)
    end
    set(gcf,'Position',figPosition)

    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Seed Convert','Mask Convert','Add Polygon Seed','Add Line Seed',...
        'Add Polygon Mask','Remove Polygon','Shift left','Undo','Save+Quit','Quit');
    
    % hold current seeds and mask in case changes need to be undone
    
    switch choice

        % Continue
        case 1
            unsatisfactory = false;
            
        % Mask to Seed    
        case 2
            last_seeds = seeds;
            last_mask = mask;
            
            
            imshow(imgrgb); 
            set(gcf,'Position',figPosition)
            
            [i2seeds] = im2seeds(I,seedfilter,0);
            imageMasked = imimposemin(Ifilt, i2seeds);
            % perform watershed
            L2 = watershed(imageMasked,8);
            
            title(gca,'Select Mask Region to Convert to Seed','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)
                 if mask(x(k),y(k)) % if mask region is clicked dilate it only
                    newseglabel = L2(x(k),y(k));
                    newsegmask  = L2 == newseglabel;
                    
                    BWdilate = imdilate(newsegmask,SEerodeMask);
                    newsegmask = bwmorph(newsegmask,'shrink',shrinkVal);
                    %BWerode = imerode(newsegmask,strel('disk',4));
                    mask(BWdilate) = false;
                    seeds(newsegmask) = true;
                 end
%                 clickedMask  = L2 == L2(x(k),y(k));
% 
%                 BWdilate = imdilate(clickedMask,SEerodeMask);
%                 clickedMask = bwmorph(clickedMask,'shrink',shrinkVal);
%                 %BWerode = imerode(newsegmask,strel('disk',4));
%                 mask(BWdilate) = false;
%                 seeds(BWdilate) = false;
%                 seeds(clickedMask) = true;

            end            
           
        % convert Region to Mask
        case 3
            last_seeds = seeds;
            last_mask = mask;
            
            maskLabels = unique(L(mask));
            
            imshow(imgrgb); 
            set(gcf,'Position',figPosition)
            title(gca,'Select Regions to Convert to Mask','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            for k=1:length(x)


                label = L(x(k),y(k));
                seeds(L==label) = false;
                maskLabels = [maskLabels;label];

            end
            mask = false(size(mask));
            for k=1:numel(maskLabels)
                mask(L==maskLabels(k)) = true;
            end
            mask = imclose(mask,strel('disk',1));
            mask = imerode(mask,SEerodeMask);

        % Add Polygon Seed    
        case 4
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);                
            seeds = seeds | BW;
            %seeds = bwmorph(seeds,'shrink',10);
            mask(imdilate(BW,SEerodeMask)) = false;

        % Add line seed
        case 5
            h = imline;
            position = wait(h);
            x1 = position(1,2);
            y1 = position(1,1);
            x2 = position(2,2);
            y2 = position(2,1);
            BW = lineMask(size(I), [x1 y1], [x2 y2],1);
            seeds = seeds | BW;
            mask(imdilate(BW,SEerodeMask)) = false;
            
        % Add Polygon Mask
        case 6 
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            mask = mask | BW; 
            seeds(imdilate(BW,SEerodeMask)) = false;
            
        % Remove Polygon (removes polygon from both seeds and mask)    
        case 7
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            seeds(BW) = false;
            mask(BW) = false;   
            
        % shift left    
        case 8
            last_seeds = seeds;
            last_mask = mask;
            
            seeds = circshift(seeds,-5,2);
            mask = circshift(mask,-5,2);

        % Undo last change
        case 9
            seeds = last_seeds;
            mask  = last_mask;

        % Save and quit function    
        case 10
            eflag2 = true;
            return

        % Exit Function
        case 11
            eflag = true;
            return

              
    end % switch

    
end
 
 
% close everything
close all

end

function  [seeds] = im2seeds(I,filtersize,viewResults)
%IM2SEEDS produces point-like watershed seeds for a 2D image based on
%regional minima of the filtered image.
%
% Input:    I:              2D image.
%           filtersize:     Gaussian filter size. Must have the function
%                           "filterImage3DpaddedEdges.m"
%           viewResults:    If set to 1 an overlay of the seeds will be
%                           shown
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
SE2 = strel('disk',7);
seeds = imdilate(IRM,SE2);
 
 
% if viewResults = 1 show overlay of the seeds
if nargin >2
    if viewResults
        imoverlay(I,seeds,[1 0 0]);
    end
end
end



