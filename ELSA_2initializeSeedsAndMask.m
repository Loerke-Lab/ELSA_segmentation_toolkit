function [] = ELSA_2initializeSeedsAndMask(data,timeFrame)
% This function initializes the seeds and mask for the watershed segmentation.
% Usually initialization is done just once at the first time frame (t=1).
% The segmentation algorithm then propagates the seeds and mask for all
% subsequent frames. This function saves the seeds and mask to the
% SegmentationData folder associated with the input time frame (there is no
% output). If seeds have already been initialized then they will be loaded 
% and can be edited.
%
% INPUT:
%           data: structure that contains the fields 'ImageFileList' and
%                'Source'.
%           timeFrame: the time frame to initialize (default is 1)


% record original directory (and return to it at the end)
od = cd;


% if the initialization time frame is not given then set to 1
if nargin < 2
    timeFrame = 1;
end

% load the image 
image = imread(data.ImageFileList{timeFrame});
image = imadjust(image);


% move to directory for loading/saving seeds and mask
cd(data.Source)
cframefoldername = sprintf('SegmentationData/frame%04d',timeFrame);
cd(cframefoldername);



% if a seeds.mat file is already present then load it and mask.mat,
% else initialize new seeds and mask
if exist('seeds.mat','file')   
    loadseeds = load('seeds.mat');
    seeds = loadseeds.seeds;
    loadmask = load('mask.mat');
    mask = loadmask.mask;
else % initialize new (all false) seeds and mask
    seeds = false(size(image));
    mask  = false(size(image));
end


% if the seeds are all false then initialize automated seeds and mask
if all(~seeds,'all')
    [seeds,mask] = im2seedsAndMask(image,9);

end

% edit seeds and mask if necessary
[seeds,mask] = modifySeedsAndMask_ForInitialization(image,seeds,mask);


% save to current directory
save('seeds','seeds')
save('mask','mask')


% switch back to original directory
cd(od)
end


function  [seeds,mask] = im2seedsAndMask(I,filtersize)
%this function generates automated watershed seeds for a 2D image based on
%regional minima of the filtered image. It also generates a mask based on
%all the seeds that touch the image border.
%
% Input:    I:              2D image.
%           filtersize:     Filter size for use with imregionalmin so that each cell has
%                           approximately 1 minima. Must have the function
%                           "filterImage3DpaddedEdges.m"
%
% Output:   seeds = mask of points of minima in the image for seeded
%                   watershed.
%           mask  = mask of background (typically cells on the image
%                   boundary.
%
% updated 11/30/21 by Tim Vanderleest



% filter image for seeding and for segmentation
[Iseedfiltered] = filterImage3DpaddedEdges(I, 'Gauss', filtersize);

% find the regional minimum in the filtered image
imregmin = imregionalmin(Iseedfiltered);

% slightly dilate these points
SEdisk5 = strel('disk',5);
seeds = imdilate(imregmin,SEdisk5);

% find all segments that are close to the boundary and remove these seeds
[Isegfiltered] = filterImage3DpaddedEdges(I, 'Gauss', 2);
imageMasked = imimposemin(Isegfiltered, seeds);
imageSegmented = watershed(imageMasked,8);
edgepixels1 = imageSegmented(:,1:4);
edgepixels2 = imageSegmented(1:4,:);
edgepixels3 = imageSegmented(:,(end-3):end);
edgepixels4 = imageSegmented((end-3):end,:); 
edgeLabels = unique( [edgepixels1(:),edgepixels2(:),edgepixels3(:),edgepixels4(:)] );
for ie=1:length(edgeLabels) % loop through edge labels and set to zero
    labelvalue = edgeLabels(ie);

    imageSegmented(imageSegmented==labelvalue) = 0;

end

% use the segmented label image to create seeds via shrinking
seeds = bwmorph(imageSegmented>0,'shrink',7);

% a mask of the background is created by the negation of the dilated
% segments
mask = ~imdilate(imageSegmented>0,strel('disk',7));

end


function [seeds_out,mask_out] = modifySeedsAndMask_ForInitialization(I,seeds_in,mask_in)
%Allows a user to easily make changes to watershed seeds
%and background mask.
%
% INPUTS:       I:      The image for verifying seeds and masks.
%               seeds_in:  Watershed seeds (2D binary).
%               mask_in:   Mask of background region or regions (2D binary).
%
% OUTPUTS:      seeds_out: Modified seeds.
%               mask_out:  Modified mask.


% structuring element used below
SEdisk10 = strel('disk',10);

% imshow magnification value
magVal = 240;

% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);

% initialize the modified seeds and mask to the input seeds and mask
seeds_out = seeds_in;
mask_out = mask_in;

% generate segmentation to see segmentation overlay along with seeds and
% mask
[Ifilt] = filterImage3DpaddedEdges(I, 'Gauss', 2);
imageMasked = imimposemin(Ifilt, seeds_out|mask_out);
imageSegmented = watershed(imageMasked,8);
segboundarylines = imageSegmented<1;
%set mask regions to zero
maskLabels = unique(imageSegmented(mask_out));
for ii=1:length(maskLabels)
    imageSegmented(imageSegmented==maskLabels(ii)) = 0;
end


unsatisfactory = true;
while unsatisfactory
    
    % generate an image to show the current state of the seeds, mask, and
    % segmentation. Seeds in green, mask in blue, and segmentation lines in
    % red.
    imgRGB = repmat(I,1,1,3);
    mask1 = cat(3,segboundarylines,seeds_out,mask_out); 
    mask2 = cat(3,false(size(segboundarylines)),segboundarylines,segboundarylines);
    imgRGB(mask1) = imgRGB(mask1)*0.3 + 0.7; 
    imgRGB(mask2) = 0; 
    imshow(imgRGB,'InitialMagnification',magVal)
    

    % Choose an option for modifying current layer
    choice = menu('Modification:','Done',...
        'Add Polygon Seed','Add Polygon Mask','Remove Polygon',...
        'Remove by Clicking','Re-initialize Seeds','Shrink Seeds','Clear All');
    
    switch choice
        case 1 % Done
            unsatisfactory = false;
            
            
        case 2 % Add Polygon Seed
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);                
            seeds_out = seeds_out | BW;
            
            % also if new seed overlaps with mask convert mask region to false 
            mask_out(imdilate(BW,SEdisk10)) = false;
            
        case 3 % Add Polygon Mask    
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            mask_out = mask_out | BW;
            
            % also if new mask overlaps with seed convert seed region to false 
            seeds_out(imdilate(BW,SEdisk10)) = false;
            
        case 4 % Remove Polygon (removes polygon from both seeds and mask)
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            seeds_out(BW) = false;
            mask_out(BW) = false;
            
        case 5  % Remove Regions 
            title(gca,'Remove Regions','FontSize',16);
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

         
        case 6 % Re-initialize seeds and mask
            [seeds_out,mask_out] = im2seedsAndMask(I,9);
        case 7 % re-shrink seeds
            seeds_out = bwmorph(imageSegmented>0,'shrink',7);
    
            
            
        case 8 % clear all
            seeds_out = false(size(I));
            mask_out  = false(size(I));
            
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
close all
end
