function []= ELSA_initialize3Dseeds_SingleLayer_LS(data,t,layer,seedfilter)
%ELSA_initialize3Dseeds_SingleLayer_LS allows the user to initialize the
%seeds and mask for a single layer for 3D image data.
%
% INPUT:    data: structure that contains the fields 'ImageFileList' and
%                 'Source'.
%           t: the time frame to initialize the seeds
%           layer: the z layer to initialize the seeds (if movie is only 2D
%                  this will be 1)
%           seedfilter: Gaussian filter size. Must have the function
%                       "filterImage3DpaddedEdges.m"
%
% OUTPUT:   No output but the results are saved to the SegmentationData
%           folder associated with the input time frame
 
% get original directory
od = cd;

% get information from data structure
list = data.ImageFileList;
pathimages = data.Source; 

% move to the image path
cd(pathimages);

% get image data to use for initializing seeds
Ifull = tiffreadVolume(list{t});
for z = 1:size(list,2)
    Ifull(:,:,z) = imread(list{t,z});
end

% move to frame specific segmentation folder
cd('SegmentationData');
cframefoldername = sprintf('frame%04d',t);
cd(cframefoldername);

% load seeds, mask, ImageBWlabel, and ImageSegment
seeds = load('seeds.mat').seeds;
mask = load('mask.mat').mask;
ImageBWlabel = load('ImageBWlabel.mat').ImageBWlabel;
ImageSegment = load('ImageSegment.mat').ImageSegment;


% ImBoundarySeed = Ifull==0;
% 
% ImBoundarySeed(:,1) = true;
% ImBoundarySeed(1,:) = true;
% ImBoundarySeed(:,end) = true;
% ImBoundarySeed(end,:) = true;

% filter the image
Ifull_Filt = imgaussfilt3(Ifull,[2 2 1]);

% get number of layers
Nlayers = size(Ifull,3);


% get specific layer for image and filtered image specified by layer input   
Ilayer= Ifull(:,:,layer);
Ilayer_Filt= Ifull_Filt(:,:,layer);

% create title string to display during seed initialization
titlestr = strcat('Layer: ',num2str(layer),', Depth: ',num2str((layer-1)*0.5));

% initialize the seeds and mask
[seeds(:,:,layer),mask(:,:,layer),ImageBWlabel(:,:,layer),ImageSegment(:,:,layer)] = modifySeedsAndMask_LS(Ilayer,Ilayer_Filt,seeds(:,:,layer),mask(:,:,layer),titlestr,seedfilter);
    

% save seeds, mask, ImageBWlabel, and ImageSegment   
save('ImageSegment', 'ImageSegment');
save('ImageBWlabel', 'ImageBWlabel');
save('seeds','seeds');
save('mask','mask');

% return to original directory
cd(od)
 
% close figure 
close(gcf)

end


 
function [seeds,mask,imageBWlabel,imageSegmented] = modifySeedsAndMask_LS(I,Ifilt,seeds,mask,titlestr,seedfilter)
% this function allows the user to make easy changes to watershed seeds and
% background mask
%
% INPUT:    I:  The image to display while modifications are made with RGB
%               values
%           Ifilt:  The filtered image used for segmentation
%           seeds:  Logical matrix of seed regions for cells
%           mask:   Logical matrix of background mask regions
%           titlestr: The title to diplayed on the figure
%           seedfilter: Gaussian filter size. Must have the function
%                       "filterImage3DpaddedEdges.m"
%
% OUTPUT:   seeds: Logical matrix of new seed regions for cells
%           mask: Logical matrix of new background mask regions
%           imageBWlabel: An integer labeled image where cells have unique
%                         integer labels, interfaces are zeros, and background is -1
%           imageSegmented: Logical matrix where segmented regions are true and
%                           boundary lines (interfaces) are false




% set structuring elements
SEpoints = strel('disk',7);
SEmorph  = strel('disk',3);
SEseg    = strel('disk',1);
SEcloseMask   = strel('disk',10);
SEerodeMask = strel('disk',10);

% set magnification value
magVal = 200;

% set exit flag to false
eflag = false;
 
% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);

 
% initialize the modified seeds and mask to the input seeds and mask
seg = false(size(seeds));
 
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

    imshow(imgrgb)%,'InitialMagnification',magVal)
    if nargin > 3
        title(titlestr,'FontSize',18)
    end
    
    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Seed Convert','Mask Convert','Add Polygon Seed',...
        'Add Polygon Mask','Remove Polygon','Add Point Seeds','Initialize','Undo','Quit');
    
    % hold current seeds and mask in case changes need to be undone
    
    switch choice

        % Continue
        case 1
            unsatisfactory = false;
        
        % Mask to Seed    
        case 2
            last_seeds = seeds;
            last_mask = mask;
            
            imshow(imgrgb,'InitialMagnification',magVal); 
            
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
                    
                    BWdilate = imdilate(newsegmask,strel('disk',4));
                    newsegmask = bwmorph(newsegmask,'shrink',10);
                    %BWerode = imerode(newsegmask,strel('disk',4));
                    mask(BWdilate) = false;
                    seeds(newsegmask) = true;
                elseif seeds(x(k),y(k)) % if seed region is clicked dilate it only
                    L2 = bwlabel(seeds);
                    label = L2(x(k),y(k));
                    BW = L2 == label;
                    BWdilate = imdilate(BW,strel('disk',1));
                    seeds(BWdilate) = true;
                else % if no region is clicked
                    newseglabel = L2(x(k),y(k));
                    newsegmask  = L2 == newseglabel;
                    
                    BWdilate = imdilate(newsegmask,strel('disk',4));
                    newsegmask = bwmorph(newsegmask,'shrink',10);
                    %BWerode = imerode(newsegmask,strel('disk',4));
                    mask(BWdilate) = false;
                    seeds(BWdilate) = false;
                    seeds(newsegmask) = true;
                end
            end
            %bwmorph(seeds,'shrink',10);
            
        % convert Seed to Mask    
        case 3
            last_seeds = seeds;
            last_mask = mask;
            
            imshow(imgrgb,'InitialMagnification',magVal); 
            title(gca,'Select Regions to Convert','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            BW = false(size(I));
            for k=1:length(x)
                if seeds(x(k),y(k)) % 

                    label = L(x(k),y(k));
                    BW(L == label) = true;
                    seeds(L==label) = false;
                end
            end
            mask = imclose(BW | mask,SEcloseMask);
            mask = imerode(mask,SEerodeMask);

        % Add Polygon Seed    
        case 4
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);                
            seeds = seeds | BW;
            seeds = bwmorph(seeds,'shrink',10);
            mask(imdilate(BW,SEerodeMask)) = false;
            
        % Add Polygon Mask    
        case 5 
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            mask = mask | BW; 
            seeds(imdilate(BW,SEerodeMask)) = false;
        
        %  Remove Polygon (removes polygon from both seeds and mask)    
        case 6
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            seeds(BW) = false;
            mask(BW) = false;            
        
        % Add Point Seeds (can continue adding until user hits enter)    
        case 7
            last_seeds = seeds;
            last_mask = mask;
            
            title(gca,'Add Point Seeds','FontSize',16);
            [y, x] = ginput;
            x = round(x); y = round(y);
            BW = false(size(I));
            for k=1:length(x)
                BW(x(k),y(k)) = true;
            end
            % slightly dilate so seed is not single pixel
            BW = imdilate(BW,SEpoints);
            seeds = seeds | BW;
            
        % Re-initialize
        case 8
            last_seeds = seeds;
            last_mask = mask;
            
            seeds = im2seeds(I,seedfilter,0);
            
            imageMasked = imimposemin(Ifilt, seeds);
            Lseeds = watershed(imageMasked,8);
            
            edgepixels1 = Lseeds(:,1:4);
            edgepixels2 = Lseeds(1:4,:);
            edgepixels3 = Lseeds(:,(end-3):end);
            edgepixels4 = Lseeds((end-3):end,:);
            %contiguous areas that touch the edge    
            edgeLabels = unique( [edgepixels1(:);edgepixels2(:);edgepixels3(:);edgepixels4(:)] );
            edgeLabels(edgeLabels<1) = [];
            
            mask = false(size(mask));
            for ie=1:length(edgeLabels)
                mask(Lseeds==edgeLabels(ie)) = true;
                seeds(Lseeds==edgeLabels(ie)) = false;
                Lseeds(Lseeds==edgeLabels(ie)) = 0;
            end
            mask = imclose(mask,strel('disk',2));
            mask = imerode(mask,SEerodeMask);
            
            seedsEr = imerode(Lseeds>0,SEmorph);
            seeds = bwmorph(Lseeds>0,'shrink',10);
            
            
            LseedsPerim = Lseeds;
            LseedsPerim(seedsEr) = 0;
            LseedsCenter = Lseeds;
            LseedsCenter(~seedsEr) = 0;
            StatsP = regionprops(LseedsPerim,Ifilt,'MeanIntensity');
            StatsC = regionprops(LseedsCenter,Ifilt,'MeanIntensity');
            MIP = [StatsP.MeanIntensity];
            MIC = [StatsC.MeanIntensity];
            
            MeanIRatio = MIP./MIC;
            idxWeakRatio = find(MeanIRatio< 1.2);

            for ii=1:numel(idxWeakRatio)
                seeds(Lseeds==idxWeakRatio(ii)) = 0;
            end
        
        % Undo last change
        case 9
            seeds = last_seeds;
            mask  = last_mask;

        % Exit Function    
        case 10
            eflag = true;
            return

                
 
    end % switch


    % after modifications have been made use the new seeds and mask to
    % generate a new segmentation and new "shrunken" seeds
    imageMasked = imimposemin(Ifilt, seeds|mask);
    L = watershed(imageMasked,8);
    imageSegmented = L;
    imageSegmented(imageSegmented>0) = 1;
    imageSegmented = logical(imageSegmented);
    seg = imdilate(L==0,SEseg);


    maskLabels = unique(L(mask));
    for ie=1:length(maskLabels)
        mask(L==maskLabels(ie)) = true;
        L(L==maskLabels(ie)) = false;
    end
    mask = imclose(mask,strel('disk',2));
    mask = imerode(mask,SEerodeMask);
    seeds = L > 0;

    seeds = bwmorph(seeds,'shrink',10);
    
    imageBWlabel = bwlabel(imageSegmented,8);


    % Note that any areas that touch the edge of the image are set to value
    % -1 to mark them as background - this will ensure in the subsequent
    % analysis that these cells are not used for direct characterization of
    % nodes, since they do not yield full neighborhood information
    edgepixels1 = imageBWlabel(:,1:4);
    edgepixels2 = imageBWlabel(1:4,:);
    edgepixels3 = imageBWlabel(:,(end-3):end);
    edgepixels4 = imageBWlabel((end-3):end,:);
    %contiguous areas that touch the edge    
    edgeAreas = unique( [edgepixels1(:);edgepixels2(:);edgepixels3(:);edgepixels4(:)] );
    % set these areas to value -1
    for ie=1:length(edgeAreas)
        cvalue = edgeAreas(ie);
        if cvalue>0
            imageBWlabel(imageBWlabel==cvalue) = -1;
        end
    end

    % Besides the edges there may be masks that do not touch edges, so find
    % them from the input mask and set those regions to -1 also.
    maskAreas = unique(imageBWlabel(mask));
    for ie=1:length(maskAreas)
        cvalue = maskAreas(ie);
        if cvalue>0
            imageBWlabel(imageBWlabel==cvalue) = 0;
        end
    end
    
    imageBWlabel = bwlabel(imageBWlabel>0);
    
    
end
 
 
% close everything 
close all

end



% function [seeds,mask] = removeBoundaryCells(ImageBoundarySeed,L,seeds,mask)
% DistThreshBoundary = 40;
% DistThreshSegLines = 6;
% ImageBoundarySeed = bwareaopen(ImageBoundarySeed,10);
% maskLabels = unique(L(mask));
% for ii=1:numel(maskLabels)
%     L(L==maskLabels(ii)) = 0;
% end
% 
% uniqueL = unique(L);
% uniqueL(uniqueL==0) = [];
% NumL = max(uniqueL);
% Lseeds = L;
% Lseeds(~seeds) = 0;
% 
% 
% DistFromBoundary = bwdist(ImageBoundarySeed);
% DistFromSegLines = bwdist(L==0);
% 
% StatsBoundary = regionprops(L,DistFromBoundary,'MinIntensity');
% StatsSeg = regionprops(Lseeds,DistFromSegLines,'MinIntensity');
% 
% 
% MIB = [StatsBoundary.MinIntensity];
% MinIntensityB = NaN(NumL,1);
% MinIntensityB(uniqueL) = MIB;
% MIS = [StatsSeg.MinIntensity];
% MinIntensityS = NaN(NumL,1);
% MinIntensityS(uniqueL) = MIS;
% 
% idx= find(MinIntensityB < DistThreshBoundary & MinIntensityS < DistThreshSegLines);
% 
% if ~isempty(idx)
%     for ii = 1:numel(idx)
%         seeds(L==idx(ii)) = false;
%         mask(L==idx(ii)) = true;
%     end
% end
% 
% 
% end
 
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


