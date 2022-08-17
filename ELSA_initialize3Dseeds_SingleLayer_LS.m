function []= ELSA_initialize3Dseeds_SingleLayer_LS(data,t,layer,seedfilter)
%ELSA_initialize3Dseeds_SingleLayer_LS allows the user to initialize the
%seeds and mask for a single layer for 3D image data.
 
% get original directory
od = cd;
 
list = data.ImageFileList;
pathimages = data.Source; 
cd(pathimages);
Ifull = tiffreadVolume(list{t});
for z = 1:size(list,2)
    Ifull(:,:,z) = imread(list{t,z});
end


cd('SegmentationData');
cframefoldername = sprintf('frame%04d',t);
cd(cframefoldername);
% load seeds and mask

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

Ifull_Filt = imgaussfilt3(Ifull,[2 2 1]);

Nlayers = size(Ifull,3);



    
Ilayer= Ifull(:,:,layer);
Ilayer_Filt= Ifull_Filt(:,:,layer);
titlestr = strcat('Layer: ',num2str(layer),', Depth: ',num2str((layer-1)*0.5));
[seeds(:,:,layer),mask(:,:,layer),ImageBWlabel(:,:,layer),ImageSegment(:,:,layer)] = modifySeedsAndMask_LS(Ilayer,Ilayer_Filt,seeds(:,:,layer),mask(:,:,layer),titlestr,seedfilter);
    

    
save('ImageSegment', 'ImageSegment');
save('ImageBWlabel', 'ImageBWlabel');
save('seeds','seeds');
save('mask','mask');

% return to original directory
cd(od)
 
 
close(gcf)
end


 
function [seeds,mask,imageBWlabel,imageSegmented] = modifySeedsAndMask_LS(I,Ifilt,seeds,mask,titlestr,seedfilter)

SEpoints = strel('disk',7);
SEmorph  = strel('disk',3);
SEseg    = strel('disk',1);
SEcloseMask   = strel('disk',10);
SEerodeMask = strel('disk',10);

magVal = 200;
eflag = false;
 
% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);

 
% initialize the modified seeds and mask to the input seeds and mask
seg = false(size(seeds));
 

unsatisfactory = true;
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
        case 1 % Continue
            unsatisfactory = false;
            
        case 2 % Mask to Seed
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
            
            
        case 3 % convert Seed to Mask
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

            
        case 4 % Add Polygon Seed
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);                
            seeds = seeds | BW;
            seeds = bwmorph(seeds,'shrink',10);
            mask(imdilate(BW,SEerodeMask)) = false;
            
        case 5 % Add Polygon Mask 
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Mask','FontSize',16);                
            mask = mask | BW; 
            seeds(imdilate(BW,SEerodeMask)) = false;
            
        case 6 % Remove Polygon (removes polygon from both seeds and mask)
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Remove Polygon','FontSize',16);                
            seeds(BW) = false;
            mask(BW) = false;            
            
        case 7 % Add Point Seeds (can continue adding until user hits enter).
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
            

        case 8 % Re-initialize
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
            
        case 9 % Undo last change
            seeds = last_seeds;
            mask  = last_mask;
        case 10 % Exit Function
            eflag = true;
            return

                
 
    end % switch

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


