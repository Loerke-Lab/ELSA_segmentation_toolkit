function [] = ELSA_segmentCells_PropagateInZ_LS_Edits(data,MovieNum,LayerVec,t,cellNum)
% record original directory (and return to it at the end)
od = cd;
close all


erosionSE = strel('disk',1);
%filtering parameter
filtersize = [2 2 1];
pathImages = data(MovieNum).Source;

% get list of images
list = data(MovieNum).ImageFileList;
StartLayer = LayerVec(1);
LastLayer = LayerVec(end);
interval = LayerVec(2)-LayerVec(1);


% load Current image stack
cd(pathImages);
% I = tiffreadVolume(list{t});
for z = 1:size(list,2)
    I(:,:,z) = imread(list{t,z});
end
Ifilt = imgaussfilt3(I,filtersize);

% load Current Seeds, Mask, and Label matrix
%cd ..
cframefoldername = strcat('SegmentationData/',sprintf('frame%04d',t));
cd(cframefoldername);
seeds = load('seeds.mat').seeds;
mask = load('mask.mat').mask;
ImageSegment = load('ImageSegment.mat').ImageSegment;
ImageBWlabel = load('ImageBWlabel.mat').ImageBWlabel;


% imshow(ImageBWlabel(:,:,StartLayer)==cellNum);
figure('Units','normalized','Position',[0 0 1 .9])
    imshow(ImageBWlabel(:,:,StartLayer)==cellNum,'InitialMagnification','fit')
pause;
close all;
finishflag = false;
exitflag = false;

% Go down the stack starting from segmented layer
for zlayer = (StartLayer+interval):interval:LastLayer
    seedsZ = seeds(:,:,zlayer-interval);
    maskZ  = mask(:,:,zlayer-interval);
    titlestr = strcat('Layer ',num2str(zlayer));
    if ~finishflag
        [seedsZ,maskZ,finishflag,exitflag] = modifySeedsAndMask_LS(I(:,:,zlayer),Ifilt(:,:,zlayer),seedsZ,maskZ,titlestr);
    end
    if exitflag
        return
    end
    
    
    [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),seedsZ,maskZ);
    [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabel(:,:,zlayer-interval));
    
    Areas_PrevZ = cell2mat({regionprops(ImageBWlabel(:,:,zlayer-interval),'Area').Area});
    Areas = cell2mat({regionprops(ImageBWlabel_z,'Area').Area});
    Areas(Areas==0) = NaN;
    maxN = max([numel(Areas_PrevZ),numel(Areas)]);
    Areas((end+1):maxN) = NaN;
    Areas_PrevZ((end+1):maxN) = NaN;
    DiffZ = (Areas - Areas_PrevZ)./Areas_PrevZ*100;
    BigChangeOrTooSmall = find(abs(DiffZ) >35 | Areas < 500);
    % remove seeds and set ImageBWlabel to zero for regions that get too
    % small or have a large percent change from one z-layer to the next
%     if ~isempty(BigChangeOrTooSmall) 
%         for ii=1:numel(BigChangeOrTooSmall)
%             seedsZ(ImageBWlabel_z==BigChangeOrTooSmall(ii)) = false;
%             ImageBWlabel_z(ImageBWlabel_z==BigChangeOrTooSmall(ii)) = 0;
%         end
%         maskZ = imerode(ImageBWlabel_z==0,erosionSE);
%         [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),seedsZ,maskZ);
%         [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabel(:,:,zlayer-interval));
%     end
% I took this out because it bugged me -Noah    

    seeds(:,:,zlayer) = seeds_z;
    mask(:,:,zlayer) = mask_z;
    ImageBWlabel(:,:,zlayer) = ImageBWlabel_z;
    ImageSegment(:,:,zlayer) = ImageSegment_z;
    
end


% save results into the current results folder
save('ImageSegment', 'ImageSegment');
save('ImageBWlabel', 'ImageBWlabel');
save('seeds','seeds');
save('mask','mask');





cd(od);

end % of function

%% =======================================================================
%
%                               subfunctions
%
%  =======================================================================

function [ImageBWlabelNew_Track] = labelTracker(ImageBWlabelNew,ImageBWlabelOld)

maxOldLabel = max(unique(ImageBWlabelOld(:)));


ImageBWlabelNew_Track = ImageBWlabelNew;

ulabelsNew = unique(ImageBWlabelNew(:));
ulabelsNew(ulabelsNew<1) = [];

NumNew = numel(ulabelsNew);

for ii=1:NumNew
    maskCell = ImageBWlabelNew == ulabelsNew(ii);
    
    overlapOld = ImageBWlabelOld(maskCell);
    OldLabel = mode(overlapOld);
    
    if OldLabel < 1 % Then tracking new cell
        ImageBWlabelNew_Track(maskCell) = maxOldLabel + 1;
        maxOldLabel = maxOldLabel + 1;
    else
        ImageBWlabelNew_Track(maskCell) = OldLabel;
        
    end
end
end % labelTracker

function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageV5(image,seeds_in,mask_in)
% Perform seeded watershed segmentation on image given the input seeds and
% mask logical arrays. This function will output two different types of
% seeds to be saved for segmentation of the next frame
% Inputs:
%           image: an image that has already been filtered in the parent
%           function.
%           seeds_in: logical matrix of seed regions (same size as image).
%           mask_in: logical matrix of mask regions (same size as image).



% Structuring element for the erosion of seeds and mask
SEerosionMask = strel('disk',10);
%SEclose   = strel('disk',2);

% impose minima for both the watershed seeds and the mask
imageMasked = imimposemin(image, seeds_in | mask_in);

% perform watershed, set all labeled regions to 1, convert to logical (BW)
imageSegmented = watershed(imageMasked,8);
imageSegmented(imageSegmented>0) = 1;
imageSegmented = logical(imageSegmented);

% convert BW to label image
imageBWlabel = bwlabel(imageSegmented,8);

% any regions that touch the background or are within a certain distance to
% it
background = false(size(mask_in));
background(:,1) = true;
background(1,:) = true;
background(:,end) = true;
background(end,:) = true;
backdist = bwdist(background);
edgepixels = backdist <= 3;
edgepixels(mask_in) = true;
edgeAreas = unique(imageBWlabel(edgepixels));
for ie=1:length(edgeAreas)
    cvalue = edgeAreas(ie);
    if cvalue>0
        imageBWlabel(imageBWlabel==cvalue) = 0;
    end
end

% shrink seeds for propagation
seeds_out = bwmorph(imageBWlabel > 0,'shrink',10);

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


% as segmentation progresses towards the more basal regions of the tissue
% seeds tend to collapse and the mask takes over. This is evident by a
% rapid reduction in size and possibly also a reduced signal surrounding
% the seed. Determine a method of removing these seeds to prevent
% segmentation errors. They also tend to be edge cells.




mask_out = ~imdilate(imageBWlabel>0,SEerosionMask);
end % of subfunction

function [seeds,mask,finishflag,exitflag] = modifySeedsAndMask_LS(I,Ifilt,seeds,mask,titlestr)

SEpoints = strel('disk',7);
SEmorph  = strel('disk',3);
SEseg    = strel('disk',1);
SEcloseMask   = strel('disk',10);
SEerodeMask = strel('disk',10);
seedfilter = 8;
magVal = 240;
finishflag = false;
exitflag = false;

% rescale image (if not already the case)
I = mat2gray(I);
I = imadjust(I);


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
    
%     imshow(imgrgb,'InitialMagnification',magVal)
    figure('Units','normalized','Position',[0 0 1 .9])
    imshow(imgrgb,'InitialMagnification','fit')
    if nargin > 3
        title(titlestr,'FontSize',18)
    end
    
    % Choose an option for modifying current layer
    choice = menu('Modification:','Continue','Seed Convert','Mask Convert','Add Polygon Seed',...
        'Add Polygon Mask','Remove Polygon','Add Point Seeds','Initialize','Line Seed','Undo','Finish stack','Quit');
    
    % hold current seeds and mask in case changes need to be undone
    
    switch choice
        case 1 % Continue
            unsatisfactory = false;
            
        case 2 % Mask to Seed
            last_seeds = seeds;
            last_mask = mask;
            
%             imshow(imgrgb,'InitialMagnification',magVal);
            figure('Units','normalized','Position',[0 0 1 .9])
    imshow(imgrgb,'InitialMagnification','fit')
            
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
            
%             imshow(imgrgb,'InitialMagnification',magVal);
figure('Units','normalized','Position',[0 0 1 .9])
    imshow(imgrgb,'InitialMagnification','fit')
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
            %mask = imerode(mask,SEerodeMask);
            
            
        case 4 % Add Polygon Seed
            last_seeds = seeds;
            last_mask = mask;
            
            BW = roipoly();
            title(gca,'Add Polygon Seed','FontSize',16);
            seeds = seeds | BW;
            %seeds = bwmorph(seeds,'shrink',10);
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
            
        case 9
            last_seeds = seeds;
            last_mask = mask;
            
            h = imline;
            position = wait(h);
            x1 = position(1,2);
            y1 = position(1,1);
            x2 = position(2,2);
            y2 = position(2,1);
            BW = lineMask(size(I), [x1 y1], [x2 y2],1);
            seeds = seeds | BW;
            
        case 10 % Undo last change
            seeds = last_seeds;
            mask  = last_mask;
        case 11
            finishflag = true;
            return
        case 12
            exitflag = true;
            return
            
            
            
    end % switch
end
close all
end

