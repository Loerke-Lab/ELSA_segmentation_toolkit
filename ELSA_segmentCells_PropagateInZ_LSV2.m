function [] = ELSA_segmentCells_PropagateInZ_LSV2(data,MovieNum,StartLayer,tvec)
%ELSA_segmentCells_PropagateInZ_LSV2 propagates the segmentation from
%"StartLayer" both down to the "LastLayer" (Note: the last layer must be stored
%in data under the field "Nlayers", and up to layer 1.

% record original directory (and return to it at the end)
od = cd;

erosionSE = strel('disk',10);
SEdarkzone = strel('disk',24);
% zlayers to be processed, path to images, and list of image file names
FirstLayer = 1;
LastLayer = data(MovieNum).Nlayers;
pathImages = data(MovieNum).Source;
list = data(MovieNum).ImageFileList;


filtersize = [2 2 1];%filtering parameter



% loop over all timepoints
for ti = 1:length(tvec)
    
    t = tvec(ti);
    % display current progress of processing
    fprintf('segmentation @ timepoint %04d',t);
    
    % load Current image stack and 3D filter
    cd(pathImages);
    I = tiffreadVolume(list{t});
%     for z = 1:size(list,2)
%         I(:,:,z) = imread(list{t,z});
%     end
    Ifilt = imgaussfilt3(I,filtersize);
    
    % load Current Seeds, Mask, and Label matrix
    %cd ..
    cframefoldername = strcat('SegmentationData/',sprintf('frame%04d',t));
    cd(cframefoldername);
    seeds = load('seeds.mat').seeds;
    mask = load('mask.mat').mask;
    ImageSegment = load('ImageSegment.mat').ImageSegment;
    ImageBWlabel = load('ImageBWlabel.mat').ImageBWlabel;


    % Go down the stack starting from segmented layer
    for zlayer = (StartLayer+1):LastLayer
        seedsZ = seeds(:,:,zlayer-1);
        maskZ  = mask(:,:,zlayer-1);
        
        % the darkzone is the name I call the region of zeros in light
        % sheet data that are a result of deskewing. As you go down the
        % z-stack the darkzone increases and it is essentially the new edge
        % of the field of view. I dilate this region and use its mask to
        % remove cells that intersect the darkzone (as we would normally
        % remove cells that intersect the image boundary). If you don't
        % have a darkzone then set it to false (e.g. darkzone =
        % darkzone*false;).
        darkzone = imdilate(Ifilt(:,:,zlayer)==0,SEdarkzone);
        
       
        [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),seedsZ,maskZ,darkzone);
        [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabel(:,:,zlayer-1));
        
        imoverlay(imadjust(Ifilt(:,:,zlayer)),~ImageSegment_z);
        title(strcat('Frame ',num2str(t),', Layer ',num2str(zlayer)))
        
        % The following commented lines were intended to remove weird cells
        % that changed area by a large amount or were too small.
%         Areas_PrevZ = cell2mat({regionprops(ImageBWlabel(:,:,zlayer-1),'Area').Area});
%         %Areas_PrevT = cell2mat({regionprops(ImageBWlabel_Prev(:,:,zlayer),'Area').Area});
%         Areas = cell2mat({regionprops(ImageBWlabel_z,'Area').Area});
%         Areas(Areas==0) = NaN;
%         maxN = max([numel(Areas_PrevZ),numel(Areas)]);
%         Areas((end+1):maxN) = NaN;
%         Areas_PrevZ((end+1):maxN) = NaN;
%         %Areas_PrevT((end+1):maxN) = NaN;
%         %DiffT = (Areas - Areas_PrevT)./Areas_PrevT;
%         DiffZ = (Areas - Areas_PrevZ)./Areas_PrevZ*100;
% %         subplot(2,2,1);imshow(label2rgb(ImageBWlabel_z));title('Current')
% %         subplot(2,2,2);imshow(label2rgb(ImageBWlabel_Prev(:,:,zlayer)));title('Prev. T')
% %         subplot(2,2,3);imshow(label2rgb(ImageBWlabel(:,:,zlayer-1)));title('Prev. Z')
% %         subplot(2,2,4);plot(DiffZ);hold on;plot(DiffT);legend('DiffZ','DiffT');hold off
%         
%         BigChangeOrTooSmall = find(abs(DiffZ) > 20 | Areas < 500);
%         % remove seeds from
%         if ~isempty(BigChangeOrTooSmall)
%             for ii=1:numel(BigChangeOrTooSmall)
%                 seedsZ(ImageBWlabel_z==BigChangeOrTooSmall(ii)) = false;
%                 ImageBWlabel_z(ImageBWlabel_z==BigChangeOrTooSmall(ii)) = 0;
%             end
%             maskZ = imerode(ImageBWlabel_z==0,erosionSE);
%             
% 
%             [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),seedsZ,maskZ,darkzone);
%             [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabel(:,:,zlayer-1));
%         end
            

        % store z-layer results for seeds, mask, etc.
        seeds(:,:,zlayer) = seeds_z;
        mask(:,:,zlayer) = mask_z;
        ImageBWlabel(:,:,zlayer) = ImageBWlabel_z;
        ImageSegment(:,:,zlayer) = ImageSegment_z;

    end
    % Go up the stack starting from segmented layer
    for zlayer = (StartLayer-1):-1:FirstLayer
        seedsZ = seeds(:,:,zlayer+1);
        maskZ  = mask(:,:,zlayer+1);
        
        darkzone = imdilate(Ifilt(:,:,zlayer)==0,SEdarkzone);
        
        [ImageSegment_z,ImageBWlabel_z,seeds_z,mask_z] = wsSegmentSingleImageV5(Ifilt(:,:,zlayer),seedsZ,maskZ,darkzone);
        [ImageBWlabel_z] = labelTracker(ImageBWlabel_z,ImageBWlabel(:,:,zlayer+1));
        
        imoverlay(imadjust(Ifilt(:,:,zlayer)),~ImageSegment_z);
        title(strcat('Frame ',num2str(t),', Layer ',num2str(zlayer)))
        

        
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
  
       
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for t

fprintf('\n');
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




function [imageSegmented,imageBWlabel,seeds_out,mask_out] = wsSegmentSingleImageV5(image,seeds_in,mask_in,darkzone)
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
% it will be set to zero. Also we're setting any regions taht overlap with
% the "darkzone" (i.e. the empty regions in light sheet data due to
% deskewing).
edgepixels = false(size(mask_in));
edgepixels(:,1:5) = true;
edgepixels(1:5,:) = true;
edgepixels(:,(end-4):end) = true;
edgepixels((end-4):end,:) = true;
edgepixels(darkzone|mask_in) = true;
edgeAreas = unique(imageBWlabel(edgepixels));
for ie=1:length(edgeAreas)
    cvalue = edgeAreas(ie);
    if cvalue>0
        imageBWlabel(imageBWlabel==cvalue) = 0;
    end
end

% shrink seeds for propagation
seeds_out = bwmorph(imageBWlabel > 0,'shrink',5);
mask_out = ~imdilate(imageBWlabel>0,SEerosionMask);
end % of subfunction

