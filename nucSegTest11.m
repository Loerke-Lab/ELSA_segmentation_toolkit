function [] = nucSegTest11(data,MovieNum,tRange,filename,solthresh,zscale)
% Automatically segments nuclei
%
%% Inputs
% data: The data structure used by all segmentation functions
% MovieNum: Typically 1
% tRange: Vector of frames to segment
% filename: Name of cell  labe matrix, such as 'ImageBWlabel.mat'
%
% solthresh: Solidity threshold - nuclei below the threshold are thrown
% out, I typically use solthresh = 0.8, but you can go lower if too many
% nuclei are being removed.
%
% zscale: Scaling factor to resize images in z so that voxels have the
% same length. For our images, 1 pixel is about 1/6 of a micron, so if the
% z-spacing is 1 micron, zscale = 6.
%
%% Outputs
% nuclei: A volumetric label matrix containing the segmented nuclei,
% similar to BWlabel or equivalent.
%
% nucleusData: A structure containing cropped and resized 3D 

%% Function
% Initialize variables
od = cd;
%listNuc = data(MovieNum).ImageFileListNuc;
listNuc = data(MovieNum).ImageFileListCh2;
cd(data(MovieNum).Source)
s4 = strel('sphere',4); s3 = strel('sphere',3); s2 = strel('disk',2);

%% Loop over the specified frames
for t = tRange % each time point
    cd('SegmentationData');
    temp = imread(listNuc{t,1});

    % This section creates two versions of the volumetric nucleus image
    % stack, one roughly filtered and one finely filtered. Could be done
    % without the for loop by using a 3D filter
    iNucRough = zeros(size(temp,1),size(temp,2),size(listNuc,2));
    iNucFine = iNucRough;
    for z = 1:size(listNuc,2) % each plane
        imageNuc = imread(listNuc{t,z});
        imageNuc = double(imageNuc(:,:,1));
        iNucFine(:,:,z) = imgaussfilt(imageNuc,0.5);
        iNucRough(:,:,z) = imgaussfilt(imageNuc,4);
    end
    iNucFine = imresize3(iNucFine,[size(iNucFine,1) size(iNucFine,2) zscale*size(iNucFine,3)]);
    iNucRough = imresize3(iNucRough,[size(iNucRough,1) size(iNucRough,2) zscale*size(iNucRough,3)]);

    % This section pulls the cell segmentation for the current time point,
    % and gets the list of cell indices that will be looped over.
    cframefoldername = sprintf('frame%04d',t);
    cd(cframefoldername);
    imageSeg = load(filename);
    imageSeg = imageSeg.ImageBWlabel;
    imageSeg(imageSeg < 0) = 0;
    cellNumbers = unique(imageSeg(:));
    cellNumbers(cellNumbers == 0) = [];

    %% Loop over cell index numbers
    nuclei = zeros(size(imageSeg,1),size(imageSeg,2),zscale*size(imageSeg,3));
    for c = 1:size(cellNumbers,1) % each cell

        % This section finds the current cell, shrinks it to look at just
        % the cell interior, rescales it in z, then finds its bounding box
        cc = cellNumbers(c,1);
        cCell = imerode(imageSeg == cc,s2);
        cCell = imresize3(cCell,[size(cCell,1) size(cCell,2) zscale*size(cCell,3)]);
        isoCell = false(size(cCell,1),size(cCell,2),size(cCell,3));
        cellProps = regionprops3(cCell,'Centroid','BoundingBox');
        cellBox = floor(cellProps.BoundingBox);

        if size(cellBox,1) > 1
            CC = bwconncomp(cCell);
            numPixels = cellfun(@numel,CC.PixelIdxList);
            [~,idx] = max(numPixels);
            cCell = false(size(cCell));
            cCell(CC.PixelIdxList{idx}) = true;
            cellBox = cellBox(idx,:);
        elseif isempty(cellBox)
            continue
        end

        cCell = cCell(cellBox(2):(cellBox(2)+cellBox(5)),... % get cell chunk
            cellBox(1):(cellBox(1)+cellBox(4)),cellBox(3)+1:(cellBox(3)+cellBox(6)));
        cCell = bwmorph3(cCell,'majority');

        if size(cCell,3) < 3
            continue
        end

        % Then the same chunk is pulled out of both the fine and rough
        % nucleus volume images
        imageNucFineTrim = iNucFine(cellBox(2):(cellBox(2)+cellBox(5)),... % get nuclear chunk
            cellBox(1):(cellBox(1)+cellBox(4)),cellBox(3)+1:(cellBox(3)+cellBox(6)));

        imageNucRoughTrim = iNucRough(cellBox(2):(cellBox(2)+cellBox(5)),... % get nuclear chunk
            cellBox(1):(cellBox(1)+cellBox(4)),cellBox(3)+1:(cellBox(3)+cellBox(6)));

        isoNuc = false(cellBox(5)+1,cellBox(4)+1,cellBox(6));
        cNucXY = isoNuc; cNucXZ = isoNuc; cNucYZ = isoNuc;

        cNuc = imageNucRoughTrim.*cCell;
        cNuc(cNuc == 0) = NaN;

        %% Median threshold the rough nucleus in each direction
        for z = 1:size(cCell,3) % each xy
            cNucZ = squeeze(cNuc(:,:,z));
            cNucMedian = median(cNucZ(:),'omitnan');
            cNucXY(:,:,z)=imbinarize(cNucZ,cNucMedian);
        end

        for y = 1:size(cCell,2) % each xz
            cNucY = squeeze(cNuc(:,y,:));
            cNucMedian = median(cNucY(:),'omitnan');
            if isnan(cNucMedian)
                cNucXZ(:,y,:) = false(size(cNucY,1),size(cNucY,2));
            else
                cNucXZ(:,y,:) = imbinarize(cNucY,cNucMedian);
            end
        end

        for x = 1:size(cCell,1) % each yz
            cNucX = squeeze(cNuc(x,:,:));
            cNucMedian = median(cNucX(:),'omitnan');
            if isnan(cNucMedian)
                cNucYZ(x,:,:) = false(size(cNucX,1),size(cNucX,2));
            else
                cNucYZ(x,:,:) = imbinarize(cNucX,cNucMedian);
            end
        end
        % The rough segmentation is comprised of voxels where all the
        % thresholds agree
        nucSeg = (cNucXY & cNucXZ & cNucYZ);

        CC = bwconncomp(nucSeg);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        if isempty(numPixels)
            continue
        end
        [~,idx] = max(numPixels);
        nucSeg = false(size(nucSeg));
        nucSeg(CC.PixelIdxList{idx}) = true;

        %% Refine segementation
        % This section finds the average image intensities inside and
        % outside the rough segmentation in the fine nucleus image, then
        % creates a new segmentation threshold at the 25th percentile of
        % the linear space between those two values, and segments the fine
        % image at this value.
        cCelSeg = cCell;
        cNucFine = imageNucFineTrim;
        cNucRough = nucSeg;

        outsideImage = cNucFine.*cCelSeg.*(imdilate(cNucRough,s4)-cNucRough);
        outsideImage(outsideImage == 0) = NaN;
        outsideVal = mean(outsideImage(:),'omitnan');

        insideImage = cNucFine.*imdilate(cNucRough,s3);
        insideImage(insideImage == 0) = NaN;
        insideVal = mean(insideImage(:),'omitnan');

        newThresh = prctile(linspace(outsideVal,insideVal),25);
        nucSeg = (cCelSeg.*cNucFine) >= newThresh;

        %% Clean up the segmentation
        % The segmentation is cleaned up using 2D and 3D morphological
        % operations
        for z = 1:size(nucSeg,3)
            nucSeg(:,:,z) = bwmorph(nucSeg(:,:,z),'majority');
            nucSeg(:,:,z) = bwmorph(nucSeg(:,:,z),'fill');
        end

        CC = bwconncomp(nucSeg);
        numPixels = cellfun(@numel,CC.PixelIdxList);
        [~,idx] = max(numPixels);
        nucSeg = false(size(nucSeg));
        if size(CC.PixelIdxList,2) == 0
            continue
        else
            nucSeg(CC.PixelIdxList{idx}) = true;
        end
        nucSeg = logical(smooth3(nucSeg));
        nucSeg = bwmorph3(nucSeg,'majority');

        cObjects = regionprops3(nucSeg,'Volume','VoxelIdxList');
        [~,objInd] = max(cObjects{:,1});

        if isempty(objInd)
            continue
        elseif ~iscell(cObjects.VoxelIdxList)
            continue
        end

        isoNuc(cObjects.VoxelIdxList{objInd,1}) = true;
        isoNuc = bwmorph3(isoNuc,'majority');
        isoNuc = bwmorph3(isoNuc,'clean');
        isoNuc(isoNuc == 1) = cc;

        %% Store segmentation results
        % The solidity of the segmentation is checked, and if the quality
        % is good enough the results are stored.
        solTest = table2array(regionprops3(isoNuc,'Solidity'));
        if solTest <= solthresh
            isoNuc(isoNuc == 1) = 0;
        end

        isoCell(cellBox(2):(cellBox(2)+cellBox(5)),...
            cellBox(1):(cellBox(1)+cellBox(4)),cellBox(3)+1:(cellBox(3)+cellBox(6))) = isoNuc;

        nInds = (isoCell(:) == 1);
        nuclei(nInds) = cc;

        nucProps = regionprops3(isoCell,'Centroid');

        nucleusData(c).ID = cc;
        nucleusData(c).CellCentroid = cellProps.Centroid;
        nucleusData(c).NucCentroid = nucProps.Centroid;
        nucleusData(c).Cell = cCell;
        nucleusData(c).Nuc = isoNuc;
    end
    nucleusData = nucleusData(all(~cellfun(@isempty,struct2cell(nucleusData))));

    save('nuclei','nuclei')
    save('nucleusData','nucleusData')
    clear nucleusData

    cd(od);
end

end

