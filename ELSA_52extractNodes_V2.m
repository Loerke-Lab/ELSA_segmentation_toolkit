function [ ] = ELSA_52extractNodes_V2(data, tvec)
% this function segments the adhesions between cells and the connecting
% nodes, using the tracking information
%
% Input:    data: structure that contains lists of image files (as created by the
%                   function ELSA_1loadImageList) and the source (or path) to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%           tvec (optional): timevec specifying the frames to be analyzed
%              (e.g. [1:20])(optional) 
%
% Output:   no function output - results are written into the specified
%           directory
%
% 03/05/2011 Dinah Loerke
% 01/06/2013 minor modifications, Tim Vanderleest
% Changes:
% 1. Changed the input to the data structure which contains the image file
%   list
% 2. I have been calling my tracked cells file 'ImageBWlabel_trackT.mat' so
%   that doesn't need to be an input

% Note: this function still can handle volumetric images  

%% 
% record the file name
filename = 'ImageBWlabel_trackT.mat';

% record original directory (to return to it at the end)
od = cd;

% collect information from the data structure
list = data.ImageFileList;
spath = data.Source;


% create a folder for each time point in the Segmentation Data folder
numframes = size(list,1);
numsections = size(list,2);


% DEFAULT time and section numbers are set to min and max available,
% unless different input is specified
tstart = 1;
tend = numframes;

fstart = 1;
fend = numsections;

tlvec = [tstart:1:tend];
if nargin>1
    if ~isempty(tvec)
        tlvec = tvec;
    end
end


% loop over all timepoints
for ti = 1:length(tlvec)
    
    % specify t inside the loop
    t = tlvec(ti);
    
    % display current progress of processing
    fprintf('node extraction @ timepoint %04d',t);
    
    % move to appropiate timepoint subfolder...
    cd(spath);
    cd('SegmentationData');
    cframefoldername = sprintf('frame%04d',t);   
    cd(cframefoldername);
    
    % upload ImageMatrixSegment 
    loadmat = load('ImageSegment.mat');
    ImageMatrixSegment = loadmat.ImageSegment;

    % ...and upload ImageMatrixBWlabel_track if it exists
    if exist(filename)==2
        loadmat = load(filename);
        ImageMatrixBWlabel_use = loadmat.ImageBWlabel_trackT;
    end
    
    % return to original directory
    cd(od);
    
    % extract position centroids from BWlabel
    for z = fstart:fend

        % current image
        cimage = ImageMatrixBWlabel_use(:,:,z);
        
        % create a copy and set regions where cimage is -1 to 0 instead
        cimage_copy = cimage;
        cimage_copy(cimage==-1) = 0;

        % extract centroid positions of remaining separate areas
        rprops = regionprops(cimage_copy, 'Centroid');
        centroids = cat(1, rprops.Centroid);
        dstruct_cellCentroids(z).positions = centroids;   
        
    end
   
    
    % based on the segmentation results, we now determine nodes in the grid
    % and their cell-neighbor relationships; this function also 
    % automatically saves the results into the current folder
    [dstruct_nodes,dstruct_nodeNodeMat,dstruct_cellCellMat,dstruct_nodeCellMat] = nodeAnalysis(ImageMatrixSegment,dstruct_cellCentroids,ImageMatrixBWlabel_use);
    
    % move to appropiate timepoint subfolder
    cd(spath);
    cd('SegmentationData');
    cd(cframefoldername);
    
    % save results into the current results folder
    save('dstruct_cellCentroids', 'dstruct_cellCentroids');
    save('dstruct_nodes', 'dstruct_nodes');
    save('dstruct_nodeNodeMat', 'dstruct_nodeNodeMat');
    save('dstruct_cellCellMat', 'dstruct_cellCellMat');
    save('dstruct_nodeCellMat', 'dstruct_nodeCellMat');
    
    % return to original directory
    cd(od);    
       
    % delete message 'node extraction @ timepoint %04d' before next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for t

% enter new line
fprintf('\n');

    
end % of function