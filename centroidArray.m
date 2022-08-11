function [] = centroidArray(data,MovieNum)
% centroidArray gathers the centroid positions of all the cells for all time
% points and stores them in an array. It gets the centroids from 
% dstruct_cellCentroids.mat. Rather than output results it just saves them
% to the parent directory of data.Source.
%
% INPUT:        data: structure that contains the fields for ImageFileList and
%                     Source.
%               MovieNum: index of the data structure because I store
%               multiple movies in a data structure.

% OUTPUT:       No output but results are saved to the parent directory of
%               the data.Source path. The saved mat file contains an array
%               called CentroidArray which contains both the X and Y
%               coordinates of the cells. Each row corresponds to a cell
%               and every two columns correspond to X and Y for a time
%               point, to extract X and Y use the following:
%                   Cy = CentroidArray(:,1:2:end);
%                   Cx = CentroidArray(:,2:2:end); 
%               The vector CellNumbers gives the cell ID for each row.
%
% T. Vanderleest 8/7/2014

% record original directory (to return to at the end)
od = cd;

% go to segmentation data directory where all the frame folders are
cd(data(MovieNum).Source)
cd('SegmentationData/')

%% first get all unique cell ID numbers
% initialize results matrix
idx_unique = [];

% enter a loop through all frame folders and load cell centroid matrices
t = 1;
cframefoldername = sprintf('frame%04d',t);
while exist(cframefoldername)==7

    % change directory to specific frame folder
    cd(cframefoldername);   

    % display current progress of processing
    fprintf('extracting @ timepoint %04d',t);
    
    % load cell centroid data
    loaddata = load('dstruct_cellCentroids.mat');

    % centroid position array (x,y) where the row # is the cell ID #
    positions = loaddata.dstruct_cellCentroids.positions;
    
    % get one column of positions 
    % find which ones are finite
    % the finite indices correspond to actual cell tracking numbers
    pos = positions(:,1); 
    idx_t = find(isfinite(pos));

    % save unique cell ID numbers for all of the loop
    idx_unique = unique([idx_unique;idx_t]);
    
    % update t and folder name
    t=t+1;
    cframefoldername = sprintf('frame%04d',t);
    
    % return to upper directory
    cd ..

    % delete the message 'extracting @ timepoint %04d' before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');

end

% the total number of frames will be t - 1
Nframes = t - 1;

% save cell ID numbers in a matrix
CellNumbers = idx_unique;


%% Initialize a centroidArray with a row for each unique cell tracking number
CentroidArray = NaN(length(idx_unique),Nframes*2);

% now loop through each frame and store centroid positions in array
for t=1:Nframes
    
    % change directory to specific frame folder
    cframefoldername = sprintf('frame%04d',t);
    cd(cframefoldername);      

    % display current progress of processing
    fprintf('extracting centroids @ timepoint %04d',t);
    
    % load cell centroid data
    loaddata = load('dstruct_cellCentroids.mat');

    % centroid position array (x,y) where the row # is the cell ID #
    positions = loaddata.dstruct_cellCentroids.positions;
    
    % find which ones are finite
    idx = find(isfinite(positions(:,1)));
    

    idxArray = ismember(CellNumbers,idx);

    % save positions in centroid array
    CentroidArray(idxArray,2*t+[-1 0]) = positions(idx,:);
    
    % return to upper directory
    cd ..

    % delete the message 'extracting centroids @ timepoint %04d' before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
end
    
 
% save the results to the parent directory of the raw image source directory
cd(data(MovieNum).Source)
cd ..
save('CentroidArray','CellNumbers','CentroidArray')
    
% return to original directory
cd(od)

end

