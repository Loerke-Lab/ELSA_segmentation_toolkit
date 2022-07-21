function [ ] = ELSA_42trackCellsInT_SingleZ(data);
% this function tracks cells from one t-frame to the next in 4D image
% stacks; the function assumes that cells have already been tracked in
% z-direction, and uses a specific z-section as reference
% Input:    data: structure that contains lists of image files (as created by the
%                   function ELSA_1loadImageList) and the source (or path) to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%
% Output:   no function output - results are written into the specified
%           directory
%
% 04/14/2011 Dinah Loerke
% 01/06/2013 modified by Tim Vanderleest
% Changes:
% 1. The only input is the data structure that contains the path.
% 2. I had been doing single-Zslice movies in the beginning so I skipped
%    the step of tracking in Z. So this function loads 'ImageBWlabel.mat'
%    instead of the z-tracked output which was named 'ImageMatrixBWlabel_track.mat'
% 3. This function saves it's output as 'ImageMatrixBWlabel_trackT.mat'
%    


%% 

% record original directory (and return to it at the end)
od = cd;


spath = data.Source;


   
cd(spath);
cd('SegmentationData');

% loop over all existing time points
t = 1;
cframefoldername = sprintf('frame%04d',t);
while exist(cframefoldername)==7
    
    cd(cframefoldername);      
    % display current progress of processing
    fprintf('extracting @ timepoint %04d',t);
    
    % ...and upload ImageMatrixBWlabel
    loadmat = load('ImageBWlabel.mat');
    ImageBWlabel_track = loadmat.ImageBWlabel;
    
    if t==1
        % initialize results matrices
        [bwsx,bwsy] = size(ImageBWlabel_track);
        ImageMatrixBWlabel_T = nan*zeros(bwsx,bwsy);
    end
    
    % enter the appropriate z-section for this t into the tracking matrix
    ImageMatrixBWlabel_T(:,:,t) = ImageBWlabel_track;
    
    % update t and folder name
    t=t+1;
    cframefoldername = sprintf('frame%04d',t);
    
    % return to upper directory
    cd(spath);
    cd('SegmentationData');
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for t-loop
fprintf('\n');

% return to start directory   
cd(od); 
% track the subsequent image frame 
ImageMatrixBWlabel_trackAllT = trackArea2( ImageMatrixBWlabel_T );

% after the individual frames are linked for this z, return to all the 
% individual ImageMatrixBWlabel matrices for all the values of t, and
% rename them based on the tracking results

cd(spath);
cd('SegmentationData');
% loop over all existing time points
t = 1;
cframefoldername = sprintf('frame%04d',t);
while exist(cframefoldername)==7
    
    cd(cframefoldername);   
    % display current progress of processing
    fprintf('renaming cells @ timepoint %04d',t);

    
    % the corresponding z-section in the t-tracked matrix is
    ImageBWlabel_trackT = ImageMatrixBWlabel_trackAllT(:,:,t);
    

    % save results into the current results folder
    save('ImageBWlabel_trackT','ImageBWlabel_trackT');
    cd(od);    
    
    % update t and folder name
    t=t+1;
    cframefoldername = sprintf('frame%04d',t);
    % return to upper folder
    cd(spath);
    cd('SegmentationData');
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');

end % of for t-loop


fprintf('\n');
cd(od);
    
end % of function



