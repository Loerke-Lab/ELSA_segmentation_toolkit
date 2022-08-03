function [imageFileList] = ELSA_1loadImageList(path, switchvar)
% this function reads all images in a 4-D image stack, outputs them as a
% filelist for later reference, and prepares the appropriate folder
% Input:    path: (optional) directory path where to look for images, if 
%                  this is not specified, the function uses the current
%                  directory by default
%           switchvar: the function assumes a double index in the image 
%                  file names, and the default parsing is based on names of
%                  the type 'xxx_taaa_fbbb', i.e. with the time index 
%                  first and the section index last. If the indices are 
%                  switched (or you as a user wish to switch them), 
%                  indicate this by setting switchvar=1
% Output:   imageFileList: cell array of complete image names
%                   the results has a double index, where in the
%                   position outputFileList{n,k}
%                   n corresponds to the first index (time in the example
%                   above, and k corresponds to the second index (zpos in
%                   the example above)
%
% 03/05/2011 Dinah Loerke

%% PREPARATION AND INITIALIZATION

% specify the default mode to be 0
mode = 0;

% check if there was more than one input in the function
if nargin>1
    % if switchvar was an input, check if it was inputted as 1 and change the mode to 1 instead of 0
    if ~isempty(switchvar)
        if switchvar==1
            mode = 1;
        end
    end
end

% record original directory (to return to it at the end)
od = cd;

% specify the default searchpath to be the original directory the function is run from 
searchpath = od;

% if path was an input, change the searchpath to be the inputted path
if nargin>0
    searchpath = path;
end
   
% read list of raw image files in the stack with double index; the
% assumption is that images are named sequentially along the lines of 
% 'xxxxt001_f001.tif' with a numerical index for time t and for section f
[imageFileList, filepath] = loadImageNames_dualIndex(searchpath);

% create a folder for each time point in the Segmentation Data folder
numframes = size(imageFileList,1);
numsections = size(imageFileList,2);

% invert list (rows vs columns) if necessary (if switchvar is set to 1)
if mode==1
    for i1 = 1:numframes
        for i2 = 1:numsections
            centry = imageFileList{i1,i2};
            imageFileList_copy{i2,i1} = centry;
        end
    end
    imageFileList = imageFileList_copy;
    numframes = size(imageFileList,1);
    numsections = size(imageFileList,2);
end

% change directory to the filepath (the folder where images are stored)
cd(filepath);

% if a segmentation results folder doesn't already exist, create one
if  ~ (exist('SegmentationData')==7) 
    [success, msg, msgID] = mkdir(filepath,'SegmentationData');
end

% move to segmentation data folder
cd('SegmentationData');

% record the current directory under SegmentationData
od2 = cd;

% DEFAULT time and section numbers are set to min and max available,
% unless different input is specified
tstart = 1;
tend = numframes;
% if nargin>0
%     if ~isempty(tstartend)
%         tstart = tstartend(1);
%         tend = tstartend(2);
%     end
% end

fstart = 1;
fend = numsections;
% if nargin>1
%     if ~isempty(fstartend)
%         fstart = fstartend(1);
%         fend = fstartend(2);
%     end
% end


% loop over desired time points and create subfolders in SegmentationData
for t=tstart:tend
    cframefoldername = sprintf('frame%04d',t);
    if  ~ (exist(cframefoldername)==7) 
        [success, msg, msgID] = mkdir(od2,cframefoldername);
    end
end

% return to the original directory
cd(od);

end % of function