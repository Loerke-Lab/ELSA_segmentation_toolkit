function [outputFileList, filepath] = loadImageNames_dualIndex(startpath)
%loadImageNames collects the names of images in a folder with the images
%from an interlaced stack, meaning that the images are numbered along two
%indices, such as time AND z-position
% SUMMARY: [outputFileList] = loadImageNames()
% INPUT:    none, user specifies image names
%           NOTE:   this function is written for a dual index in the image 
%                   names of the general type
%                   e.g. 'ProteinXGFP_time001_zpos001'
%                   where time is the first and zpos the second index
% OUTPUT:   outputFileList  = cell array of complete image names
%                   the results has a double index, where in the
%                   position outputFileList{n,k}
%                   n corresponds to the first index (time in the example
%                   above, and k corresponds to the second index (zpos in
%                   the example above)
%
% first version 09/05/2010 DLoerke


    
%% load first z-stack

if nargin>0
    cd(startpath);
end

% user specifies the first image
[fileNameI1 filePathI1] = uigetfile('.tif','choose first image of first sub-stack (e.g. t000_i000)');
% complete file name of the first image
oriImageName1 = strcat(filePathI1, fileNameI1);


%% load second z-stack for subsequent time point

% user specifies the second image
[fileNameI2 filePathI2] = uigetfile('.tif','choose first image of second sub-stack (e.g. t001_i000)');
% complete file name of the first image
oriImageName2 = strcat(filePathI2, fileNameI2);


%% compile complete file list based on the specified image names
outputFileList = getFileStackNames_two(oriImageName1,oriImageName2);

filepath = filePathI1;


end

