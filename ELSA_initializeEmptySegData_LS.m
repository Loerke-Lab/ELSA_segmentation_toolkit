function ELSA_initializeEmptySegData_LS(data,MovieNum,tvec)
%ELSA_initializeEmptySegData_LS goes into the SegmentationData directory
%frame folders and initializes the seeds, mask, ImageSegment, and
%ImageBWlabel 3D matrices. It first loads an image stack to determine the
%size of the matrices to create.

% INPUT:    data: structure that contains lists of image files and the source to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%           MovieNum: index of the data structure because
%                     multiple movies can be stored in a data structure.
%           tvec: timevec specifying the frames to be analyzed (e.g. [1:20])
%
%
% OUTPUT:   no function output - results are written into the specified
%           directory



% get original directory (to return to at the end)
od = cd;


% pull the image file list from data structure
list = data(MovieNum).ImageFileList;


% load zstack
cd(data(MovieNum).Source)
img = tiffreadVolume(list{1});

% this for loop is for 3D movies, if the movie is only 2D this can be commented out
for z = 1:size(list,2)
    img(:,:,z) = imread(list{z,1});
end


% loop over time points
for t=tvec
    
    % move to the time point segmentation folder
    cd(data(MovieNum).Source)
    cd('SegmentationData');
    cframefoldername = sprintf('frame%04d',t);
    cd(cframefoldername);

    % if seeds already exist load them, else initialize empty seeds and mask
    % initialize empty ImageSegment and ImageBWlable arrays
    if isempty(dir('seeds.mat'))
        seeds = false(size(img));
        mask  = false(size(img));
        ImageSegment = false(size(img));
        ImageBWlabel = zeros(size(img));
        

        % save seeds, mask, ImageSegment, ImageBWlabel arrays
        save('seeds','seeds')
        save('mask','mask')
        save('ImageSegment', 'ImageSegment');
        save('ImageBWlabel', 'ImageBWlabel');


    end

  
end

% return to original directory
cd(od)

end

