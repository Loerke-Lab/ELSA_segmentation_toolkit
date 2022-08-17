function ELSA_initializeEmptySegData_LS(data,MovieNum,tvec)
%ELSA_initializeEmptySegData_LS goes into the SegmentationData directory
%frame folders and initializes the seeds, mask, ImageSegment, and
%ImageBWlabel 3D matrices. It first loads an image stack to determine the
%size of the matrices to create


% get original directory
od = cd;


list = data(MovieNum).ImageFileList;


% load zstack
cd(data(MovieNum).Source)
img = tiffreadVolume(list{1});
% for z = 1:size(list,2)
%     img(:,:,z) = imread(list{z,1});
% end

for t=tvec

    % if seeds already exist load them, else initialize empty seeds and mask
    cd(data(MovieNum).Source)
    cd('SegmentationData');
    cframefoldername = sprintf('frame%04d',t);
    cd(cframefoldername);
    if isempty(dir('seeds.mat'))
        seeds = false(size(img));
        mask  = false(size(img));
        ImageSegment = false(size(img));
        ImageBWlabel = zeros(size(img));
        
        save('seeds','seeds')
        save('mask','mask')
        save('ImageSegment', 'ImageSegment');
        save('ImageBWlabel', 'ImageBWlabel');


    end


    
    
end

cd(od)
end

