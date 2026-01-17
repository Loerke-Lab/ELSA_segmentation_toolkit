function [] = matchNucleiAndCellIDs(data)

% INPUT:    data = data structure minimally containing fields Source, 
%                   ImageFileListNuc, and PixRes. This implementations
%                   assumes the file convention associated with lattice
%                   light sheet data, where each time point is saved as a
%                   volumetric image.

% OUTPUT:   relabeled nuclei matrix with nuclei numbers that match cell numbers

% choose a well segmented cell layer
layer = 15;

% get number of frames to loop over
cd(data.Source)
nFrames = numel(dir('SegmentationData'))-2;
cd('SegmentationData');

% loop through time points
for t = 1:nFrames

    % move to frame folder
    cframefoldername = sprintf('frame%04d',t); % naming convention of seg folders
    cd(cframefoldername);

    % load the cell and nuclei segmentation matrices
    nuclei = load('nuclei_cp.mat').nuclei_cp;

    cells = load('ImageBWlabel.mat').ImageBWlabel;
    cellsMax = max(cells,[],3);

    cellsLayer = cells(:,:,layer);
    

    if t > 1
        % if starting at frame > t=1, load previous frame for tracking
        cframefoldername_last = sprintf('frame%04d',t-1); % naming convention of seg folders
        cd(data.Source); cd('SegmentationData'); cd(cframefoldername_last);
        % nuclei_last = load('nuclei.mat').nuclei;
        nuclei_last = load('nuclei_cp.mat').nuclei_cp;
        nuclei_last = max(nuclei_last,[],3);

        nuc_max_proj = max(nuclei, [], 3); % max projection to get list of tracking IDs
        ind_list = unique(nuc_max_proj(:)); ind_list(ind_list==0) = [];
        counter = max(ind_list) + 1; % updating variable for tracking IDs
        for n = 1:numel(ind_list)
            fprintf(' nuc #%04d/%04d', n, length(ind_list)); % command line message
            % find pixels in last frame that overlap with nuc n in current frame
             track_list = nuclei_last(nuc_max_proj == ind_list(n)); 
             % find IDs with highest overlap between frames
             ID = mode(track_list(:));
             if ID == 0
                 ID = counter;
                 counter = counter + 1;
             end
             % re-label nuclei ID #'s to match
             nuclei(nuclei==ind_list(n)) = ID;
             fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); % clear command line message
        end
    end

    cd ..
    cd(cframefoldername)

    % display current progress of processing
    fprintf('extracting @ timepoint %04d\n',t);

    nucleiMax = max(nuclei,[],3);
    nucleiLayer = nuclei(:,:,layer);

    % figure(1)
    % vislabels(nucleiMax)
    % figure(2)
    % vislabels(cellsLayer)
    % pause

    % indList = unique(cells(:)); 
    % indList(indList==0) = [];

    % get the list of unique nuclei numbers
    indList = unique(nucleiMax(:)); 
    indList(indList==0) = [];

    % loop through nuclei numbers
    for i = 1:numel(indList) 

        % this line in case nuclei are not always ascending by one
        ii = indList(i);

        % find the cell number where the nuclei index is
        %cellNum = cellsMax(nucleiMax==ii);
        %cellNum = cellsLayer(nucleiLayer==ii);
        cellNum = cellsLayer(nucleiMax==ii);

        % if no cell is found, remove nucleus, else change nucleus number
        % to match cell number
        if sum(cellNum(:)) == 0
            nuclei(nuclei==ii) = 0;
        else
            nuclei(nuclei==ii) = mode(cellNum);
        end


    end

    nucleiMax = max(nuclei,[],3);

    % if t == 663
    %     figure(1)
    %     vislabels(nucleiMax)
    %     figure(2)
    %     vislabels(cellsLayer)
    %     pause
    % end

    % % save the new labeled nuclei
    nucleiCPlabel = nuclei;
    save('nucleiCPlabel.mat', 'nucleiCPlabel', '-v7.3');

    % move back a directory 
    cd ..

    % delete message 'extracting @ timepoint %04d' before next loop:
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b'); 

end