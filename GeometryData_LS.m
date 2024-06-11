function [] = GeometryData_LS(data,MovieNum,Nframes)
%GEOMETRYDATA gets geometric properties of each cell for each frame of a
%movie and saves them to an array

% INPUT:    data: structure that contains the fields for ImageFileList and
%                 Source.
%           MovieNum: index of the data structure because I store
%                     multiple movies in a data structure.
%           Nframes: number of frames to be analyzed
%
% OUTPUT:   No output but results are saved to the parent directory of
%               the data.Source path. The saved mat file, named GeometryData,
%               contains both an Area and Perimeter array where rows 
%               correspond to cells and columns correspond to time points. 
%               The vector CellNumbers gives the cell ID for each row.


%trackingMatrixFilename = 'ImageBWlabel_trackT.mat';
trackingMatrixFilename = 'ImageBWlabel.mat';

% record original directory (to return to at the end)
od = cd;


% load one seeds matrix to get the number of z-layers
cd(data(MovieNum).Source)
cd('SegmentationData/')
cd('frame0001')
seeds = load('seeds.mat').seeds;

Nlayers = size(seeds,3);

% initialize array for cell numbers
CellNumbers = [];

% this time loop is to get all the cell numbers
for t=1:Nframes
    
    % display current progress of processing
    fprintf('Getting cell numbers @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    % move to the frame folder location
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);

    % load the labeled tracking matrix
    loadTmatrix = load(trackingMatrixFilename);
    %Tmatrix = loadTmatrix.ImageBWlabel_trackT;
    Tmatrix = loadTmatrix.ImageBWlabel;
    
    % Find the unique cell numbers in Tmatrix
    CellNumberst = unique(Tmatrix(:));

    % store unique cell numbers from each loop
    CellNumbers = unique([CellNumbers; CellNumberst]);
    
    % change directory up one level
    cd ..

    % delete the message 'Getting cell numbers @ timepoint %04d' before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
end

% negative values and zeros do not correspond to cell numbers
CellNumbers(CellNumbers < 1) = [];

% get the total number of unique cells in the movie
Ncells = length(CellNumbers);


% initialize Geometric parameter array
Area = NaN(Ncells,Nframes,Nlayers);
Perimeter = NaN(Ncells,Nframes,Nlayers);




%% Loop through every frame in the movie and store geometric properties for each cell 
for t=1:Nframes
    
    % display current progress of processing
    fprintf('Processing @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    % move to the frame folder location
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);

    % load the labeled tracking matrix
    %Tmatrix = load(trackingMatrixFilename).ImageBWlabel_trackT;    
    Tmatrix = load(trackingMatrixFilename).ImageBWlabel;    
    
    % Now loop over all the z-layers to process
    for z=1:Nlayers 
        
        % set tracking matrix for the z layer
        TmatrixZ = Tmatrix(:,:,z);
    
        % get region props
        STATS = regionprops(TmatrixZ,'Area','Perimeter');
        A = [STATS.Area];
        P = [STATS.Perimeter];

        % first gather all cell labels
        CellNumberst = unique(TmatrixZ(:));
        CellNumberst(CellNumberst < 1) = [];

        % determine row indices of the geometry data array for the entire movie
        idx = find(ismember(CellNumbers,CellNumberst));

        % store the area and perimeter data
        Area(idx,t,z) = A(CellNumberst);
        Perimeter(idx,t,z) = P(CellNumberst);

    end
    
    % change directory up one level
    cd ..

    % delete the message 'Processing @ timepoint %04d' before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    
end

% find all cells that are nan in every z layer
allz = all(isnan(Area),3);

% find rows where all cells are nan
allr = all(allz,2);

% set those elements to empty in Area and Perimeter arrays
Area(allr,:,:) = [];
Perimeter(allr,:,:) = [];

% change back to source directory and up one directory to save data with
% all other data arrays
cd(data(MovieNum).Source)
cd ..
save('GeometryData','Area','Perimeter','CellNumbers')

% return to original directory
cd(od)

end


