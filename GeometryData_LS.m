function [] = GeometryData_LS(data,MovieNum,Nframes)
%GEOMETRYDATA gets geometric properties of each cell for each frame of a
%movie and saves them to an array

% INPUT:        data: structure that contains the fields for ImageFileList and
%                     Source.
%               MovieNum: index of the data structure because I store
%               multiple movies in a data structure.
%trackingMatrixFilename = 'ImageBWlabel_trackT.mat';
trackingMatrixFilename = 'ImageBWlabel.mat';

od = cd;


% get the total number of unique cells in the movie

%Nframes = numel(data(MovieNum).ImageFileListGap43);


% load one seeds matrix to get the number of z-layers
cd(data(MovieNum).Source)
cd('SegmentationData/')
cd('frame0001')
seeds = load('seeds.mat').seeds;

Nlayers = size(seeds,3);


CellNumbers = [];
% this time loop is to get all the cell numbers
for t=1:Nframes
    
    fprintf('Getting cell numbers @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);
    % load the labeled tracking matrix
    loadTmatrix = load(trackingMatrixFilename);
    Tmatrix = loadTmatrix.ImageBWlabel; 
    
    CellNumberst = unique(Tmatrix(:));
    
    CellNumbers = unique([CellNumbers; CellNumberst]);
    
    cd ..
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
    
    fprintf('Processing @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);
    % load the labeled tracking matrix
    Tmatrix = load(trackingMatrixFilename).ImageBWlabel;    
    
    % Now loop over all the z-layers to process
    for z=1:Nlayers 
        
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

        Area(idx,t,z) = A(CellNumberst);
        Perimeter(idx,t,z) = P(CellNumberst);

    end
    
    
    cd ..
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    
end

allz = all(isnan(Area),3);
allr = all(allz,2);

Area(allr,:,:) = [];
Perimeter(allr,:,:) = [];

% change back to source directory and up one directory to save data with
% all other data arrays
cd(data(MovieNum).Source)
save('GeometryData','Area','Perimeter','CellNumbers')

cd(od)

end


