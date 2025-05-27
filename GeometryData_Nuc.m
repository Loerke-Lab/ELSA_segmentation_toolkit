function [] = GeometryData_Nuc(data,MovieNum)%,cZ)
%GEOMETRYDATA gets geometric properties of each cell for each frame of a
%movie and saves them to an array

% INPUT:        data: structure that contains the fields for ImageFileList and
%                     Source.
%               MovieNum: index of the data structure because I store
%               multiple movies in a data structure.

od = cd;

%% first loop through each frame of the movie and find all unique cell
% tracking numbers so that the array can be initialized
% list = data(MovieNum).ImageFileList;
list = data(MovieNum).NewImageListCh0;
%Nframes = length(list);
Nframes = 750;
NucNumbers = [];
for t=1:Nframes
    
    fprintf('Getting cell numbers @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);
    % load the labeled tracking matrix
    %loadTmatrix = load('ImageBWlabel_trackT.mat');
    %Tmatrix = loadTmatrix.ImageBWlabel_trackT;
    loadTmatrix = load('ImageBWlabel.mat');
    Tmatrix = loadTmatrix.ImageBWlabel;


%     Tmatrix = Tmatrix(:,:,cZ);
    
    CellNumberst = unique(Tmatrix(:));
    
    NucNumbers = unique([NucNumbers; CellNumberst]);
    
    cd ..
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');
end
% negative values and zeros do not correspond to cell numbers
NucNumbers(NucNumbers < 1) = [];

% get the total number of unique cells in the movie
Ncells = length(NucNumbers);
Nlayers = size(Tmatrix,3);

% initialize Geometric parameter array
Area = NaN(Ncells,Nframes,Nlayers);
Perimeter = NaN(Ncells,Nframes,Nlayers);
% Intensity = NaN(Ncells,Nframes);
% RawIntensity = NaN(Ncells,Nframes);

%% Loop through every frame in the movie and store geometric properties for each cell 
for t=1:Nframes
    
    fprintf('Processing @ timepoint %04d',t);
    
    % load tracking matrix for given time frame
    cd(data(MovieNum).Source);
    cd('SegmentationData');
    framefoldername = sprintf('frame%04d',t);
    cd(framefoldername);
    % load the labeled tracking matrix
    loadTmatrix = load('nuclei.mat');
    Nmatrix = loadTmatrix.nuclei;
    
    for cZ = 1:size(Nmatrix,3)
    Tmatrix = Nmatrix(:,:,cZ);
    
    Tmatrix(Tmatrix < 1) = 0;
    
%     Imatrix = imread(data.ImageFileListNuc{t,cZ});
    
    % get region props
    STATS = regionprops(Tmatrix,'Area','Perimeter');
    A = [STATS.Area];
    P = [STATS.Perimeter];
    
    % first gather all cell labels
    CellNumberst = unique(Tmatrix(:));
    CellNumberst(CellNumberst == 0) = [];
    
%     statvec = 1:max(CellNumberst);
%     statidx = find(ismember(statvec,CellNumberst));
    
    % determine row indices of the geometry data array for the entire movie
    idx = find(ismember(NucNumbers,CellNumberst));
    
%     for l = 1:size(statvec,2)
%         temp = (Tmatrix == statvec(l));
%         temp = Imatrix(temp);
%         I(1,l) = mean(temp);
%         RI(1,l) = sum(temp);
%     end
    
    Area(idx,t,cZ) = A(CellNumberst);
    Perimeter(idx,t,cZ) = P(CellNumberst);
    end
%     Intensity(idx,t) = I(CellNumberst);
%     RawIntensity(idx,t) = RI(CellNumberst);
    
    
    cd ..
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    
end

% change back to source directory and up one directory to save data with
% all other data arrays
cd(data(MovieNum).Source)
%cd ..
% NucAreaPhase = AreaMatrixToPhaseMatrix(Area,1/15);
% save('GeometryDataNuc','Area','Perimeter','NucAreaPhase','NucNumbers')
save('GeometryDataNuc','Area','Perimeter','NucNumbers')


cd(od)

end

