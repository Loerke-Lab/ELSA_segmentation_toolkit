function [] = GeometryData(data,MovieNum)
% GEOMETRYDATA gets geometric properties of each cell, namely Area and 
% Perimeter, for the entire movie and saves them to an array

% INPUT:        data: structure that contains the fields for ImageFileList and
%                     Source.
%               MovieNum: index of the data structure because I store
%               multiple movies in a data structure.
%
%
% OUTPUT:       No output but results are saved to the parent directory of
%               the data.Source path. The saved mat file, named GeometryData,
%               contains both an Area and Perimeter array where rows 
%               correspond to cells and columns correspond to time points. 
%               The vector CellNumbers gives the cell ID for each row.
%
% T. Vanderleest 5/13/2013

% record original directory (to return to at the end)
od = cd;

%% first loop through each frame of the movie and find all unique cell IDs
% tracking numbers so that the array can be initialized
list = data(MovieNum).ImageFileList;
Nframes = length(list);
CellNumbers = [];

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
    loadTmatrix = load('ImageBWlabel_trackT.mat');
    Tmatrix = loadTmatrix.ImageBWlabel_trackT; 
    
    % Find the unique cell numbers in Tmatrix
    CellNumbersT = unique(Tmatrix(:));
    
    % store unique cell numbers from each loop
    CellNumbers = unique([CellNumbers; CellNumbersT]);
    
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
Area = NaN(Ncells,Nframes);
Perimeter = NaN(Ncells,Nframes);

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
    loadTmatrix = load('ImageBWlabel_trackT.mat');
    Tmatrix = loadTmatrix.ImageBWlabel_trackT;
    Tmatrix(Tmatrix < 1) = 0;
    
    % get region props
    STATS = regionprops(Tmatrix,'Area','Perimeter');
    A = [STATS.Area];
    P = [STATS.Perimeter];
    
    % first gather all cell labels
    CellNumbersT = unique(Tmatrix(:));
    CellNumbersT(CellNumbersT == 0) = [];
    
    % determine row indices of the geometry data array for the entire movie
    idx = find(ismember(CellNumbers,CellNumbersT));
    
    % store the area and perimeter data
    Area(idx,t) = A(CellNumbersT);
    Perimeter(idx,t) = P(CellNumbersT);
    
    % change directory up one level
    cd ..

    % delete the message 'Processing @ timepoint %04d' before the next loop
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b')
    
end

% change back to source directory and up one directory to save data with
% all other movie specific data
cd(data(MovieNum).Source)
cd ..
save('GeometryData','Area','Perimeter','CellNumbers')

% return to original directory
cd(od)

end

