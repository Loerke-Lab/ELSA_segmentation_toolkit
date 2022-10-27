function [ListCh0,ListCh1] = MakeFileList_LS(path) %,List2
% This function generates an image file list cell array from the TIF files
% in the directory 'path' or in the current directory (default).
%
% INPUT:    path: (optional) directory path where to look for images, if 
%                  this is not specified, the function uses the current
%                  directory by default
%
% OUTPUT:   ListCh0: file list of image channel one
%           ListCh1: file list of image channel two


% check if path was inputted and change directory to that path
if nargin == 1
    cd(path)
end


Mylist = dir('*.tif');

MylistNames = {Mylist.name}';

Nfiles = numel(MylistNames);

idxListCh0 = 1;
FrameVecCh0 = [];
idxListCh1 = 1;
FrameVecCh1 = [];


for ii=1:Nfiles
    NamePartsCells = strsplit(MylistNames{ii},'_');
    Channel = NamePartsCells{3}; % Use 4 for 2colorWT, Use index 3 for Y27, Use index 4 for zip, Use index 4 for nls
    Frame   = NamePartsCells{2}; % Use 6 for 2colorWT, Use index 5 for Y27, Use index 6 for zip, Use index 6 for nls
    
    Frame(1:5) = [];
    FrameInt = str2num(Frame);
    
    %if strcmpi(Channel,'ch0')
    if strcmpi(Channel,'GFP')    
        ListCh0{idxListCh0} = MylistNames{ii};
        idxListCh0 = idxListCh0 + 1;
        FrameVecCh0 = [FrameVecCh0;FrameInt];
    end
    %if strcmpi(Channel,'ch1')
    if strcmpi(Channel,'mCherry')    
        ListCh1{idxListCh1} = MylistNames{ii};
        idxListCh1 = idxListCh1 + 1;
        FrameVecCh1 = [FrameVecCh1;FrameInt];
    end


end

if any(diff(FrameVecCh0)~=1)
    error('Check frame order or missing frames in Channel 0')
end
if any(diff(FrameVecCh0)~=1)
    error('Check frame order or missing frames in Channel 1')
end
ListCh0 = ListCh0';
ListCh1 = ListCh1';

end

%%

function [] = makeSegmentationDataFiles(Nframes)

for t=1:Nframes
    mkdir(sprintf('frame%04d',t));
end
 
end

