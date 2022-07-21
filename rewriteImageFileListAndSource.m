function [dataNew] = rewriteImageFileListAndSource(data,oldString,newString)

% this function rewrites the ImageFileList and Source strings in a data structure

% INPUT:    data: a data structure that contains the ImageFileList and Source fields
%           oldString: the old portion of the string you want to replace (should be in single quotation marks)
%           newString: the new portion of the string you want to replace (should be in single quotation marks)

% OUTPUT:   dataNew: a new data structure identical to the original data
%                    structure but the ImageFileList and Source have been changed


% find the number of rows to loop over
numRows = size(data,2);

% make a copy of the data structure
dataNew = data;

% loop over the rows to change
for a = 1:numRows

    % set the ImageFileList variable
    ImageFileListToChange = dataNew(a).ImageFileList;
    % change the string in ImageFileList
    newFileList = strrep(ImageFileListToChange,oldString,newString);
    % save the new string in the new data structure
    dataNew(a).ImageFileList = newFileList;
    
    %%% UNCOMMENT if the movies have two channels
%     % set the ImageFileListCh2 variable
%     ImageFileListCh2ToChange = dataNew(a).ImageFileListCh2;
%     % change the string in ImageFileList
%     newFileListCh2 = strrep(ImageFileListCh2ToChange,oldString,newString);
%     % save the new string in the new data structure
%     dataNew(a).ImageFileListCh2 = newFileListCh2;
    
    % set the Source variable
    SourceToChange = dataNew(a).Source;
    % change the string in Source
    newSource = strrep(SourceToChange,oldString,newString);
    % save the new string in the new data structure
    dataNew(a).Source = newSource;
    
end
    