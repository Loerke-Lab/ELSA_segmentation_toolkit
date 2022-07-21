function [outputFileList]=getFileStackNames_two(firstfilename,secondfilename)
% getFileStackNames_two returns a cell array containing all file names 
% with path) belonging to an interlaced stack (i.e. with two indices, such
% as t=time and z=cross section
%
% SYNOPSIS
% [outputFileList]=getFileStackNames_two(firstfilename,secondfilename)
%
% INPUT     firstfilename: name of the first greyvalue image to be read
%                    including the full path
%                    the actual filename must consist of 
%                    - alphanumeric body
%                    - numeric number
%                    - extension
%           secondfilename: as above
%
% OUTPUT   outputFileList: names of all files belonging to the stack, in a
%               cell array of complete image names with a double index,
%               such that in position
%               outputFileList{n,k}
%               n corresponds to the first index (time in the example
%               above, and k corresponds to the second index (zpos in
%               the example above)
%                          
%
% DEPENDENCES
%
% Aaron Ponti, October 4th, 2002
% modified Dinah Loerke, Sept 2010

oldDir = [];

% Output
outputFileList = {};

% first file name segmentation
[fpath1,fname1,fno1,fext1]=getFilenameBody(firstfilename);

if(isempty(fname1) | isempty(fno1) | isempty(fext1) )
   error('invalid first filename specified');
end;

% second file name segmentation
[fpath2,fname2,fno2,fext2]=getFilenameBody(secondfilename);

if(isempty(fname2) | isempty(fno2) | isempty(fext2) )
   error('invalid second filename specified');
end;

[body,no1,no2,bridge]=getFilenameBody_two(fname1,fname2);

if strcmp(fpath1,fpath2)
    fpath = fpath1;
else
    error('paths don''t match');
end


if strcmp(fext1,fext1)
    fext = fext1;
else
    error('extensions don''t match');
end



if(~isempty(fpath))
	% change to stack directory
   oldDir = cd(fpath);
else
   %check if it is in the matlab search path
   tempName=which(firstfilename);
   if(~isempty(tempName))
      [fpath,fname,fno1,fext]=getFilenameBody(tempName);
      oldDir = cd(fpath);
	end;
end;

dirListing = dir;

% get all relevant filenames in the specified directory
iEntry = 1;
fileList = {};
for( i = 1:length(dirListing))
   if(~dirListing(i).isdir)
      fileList(iEntry) = lower({dirListing(i).name});
      iEntry = iEntry + 1;
   end;
end;

% initialize length of numeric body part
l_fno=length(num2str(fno1));
l_no=length(num2str(no1));

if(~isempty(fileList))
    
    % initialize stump index
    tEntries=1;
    % initialize stump name number index
    imIndxt = str2num(no1);
    % initialize stump search name
    searchName_trunc    = [ body, num2str(imIndxt,['%.' num2str(l_no) 'd']), bridge ];
    
    % loop over all available name stumps
    while( ~isempty(strmatch(lower(searchName_trunc),fileList)))
                
        % initialize full name index      
        nEntries=1;
        % initialize full name number index 
        imIndx = str2num(fno1);
        % initialize full search name
        searchName= [ body, num2str(imIndxt,['%.' num2str(l_no) 'd']), bridge,...
                        num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];
        
        % loop over all available full names available for this stump
        while( ~isempty(strmatch(lower(searchName),fileList)))
            
            % enter data into structure
            outputFileList(tEntries,nEntries)={strcat(fpath,filesep,searchName)};
            
            % control vector
            cvec = [tEntries,nEntries];
            
            % increase name index
            nEntries = nEntries + 1;
            index(nEntries) = imIndx;
            % increase full name number index
            imIndx = imIndx + 1;
            % update search name        
            searchName= [ body, num2str(imIndxt,['%.' num2str(l_no) 'd']), bridge,...
                        num2str(imIndx,['%.' num2str(l_fno) 'd']),fext];            
            
        end % of while
        
        % increase stump index
        imIndxt = imIndxt + 1;
        % increase stump 'name' number index
        tEntries=tEntries+1;
        % update stump name
        searchName_trunc = [ body, num2str(imIndxt,['%.' num2str(l_no) 'd']), bridge ];
    
    end % of while-loop
    
end % of if

% change back to original directory
if(~isempty(oldDir))
   cd(oldDir);
end

end % of function

