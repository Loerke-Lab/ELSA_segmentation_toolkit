function [body,no1,no2,bridge]=getFilenameBody_two(fname1,fname2)
% GETFILENAMEBODY_TWO splits a filename into its body, number and 
% extension part
%
% SYNOPSIS [path,body,no,ext]=getFilenameBody(fname)
% 
% INPUT fname : filename; the following filename structure must be 
%                         preserved:
%                         - alphanumeric body
%                         - number before extension
%                         - extension separated
%       separator: (opt) string preceding the number, e.g. '_'. Default: ''
%
% OUTPUT path : string with the path, [] if non-existent
%        body : string with body, [] if non-existent
%        no   : string with the number, [] if non-existent
%        ext  : extension, [] if non-existent
%
% SAMPLE getFileNameBody('test_1.tif') returns
%        []
%        'test_',
%        '1',
%        '.tif'
%
%        getFileNameBody('test_1.tif','_') returns
%        []
%        'test',
%        '1',
%        '.tif'
%   
%        getFileNameBody('C:\mydir\test1.tif') returns
%        'C:\mydir'
%        'test',
%        '1',
%        '.tif'
%
% SEE ALSO fileparts
%
% original author unknown. regexp formulation implemented by jonas
% modified 09/05/2010 by DLoerke

% initialize
body1 = [];
no1 = [];
bridge1 = [];

% compare filenames
bvec = (fname1 == fname2);
% rightmost position where strings are different
posr = max(find(bvec==0));

% bridge = positions right of this cutoff
bridge = fname1(posr+1:length(fname1));

% root is remainder
root1 = fname1(1:posr);
root2 = fname2(1:posr);

% in the remaining name: search for numbers 
numberPos = regexp(root1, '\d');
numberPosd = diff(numberPos);

% cut any non-subsequent positions
np_nonone = max(find(numberPosd~=1));
numberPos(1:np_nonone) = [];

posl = min(numberPos);

body1 = fname1(1:posl-1);
body2 = fname2(1:posl-1);
if strcmp(body1,body2)
    body = body1;
end

no1 = fname1(posl:posr);
no2 = fname2(posl:posr);


end % of function   
      
