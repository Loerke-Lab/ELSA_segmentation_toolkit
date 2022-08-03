function [body,no1,no2,bridge]=getFilenameBody_two(fname1,fname2)
% GETFILENAMEBODY_TWO splits two strings into their body, number and 
% bridge parts
%
% SYNOPSIS [body,no1,no2,bridge]=getFilenameBody_two(fname1,fname2)
% 
% INPUT fname1 : the following filename structure must be 
%                preserved:
%                - alpha body
%                - number after body and before bridge
%                - bridge separated after differences in strings
%       fname2 : as above
%
% OUTPUT body   : string with body, [] if non-existent
%        no1    : string with the number for fname1, [] if non-existent
%        no2    : string with the number for fname2, [] if non-existent
%        bridge : string with bridge, [] if non-existent
%
% original author unknown. regexp formulation implemented by jonas
% modified 09/05/2010 by DLoerke

% initialize arrays
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

% find the minimum (smallest) number in numberPos
posl = min(numberPos);

% specify the body in each fname before the numbers start
body1 = fname1(1:posl-1);
body2 = fname2(1:posl-1);

% check that body1 and body2 are the same, then set body
if strcmp(body1,body2)
    body = body1;
end

% specify the number for each fname that comes after the body but before the bridge
no1 = fname1(posl:posr);
no2 = fname2(posl:posr);


end % of function   
      
