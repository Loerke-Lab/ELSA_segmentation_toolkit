function [mpm_nodes] = detectNodes_ver2(imageSegm, distanceThreshold);
% detectNodes detects and classifies the nodes in a line-segmented image
%
% INPUT:    imageSegm = segmented line image
%           distanceThreshold (optional) = distance threshold for which nodes
%               are merged, default value=2 (pixels)
%
% OUTPUT:   mpm_nodes = mpm file containing detected nodes; first column =
%               x-positions, second column y-positions, third column =
%               class of the node (how many lines are intersecting)
%
% last modified: August 20, 2010 DL


%% Step 0 : (re)set segmImage such that the segmented lines have value 1, 
% and the background is value 0
[sx,sy,sz] = size(imageSegm);

f0 = find(imageSegm(:)==0);
f1 = find(imageSegm(:)==1);

imageSegmSet = imageSegm;
if length(f0)<length(f1)
    imageSegmSet(f0) = 1;
    imageSegmSet(f1) = 0;
end % of if

% find locations where segmented lines are
[fsetx,fsety] = find(imageSegmSet==1);

% set distance threshold to default 2 unless specified
if nargin>1
    dthresh = distanceThreshold;
else
    dthresh = 2;
end


% ring positions size 5 (+-2) neighborhood
rpos2 = [1:5,10,15,20,25:-1:21,16,11,6];
% ring positions size 7 (+-3) neighborhood
rpos3 = [1:7,14,21,28,35,42,49:-1:43,36,29,22,15,8];

% posiitons of 4-type pixel neighborhood
rpos1star = [2,6,8,4];


    
%% Step 1: cycle through the line points to determine whether they are
% nodes of *any* type - through the structure of the watershed segmented
% image, straight or diagonal line points always have only 2 neighbors in
% the immediate 4-pixel neighborhood, while nodes have 3 or more

imageNodeClass_test2 = 0*imageSegmSet;
imageNodeClass_test3 = 0*imageSegmSet;

for k=1:length(fsetx)
    currx = fsetx(k);
    curry = fsety(k);
    
    if (currx>1) & (curry>1) & (currx<sx) & (curry<sy)
        % type 1 neighborhood
        imageNei1 = imageSegmSet(currx-1:currx+1, curry-1:curry+1);
        % extract values in ring positions
        vecRing1star = imageNei1(rpos1star);

        fp_plus1 = find(vecRing1star==1);
        nsteps1 = length(fp_plus1);

        imageNodeClass_test1(currx,curry) = nsteps1;
    end % of if
    
end % of for k-loop

% extract the connected nodes in an 8-pixel neighborhood 
imageNode1BW = (imageNodeClass_test1>2);
imageBWlabel = bwlabel(imageNode1BW,8);

% ntotal number of regions
regmax = max(imageBWlabel(:));

% initialize mpm file for node results
mpm_nodes = zeros(regmax,3);

%% Step 2:
% loop over number of connected regions, and in those regions that consist
% of more than a single pixel, merge the connected pixels (nodes) to a
% single node of appropriate class
for k=1:regmax
    
    % positions where the BWlabel image has the value corresponding to this region
    [fregx,fregy] = find(imageBWlabel==k);
    
    % if the number of pixels is this region equals 1, enter the node into the list
    if length(fregx)==1
        mpm_nodes(k,1:3) = [fregx, fregy, 3];
    
    % else if there's multiple pixels, check each indiviual pixel belonging
    % to the region separately to determine how many distinct lines are
    % leaving its 5 (plus/minus 2) neighborhood
    else

        % initialize fregLines matrix       
        fregLines = [];

        % loop over pixels in this region
        for n=1:length(fregx)
            
            % type 2 neighborhood
            imageNei2 = imageSegmSet(fregx(n)-2:fregx(n)+2, fregy(n)-2:fregy(n)+2);
            % extract values in ring positions
            vecRing2 = imageNei2(rpos2);

            % check how many continuous pieces are in the vector 
            % (how many steps up or down)
            vecRingDiff2 = diff(vecRing2);
            fp_plus2 = find(vecRingDiff2==1);
            fp_minus2 = find(vecRingDiff2==-1);

            % set number of steps as the maximum of up- or downward steps
            nsteps2 = max(length(fp_plus2),length(fp_minus2));

            % set result for current center pixel
            fregLines(n) = nsteps2;

        end
        
        % x-position for this new cluster: center of gravity of x-positions
        % of the contributing pixels
        cpx = nanmean(fregx);

        % y-position for this new cluster: see above
        cpy = nanmean(fregy);

        % new class is maximum of contributing classes
        cclass = max(fregLines);

        % enter new node into results list
        mpm_nodes(k,1:3) = [cpx, cpy, cclass];

    end
    
end
    
% display preliminary results
% figure; imshow(imageSegm,[]);
% hold on;
% 
% fp3 = find(mpm_nodes(:,3)==3);
% fp4 = find(mpm_nodes(:,3)==4);
% fp5 = find(mpm_nodes(:,3)>=5);
% 
% plot(mpm_nodes(fp3,2),mpm_nodes(fp3,1),'r*');
% plot(mpm_nodes(fp4,2),mpm_nodes(fp4,1),'y*');
% plot(mpm_nodes(fp5,2),mpm_nodes(fp5,1),'g*');


%% Step 3:
% loop over all current nodes, and if there are neighboring nodes within a
% specified distance threshold, merge these nodes into a single one of a
% higher class (e.g. merge two nodes of type 3 into one node of type 4)

% distances of all nodes from each other
distMat = DistanceMatrix(mpm_nodes(:,1:2),mpm_nodes(:,1:2));


% loop over all currently existing nodes (eahc corresponding to one line in
% the distance matrix)
for k=1:size(distMat,1)

    % distances from current node
    cdvec = distMat(k,:);

    % positions where distances are below threshold distance
    fdist = find(cdvec<=dthresh);

    % if any nieghbors are within the threshold
    if length(fdist)>1

        % check other positions to determine whether those have additional
        % neighbors within range
        fdist_other = fdist;
        fdist_other(find(fdist==k)) = [];

        % number of neighbors (not including the point itself)
        nn_this = length(fdist_other);
        
        %==================================================================
%         % display current node merger, comment if not necessary
%         ax_x = mpm_nodes(k,1);
%         ax_y = mpm_nodes(k,2);
%         axis([ax_y-20, ax_y+20, ax_x-20, ax_x+20]);
%         pause(0.1);
        %==================================================================

        % check how many qualifying neighbors the points within range have   
        for j=1:length(fdist_other)
            cdvec2 = distMat(fdist_other(j),:);
            fdist2 = find(cdvec2<=dthresh);
            nn_other(j) = length(fdist2);
        end

        % if any of the qualifying neighbors have more neighbors
        % themselves, then skip this entry, and wait for the merged node to
        % be calculated for this other neighbor, using its full number of
        % neighbors
        if (max(nn_other)-1)>nn_this
            continue
        % if the number of qualifying neighbors of the neighbors is less or
        % equal, merge the nodes here
        else
            
            % x-position for this new cluster: center of gravity of 
            % x-positions of the contributing pixels
            cpx = nanmean(mpm_nodes(fdist,1));

            % same for y-positions
            cpy = nanmean(mpm_nodes(fdist,2));

            % merged class: add one for each additional node class above 2
            cclass = 2 + sum(mpm_nodes(fdist,3)-2);
            
            % enter merged result into results file
            mpm_nodes(k,1:3) = [cpx, cpy, cclass];

            % ...and delete the neighbors with which this point was merged
            mpm_nodes(fdist_other,1:3) = nan;
            distMat(fdist_other,:) = nan;

        end

    end % of if

end % of for
        

%% display final results
% figure; imshow(imageSegm,[]);
% hold on;
% 
% fp3 = find(mpm_nodes(:,3)==3);
% fp4 = find(mpm_nodes(:,3)==4);
% fp5 = find(mpm_nodes(:,3)>=5);
% 
% 
% plot(mpm_nodes(fp3,2),mpm_nodes(fp3,1),'r*');
% plot(mpm_nodes(fp4,2),mpm_nodes(fp4,1),'y*');
% plot(mpm_nodes(fp5,2),mpm_nodes(fp5,1),'g*');


%% define output result
fp_isfinite = find(isfinite(mpm_nodes(:,1)));
mpm_nodes = mpm_nodes(fp_isfinite,:);

 
end % of function



%%======================================================================

function [m2]=DistanceMatrix(c1,c2);
%this subfunction makes a neighbour-distance matrix
% 
% INPUT:    c1 - matrix (n1 x 2 points)
%           c2 - matrix (n1 x 2 points)
%
% OUTPUT:   m2 - matrix (n1 x n1) containing the distances of each point in 
%           c1 from each point in c2

% find the size of c1 and c2
[np1,sd1]=size(c1);
[np2,sd2]=size(c2);

% initialize results matrix
m2=zeros(np1,np2);

% loop over the rows of c1 and c2 and 
for k = 1:np1
    for n = 1:np2

        % calculate the distance
        d = sqrt((c1(k,1)-c2(n,1))^2+(c1(k,2)-c2(n,2))^2);

        % save the results in m2 matrix
        m2(k,n)=d;

    end
end


end


