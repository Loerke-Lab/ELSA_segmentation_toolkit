function [ matrixBWlabel_track ] = trackArea2( matrixBWlabel, refPlane )
%trackArea tracks individual cells over space or time using the information
%from the entire cell area (as opposed to the cell centroid position)
%INPUT:     matrixBWlabel   = matrix containing cell data; contiguous cell
%                               areas are numbered 1-n; cell outlines have
%                               value 0, and background areas (not used for
%                               tracking) have value -1
%           refPlane        = reference plane in matrix used as start point
%                               for tracking

[sx,sy,sz] = size(matrixBWlabel);
if nargin<2
    rp = 1;
else
    rp = refPlane;
end

% number of features in each plane
for n=1:sz    
    cplane_n = matrixBWlabel(:,:,n);
    cvalu_n = unique(cplane_n(:));
    clp = length(find(cvalu_n>0));
    uniqueRegionValues(n).vec = cvalu_n;
    np(n) = clp;
    np_max(n) = max(cvalu_n);
end
    
fprintf(' perform individual frame-to-frame links: ');    
% link from reference plane to top
for n=rp:sz-1
    
    fprintf('frame # %04d',n);
    cplane_n = matrixBWlabel(:,:,n);
    cvalu_n = unique(cplane_n(:));
    
    cplane_np1 = matrixBWlabel(:,:,n+1);
    cvalu_np1 = unique(cplane_np1(:));
    lp = length(find(cvalu_np1>0));
       
    % value matrix: relative overlap of different cell combinations
    [overlapMat] = linkFrame2Frame_ver2(cplane_n,cplane_np1);
    
    % evaluate matrix to determine best match
    [lapx,lapy,connectMatrix] = conflictResolution(overlapMat);
    
    PointConnections(n).connectMatrix = connectMatrix;
    PointConnections(n).lapx = lapx;
    PointConnections(n).lapy = lapy;
    GlobalConnectionMat(1:length(lapx),n) = lapx;
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');
        
end
fprintf('\n');

% stitch connectivity information together into trackInfo matrix, where
% rows are continuous trajectories, columns are frames, and the values in
% the amtrix denote the number of the cell in that frame

trackInfo = zeros(max(np),sz);
trackInfo(1:np(1),1) = 1:np(1);

%% compile a list of all continuous trajectories in the structure traj; this
% contains the start frame and a list of the cell numbers in all
% consecutive tracked frames

% counter for trajectories
ct=1;
fprintf('compile list of closed trajectories...');
% loop over all frames
for c=1:sz-1
    % connection values for this frame
    cpos_c = GlobalConnectionMat(:,c);
    %cpos_c1 = GlobalConnectionMat(:,c+1);
    
    % determine the positions that corresponds to actual features (as
    % opposed to connection placeholders), which are the positions up to
    % np(c), and the positions that still have to be tracked (as opposed to
    % those that have been 'accounted for' in a previous trajectory linking
    % step), which are those that have non-zero entries
    
%     fpos = find(cpos_c>0);
%     fpos = fpos( find(fpos<=np_max(c)) );
    regValues = uniqueRegionValues(c).vec;
    fpos = find( (cpos_c>0) & (ismember(cpos_c,regValues)) );
    
    % if any positions in this frame need to be linked - otherwise continue
    % to the next frame
    if length(fpos)>0
        % loop over all trackable positions
        for r=1:length(fpos)
            % index position of current feature
            rpos = fpos(r);
            % record startframe for this new trajectory
            traj(ct).startframe = c;
            % initialize cell number vector for this trajectory
            traj_novec = rpos;
            
            % initialize moving row index
            ind_r = rpos;
            % initialize moving column index
            ind_c = c;
            % link for this index
            cval = GlobalConnectionMat(ind_r,ind_c);
            
            % if this link points to an actual feature (value of cval is up
            % to the number of existing features for the subsequent frame),
            % record the link, update the indices and continue to next
            % frame; if the value of cval is above the number of features,
            % the link is 'empty' and the trajectory ends
            while ( cval <= np_max(ind_c+1) ) 
                % add link to cell number vector
                traj_novec = [traj_novec cval];
                % set the corresponding value in the GlobalConnectionMat to
                % zero, to mark this point as already being part of a
                % trajectory
                GlobalConnectionMat(ind_r,ind_c) = 0;
                % update row index to the value specified in cval
                ind_r = cval;
                % increase column index by one
                ind_c = ind_c+1;
                % break the loop if we have reached the maximum number of
                % linkable frames
                if ind_c > sz-1
                    break
                end
                % otherwise extract the new link for this index
                cval = GlobalConnectionMat(ind_r,ind_c);
                 
            end
            % record the cell number vector
            traj(ct).cellnumbers = traj_novec;
            % step up the index for trajectories
            ct = ct+1;
        end
    end   
end
fprintf(' - total number = %05d',ct-1);
fprintf('\n');

% write trajectories into a mtraix of the trackInfo type
trackInfo = zeros(length(traj),sz);
for c=1:length(traj)
    c_startpoint = traj(c).startframe;
    c_cellnumbers = traj(c).cellnumbers;
    trackInfo(c,c_startpoint:c_startpoint+length(c_cellnumbers)-1) = c_cellnumbers;
end
    
% overwrite numbers in matrixBWlabel with new numbers specified in
% trackInfo
matrixBWlabel_track = matrixBWlabel;
% display current progress of processing
fprintf('renumber cells based on links: ');
    
for r=1:length(traj)
    
    fprintf('traj # %05d',r);
    
    c_startpoint = traj(r).startframe;
    c_cellnumbers = traj(r).cellnumbers;
    c_length = length(c_cellnumbers);
    
    for c=1:c_length
        column_number = c_startpoint+c-1;
        cframe_oimage = matrixBWlabel(:,:,column_number);
        cframe_nimage = matrixBWlabel_track(:,:,column_number);
        fpos_pixels = find( cframe_oimage == c_cellnumbers(c) );
        if ~isempty(fpos_pixels)
            cframe_nimage(fpos_pixels) = r;
            matrixBWlabel_track(:,:,column_number) = cframe_nimage;
        end
    end
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');

end

fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b\b');

fprintf('\n');

end % of function



%% =======================================================================
%                           
%                               subfunctions
%
% =======================================================================



function [overlapMat] = linkFrame2Frame(cplane_n,cplane_np1)

uval_n = unique(cplane_n(:));
uval_np1 = unique(cplane_np1(:));

% remove 0 and -1
fpos_n = find( (uval_n==0) | (uval_n==-1) );
fpos_np1 = find( (uval_np1==0) | (uval_np1==-1) );

if length(fpos_n)>0
    uval_n(fpos_n) = [];
end

if length(fpos_np1)>0
    uval_np1(fpos_n) = [];
end

overlapMat = [];

if min(length(uval_n),length(uval_np1))>0

    for i=1:length(uval_n)

        cval_i = uval_n(i);
        binaryImage_i = (cplane_n==cval_i);

        % size of this cell area
        imageIPixels = sum(binaryImage_i(:));

        for k=1:length(uval_np1)

            cval_k = uval_np1(k);
            binaryImage_k = (cplane_np1==cval_k);

            % overlap in pixels
            overlapImage = (binaryImage_i & binaryImage_k);
            overlapPixels = sum(overlapImage(:));

            % relative overlap
            overlapRelative = overlapPixels/imageIPixels;

            overlapMat(i,k) = overlapRelative;

        end % of for k

    end % of for i
    
end % of if

end % of subfunction
        
        

function [overlapMat] = linkFrame2Frame_ver2(cplane_n,cplane_np1)

uval_n = unique(cplane_n(:));
uval_np1 = unique(cplane_np1(:));
% remove 0 and -1
fpos_n = find( (uval_n==0) | (uval_n==-1) );
if length(fpos_n)>0
    uval_n(fpos_n) = [];
end
fpos_np1 = find( (uval_np1==0) | (uval_np1==-1) );
if length(fpos_np1)>0
    uval_np1(fpos_np1) = [];
end

% set unrelevant values (background and edges at -1 and 0) in images to nan
cplane_n(find(cplane_n<1))=nan;
cplane_np1(find(cplane_np1<1))=nan;
% for tracking, we only need to consider cells that have some overlap
% between frames, and cells with zero overlap between the cell areas can be
% excluded a priori from tracking
% to determine which cells overlap for this frame combination, we first 
% combine the values of the two images into a unique value
ccomb = cplane_n + 1./(1+cplane_np1);
% for unique function, set nans back to zero
ccomb(find(isnan(ccomb)))=0;
% extract unique values (corresponding to combinations of cell overlaps
% that exist between the images)
uvec = unique(ccomb);

% now reconstruct unique combinations of cells from unique numbers
uvec(find(uvec==0))=[];

overlapMat = zeros(max(uval_n),max(uval_np1));

for u=1:length(uvec)
    % first contribution
    n1 = floor(uvec(u));
    % second contribution
    n2 = round( (1/(uvec(u)-n1)) - 1);
    uvec2(u,1:2) = [n1 n2];
    %original area for this cell
    uvec2(u,3) = length(find(cplane_n==n1)); 
    % overlap area for this combination
    uvec2(u,4) = length(find(ccomb==uvec(u)));  
    % overlap percentage
    uvec2(:,5) = uvec2(:,4)./uvec2(:,3);
    % write results into conflict resolution matrix
    overlapMat(n1,n2) = uvec2(u,5);
end
    
overlapMat = sparse(overlapMat);

end % of subfunction
        
        
function [lapx,lapy,connectMatrix] = conflictResolution(overlapMat);

% revert back to full
overlapMat_use = full(overlapMat);
[sx,sy] = size(overlapMat_use);
% set zeros (no overlap of areas) to non-link value
fpos = find(overlapMat_use==0);
% revert overlap values to represent cost
overlapMat_use = 1-overlapMat_use;
overlapMat_use(fpos) = -1;

% extend matrix for lap linking
overlapMat_extended = ones(sx+sy);
overlapMat_extended(1:sx,1:sy) = overlapMat_use;

% lap linking
[lapx, lapy] = lap(overlapMat_extended, -1);

connectMatrix_extended = 0*overlapMat_extended;
for i=1:length(lapx)
    connectMatrix_extended(i,lapx(i)) = 1;
end
connectMatrix = connectMatrix_extended(1:sx,1:sy);
connectMatrix = sparse(connectMatrix);

end % of subfunction

