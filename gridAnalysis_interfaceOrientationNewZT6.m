function [ matRes,  trackingMatrixZT, trackingVectorLI] = gridAnalysis_interfaceOrientationNewZT6( path, anglevector, filelist, tvec );
% gridAnalysis_interfaceOrientation determines the spatial orientation of 
% interfaces between cells 
% INPUT:    path:       directory with the image data, which should 
%                       contain a subfolder called SegmentationData
%           anglevector: bin vector for spatial angles (0=x-axis, 90=
%                       y-axis)
%           tvec:       vector containing desired first and last frame,
%                       e.g. [1 10]; default is from first to last
%                       available frame
%           
% OUTPUT:   matRes:    matrix with average length of interfaces into a
%                       section of spatial angle
%           position matRes(t,f,:) contains the thresholded vector for time
%           point t and section f, where the vector has the average lengths
%           for the bins corresponding to anglevector
%           trackingMatrixZT: tracking matrix contains
%           column 1:     length of interface
%           column 2:     theta (calculated from delta x and delta y)
%           column 3:     x0
%           column 4:     y0
%           column 5:     delta x
%           column 6:     delta y
%           column 7:     identity # of first cell
%           column 8:     identity # of second cell
%
%   Note: This function assumes that cells have been tracked in z and t and
%   that nodes and interfaces have been calculated based on the t-tracked
%   cells               
%        
% written 09/17/2010 DL
% modified 04/18/2011 DL
% modified 01/08/2013 DL
% modified 05/07/2021 TV (I added line 252 to fix an error with "frozen
% nodes".

% record original directory
od = cd;

if nargin>0
    if ~isempty(path)
        cd(path);
    end
end

binvector_angles_centers    = [0:10:360];
binvector_angles_binedges   = [0:10:360] + 355;
if nargin>1
    if ~isempty(anglevector)
        binvector_angles_centers = anglevector;
        width = anglevector(2)-anglevector(1);
        binvector_angles_binedges   = binvector_angles_centers + (360 - width/2 );
    end
end

% change directory to location of segmentation data and record that, too
cd('SegmentationData');
od2 = cd;

% loop over all existing time points
t_start = 1;
t_end = length(filelist);
if nargin>3
    if ~isempty(tvec)
        t_start = tvec(1);
        t_end = tvec(2);
    end
end

t = t_start;
cframefoldername = sprintf('frame%04d',t);

while exist(cframefoldername)==7
    
    % forced exit of the while-loop if the current t exceeds t_end
    if t>t_end
        break
    end
    
    % otherwise continue analysis
    cd(cframefoldername);
    
    fprintf('frame %04d',t);
        
    loadstruct = open('dstruct_nodes.mat');
    % structure with positions of nodes
    nodePosStruct = loadstruct.dstruct_nodes;
    
    loadstruct = open('dstruct_nodeNodeMat.mat');
    % structure with node-node connectivity
    nodeNodeStruct = loadstruct.dstruct_nodeNodeMat;
    
    loadstruct = open('dstruct_nodeCellMat.mat');
    % structure with node-cell connectivity
    nodeCellStruct = loadstruct.dstruct_nodeCellMat;

    % loop over all existing sections points
    for k=1:length(nodeNodeStruct)
        
        fprintf(' section %04d',k);
        
        % node-node matrix for current section
        nnmat = full(nodeNodeStruct(k).matrix);
        [sx,sy] = size(nnmat);
        
        % node-cell matrix for current section
        ncmat = full(nodeCellStruct(k).matrix);
        % reminder: rows = node numbers; columns = cell numbers        
        
        nodepos = nodePosStruct(k).positions;
        
        %initialize local results
        totalpos = length(find(nnmat>0));
        cresults_k = nan*zeros(totalpos,2);
        ct = 0;
                           
        % for every node, find its connections and the line connectors
        for n=1:sx
            
            % reference node position
            npos_x0 =  nodepos(n,1);
            npos_y0 =  nodepos(n,2);
            
            % what cells (i.e. with what numbers) meet at this ref node?
            ncellvec_n0 = ncmat(n,:);
            
            % nodes connected to this node = positions in matrix with
            % connectivity value=1
            cline = nnmat(n,:);
            cpos = find(cline>0);
            lpos = length(cpos);
                        
            % vector of node positions connected to the reference node
            npos_x = nodepos(cpos,1);
            npos_y = nodepos(cpos,2);
            
            % connection vector between reference and target nodes
            npos_xdiff = npos_x - npos_x0;
            npos_ydiff = npos_y - npos_y0;
            
            % connection distance
            nrad = sqrt( npos_xdiff.^2 + npos_ydiff.^2 );
            
            % what cells (i.e. with what numbers) meet at these connected
            % nodes?
            ncellvec_n = ncmat(cpos,:);
            
            % now determine the numbers of the two cells that meet at each
            % interface
            ncellno = zeros(1,2);
            for r=1:lpos
                rpos = find( sum([ncellvec_n0; ncellvec_n(r,:)],1) == 2 );
                if length(rpos)==2
                    ncellno(r,:) = rpos;
                else
                    ncellno(r,:) = [0 0];
                end
            end
                    
            % enter results into results matrix: first column radius
            cresults_k(ct+1:ct+lpos,1) = nrad;
            
            % enter results into results matrix: second/third column x,y  
            cresults_k(ct+1:ct+lpos,2) = npos_xdiff;
            cresults_k(ct+1:ct+lpos,3) = npos_ydiff;
            
            % for graphical representation further down: retain x0,y0
            cresults_k(ct+1:ct+lpos,4) = npos_x0;
            cresults_k(ct+1:ct+lpos,5) = npos_y0;
            
            % for tracking of individual interfaces: retain cell identities
            % leave column 6 free for theta
            cresults_k(ct+1:ct+lpos,7) = ncellno(:,1);
            cresults_k(ct+1:ct+lpos,8) = ncellno(:,2);
                       
            ct = ct+lpos;
        end % of for n-loop
               
        % extract all x and y vector positions
        xall = cresults_k(:,2);
        yall = cresults_k(:,3);
        
        % to calculate the corresponding spatial angle in a full 0-360
        % degree range, first identify the relevant quadrant:
        fpos_q1 = find( (xall>0) & (yall>=0) );
        fpos_q2 = find( (xall<=0) & (yall>0) );
        fpos_q3 = find( (xall<0) & (yall<=0) );
        fpos_q4 = find( (xall>=0) & (yall<0) );
        
        %...then calculate the angle in that quadrant using the tan
        %function
        thetaDeg_q1 = atan(yall(fpos_q1)./xall(fpos_q1)) * (360/(2*pi));
        thetaDeg_q2 = 90 + atan(abs(xall(fpos_q2))./yall(fpos_q2)) * (360/(2*pi));
        thetaDeg_q3 = 180 + atan(abs(yall(fpos_q3))./abs(xall(fpos_q3))) * (360/(2*pi));
        thetaDeg_q4 = 270 + atan(abs(xall(fpos_q4))./abs(yall(fpos_q4))) * (360/(2*pi));
%         thetaDeg_q3 = 0 + atan(abs(yall(fpos_q3))./abs(xall(fpos_q3))) * (360/(2*pi));
%         thetaDeg_q4 = 90 + atan(abs(xall(fpos_q4))./abs(yall(fpos_q4))) * (360/(2*pi));
        
        % enter the angle (in degrees) into the results matrix (fourth
        % column)
        cresults_k(fpos_q1,6) = thetaDeg_q1;
        cresults_k(fpos_q2,6) = thetaDeg_q2;
        cresults_k(fpos_q3,6) = thetaDeg_q3;
        cresults_k(fpos_q4,6) = thetaDeg_q4;
        
        % threshold the results (distance vs. angle increment) using a bin 
        % vector for the angles defined at the start of the function
        
        % double data for wraparound at +360 degrees
        ref1 = cresults_k(:,6);
        ref2 = cresults_k(:,6)+360;
        ref_compl = [ ref1; ref2 ];
        data1 = cresults_k(:,1);
        data2 = cresults_k(:,1);
        data_compl = [ data1; data2 ];
        threshVector_k = thresholdVector(ref_compl,data_compl,binvector_angles_binedges); 
        
        % enter thresholded binvector into final results matrix
        matRes(t,k,:) = threshVector_k(2,:);
        
        fprintf('\b\b\b\b\b\b\b\b\b\b\b\b\b');
        
        % maximal cell number in columns 7 and 8
        maxncell1 = max(cresults_k(:,7));
        maxncell2 = max(cresults_k(:,8));
        % for sorting, define linear index from columns 7 and 8 (add 1
        % because of possible zeros)
        linearIndex = sub2ind([maxncell2+1 maxncell1+1], cresults_k(:,8)+1, cresults_k(:,7)+1);
        cresults_k(:,9) = linearIndex;
        
        % at this point, the matrix cresults_k contains the following data
        % for each interface (for this time point and section):
        % column 1:     length of interface
        % column 2:     delta x
        % column 3:     delta y
        % column 4:     x0
        % column 5:     y0
        % column 6:     theta (calculated from delta x and delta y)
        % column 7:     identity # of first cell
        % column 8:     identity # of second cell
        % column 9:     identity # of interface 
            
        % sorted matrix containing length, angle, id numbers of
        % interfaces, and x0, y0; is sorted according to linear id#
        cresults_ksort = sortrows( cresults_k(:,[1,6,7,8,9,4,5,2,3]),5 );
        
        % since most interfaces are entered double, choose unique
        % identifiers and write only the entry with the smaller angle theta
        % into the tracking matrix of interfaces
        indices_use = unique(cresults_ksort(:,5));
        indices_use = indices_use(indices_use>1);
        
        cresults_ksort_crop = []; % This line fixes the frozen nodes problem. The variable needs to be reinitialized to empty or data from previous frames may remain.
        for ci=1:length(indices_use)
            ci_fpos = find(cresults_ksort(:,5)==indices_use(ci));
            cresults_ksort_crop(ci,1) = cresults_ksort(ci_fpos(1),1);
            cresults_ksort_crop(ci,2) = min(cresults_ksort(ci_fpos,2));
            cresults_ksort_crop(ci,3) = cresults_ksort(ci_fpos(1),3);
            cresults_ksort_crop(ci,4) = cresults_ksort(ci_fpos(1),4);
            cresults_ksort_crop(ci,5) = cresults_ksort(ci_fpos(1),5);
            cresults_ksort_crop(ci,6) = cresults_ksort(ci_fpos(1),6);
            cresults_ksort_crop(ci,7) = cresults_ksort(ci_fpos(1),7);
            cresults_ksort_crop(ci,8) = cresults_ksort(ci_fpos(1),8);
            cresults_ksort_crop(ci,9) = cresults_ksort(ci_fpos(1),9);
        end
        
        % for subsequent construction of a tracking matrix, retain the
        % length, angle, and cell id information
        InterfaceData{t,k}.mat = cresults_ksort_crop(:,[1:4, 6:9]);
        
        
        % at this point, the matrix in InterfaceData contains the following data
        % for each interface (for this time point and section):
        % column 1:     length of interface
        % column 2:     theta (calculated from delta x and delta y)
        % column 3:     identity # of first cell
        % column 4:     identity # of second cell
        % column 5:     x0
        % column 6:     y0
        % column 7:     delta x
        % column 8:     delta y
          
        
        %% display results for this time point and frame
        % read original image
        imageOriginalName = char(filelist{t,k});
        imageOriginal = imread(imageOriginalName);
        
        hold off;
        imshow(imageOriginal(:,:,1),[]); hold on;
        set(gca, 'position',[0.01 0.01  0.98  0.98]);
        
                
        % interfaces at zero degrees angles
        fpos_th0 = find( ismember( cresults_k(:,6),[0,360] ) );
        
        x1_th0 = cresults_k(fpos_th0,4);
        x2_th0 = cresults_k(fpos_th0,4) + cresults_k(fpos_th0,2);
        y1_th0 = cresults_k(fpos_th0,5);
        y2_th0 = cresults_k(fpos_th0,5) + cresults_k(fpos_th0,3);
        
        % interfaces at 90 degrees angles
        fpos_th90 = find( ismember(cresults_k(:,6),[90,270] ) );
        
        x1_th90 = cresults_k(fpos_th90,4);
        x2_th90 = cresults_k(fpos_th90,4) + cresults_k(fpos_th90,2);
        y1_th90 = cresults_k(fpos_th90,5);
        y2_th90 = cresults_k(fpos_th90,5) + cresults_k(fpos_th90,3);
        
                
        % interfaces in first angle quadrant (0-90)
        x1_q1 = cresults_k(fpos_q1,4);
        x2_q1 = cresults_k(fpos_q1,4) + cresults_k(fpos_q1,2);
        y1_q1 = cresults_k(fpos_q1,5);
        y2_q1 = cresults_k(fpos_q1,5) + cresults_k(fpos_q1,3);
        
        % interfaces in second angle quadrant (90-180)
        x1_q2 = cresults_k(fpos_q2,4);
        x2_q2 = cresults_k(fpos_q2,4) + cresults_k(fpos_q2,2);
        y1_q2 = cresults_k(fpos_q2,5);
        y2_q2 = cresults_k(fpos_q2,5) + cresults_k(fpos_q2,3);
        
%         % third quadrant points
%         x1_q3 = cresults_k(fpos_q3,4);
%         x2_q3 = cresults_k(fpos_q3,4) + cresults_k(fpos_q3,2);
%         y1_q3 = cresults_k(fpos_q3,5);
%         y2_q3 = cresults_k(fpos_q3,5) + cresults_k(fpos_q3,3);
%         
%         % fourth quadrant points
%         x1_q4 = cresults_k(fpos_q4,4);
%         x2_q4 = cresults_k(fpos_q4,4) + cresults_k(fpos_q4,2);
%         y1_q4 = cresults_k(fpos_q4,5);
%         y2_q4 = cresults_k(fpos_q4,5) + cresults_k(fpos_q4,3);
        
        plot([y1_q1 y2_q1]',[x1_q1 x2_q1]','c-');
        plot([y1_q2 y2_q2]',[x1_q2 x2_q2]','g-');
        
%         plot([y1_q3 y2_q3]',[x1_q3 x2_q3]','y.-');
%         plot([y1_q4 y2_q4]',[x1_q4 x2_q4]','r.-');
        
        plot([y1_th0 y2_th0]',[x1_th0 x2_th0]','r-');
        plot([y1_th90 y2_th90]',[x1_th90 x2_th90]','y-');
        
        text(10,20,'orientation 0 degrees','Color','red','Fontsize',16);
        text(10,50,'orientation 90 degrees','Color','yellow','Fontsize',16);
        
        
        hold off
        pause(0.1);
           
    end % of for k-loop (sections)
    
    % return to parent directory
    cd(od2);
    % step index forward
    t=t+1;
    % update folder name for while-loop
    cframefoldername = sprintf('frame%04d',t);
    
    if nargin>3
        if ~isempty(tvec)
            tmax=tvec(2);
            if t>tmax
                break
            end % of if
        end % of if
    end
   
    
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    
end % of while loop over frames

fprintf('\n');

cd(od);

% determine highest available cell number for all available t's and z's
[six,siy]=size(InterfaceData);
for t=1:six
    for z=1:siy
        cmat = InterfaceData{t,z}.mat;
        cmax = max( max(cmat(:,3)), max(cmat(:,4)) );
        cmaxmat(t,z) = cmax;
    end
end
totalmax = max(cmaxmat(:));

% continue here: construct 3-D tracking matrix (layers = z-sections)
disp('constructing 3-D tracking matrix:');

for z=1:siy
          
    fprintf('layer %04d',z);
    
    for t=1:six
        
        fprintf(' frame %04d',t);
        
        % at this stage, the matrix in InterfaceData contains the following
        % data for each interface (for this time point and section):
        % column 1:     length of interface
        % column 2:     theta (calculated from delta x and delta y)
        % column 3:     identity # of first cell
        % column 4:     identity # of second cell
        % column 5:     x0
        % column 6:     y0
        % column 7:     delta x
        % column 8:     delta y
        
        cmat = InterfaceData{t,z}.mat;
        % first column: length of interface in pixels
        cc1 = cmat(:,1);
        % second column: angle of interface in degree
        cc2 = cmat(:,2);
        % third column: identifier of first cell
        cc3 = cmat(:,3);
        % fourth column: identifier of second cell
        cc4 = cmat(:,4);
        % fifth column: x0 position
        cc5 = cmat(:,5);
        % sixth column: y0 position
        cc6 = cmat(:,6);
        % seventh column: delta x 
        cc7 = cmat(:,7);
        % eighth column: delta y
        cc8 = cmat(:,8);
        % as a result, the linearIndex (a combination of cell numbers) is a
        % unique identifier of this particular cell combination
        linearIndex = sub2ind([totalmax totalmax], cc4, cc3);
        %linIndex{t,z}.vec = linearIndex;
        
        % for the tracking matrix, we are only interested in interfaces
        % that have trajectories of more than one frame (for whcih we can
        % determine a delta-d or delta-theta); thus, we will save space by
        % ignoring interfaces that appear in only one frame
        % this is achieved by determining which of the interfaces in this
        % frames also appeared in the previous or subsequent frame
        
        % for the first frame, the initial ids in frame n-1 are empty, and
        % thus the intersection of interfaces present in frame n-1 and n is
        % empty,too
        if t==1
            goodIDs_A = [];
        % for all other frames, since we're counting up t, the intersection
        % of n-1 and n is identical to the intersection of n and n+1 from
        % the previous value of t
        else
            goodIDs_A = goodIDs_B;
        end % of if
            
        % for frames from t=1 up to t=(six-1), determine ids in frame n+1
        if t<=six-1
            cmat_np1 = InterfaceData{t+1,z}.mat;
            % third column: identifier of first cell
            cc3_np1 = cmat_np1(:,3);
            % fourth column: identifier of second cell
            cc4_np1 = cmat_np1(:,4);
            % again linearIndex (a combination of cell numbers) is a
            % unique identifier of this particular cell combination
            linearIndex_np1 = sub2ind([totalmax totalmax], cc4_np1, cc3_np1);
            %linIndex{t,z}.vec = linearIndex;
            
            % intersection of n-1 and n
            goodIDs_B = intersect(linearIndex, linearIndex_np1);
        end % of if
        
        % for final frame, ids in frame n+1 (and thus intersection) are empty 
        if t==six
            goodIDs_B = [];
        end % of if
        
        % thus, the good IDs for this frame are the union of (n intersects 
        % with n-1) and of (n intersects with n+1)
        goodIDs_total = union(goodIDs_A, goodIDs_B);
        
        goodPos = find(ismember(linearIndex,goodIDs_total));
        
        % now enter appropriate values into the tracking matrix
        
        % initialize results for t=1 and z=1
        if t==1 & z==1
            trackingMatrixZT = [cc1(goodPos), cc2(goodPos), cc5(goodPos), cc6(goodPos)...
                cc7(goodPos), cc8(goodPos), cc3(goodPos), cc4(goodPos) ];
            trackingVectorLI = [cc3(goodPos) cc4(goodPos) linearIndex(goodPos)];
        
        % else determine appropriate position (in matrix if linear index
        % value already exists, add to bottom otherwise)
        else
            % loop over all entries
            for li=1:length(goodPos)
                % linear index (ID) value of every 'good' trajectory
                cli = linearIndex(goodPos(li));
                % check if any entries (in previous frames or in previous
                % layers) already exist for this ID
                epos = find(trackingVectorLI(:,3)==cli);
                % if no previous entries exist, add the current ID to the
                % bottom of the trackingVectorLI list, and set that
                % position as epos
                if isempty(epos)
                    clen_t = length(trackingVectorLI(:,3));
                    trackingVectorLI(clen_t+1,[1:3]) = [ cc3(li) cc4(li) cli];
                    epos = clen_t+1;
                end
                
                % add interface length, angle orientation, and x0, y0 to 
                % the row epos (and column and layer as appropriate)
                trackingMatrixZT(epos,8*t-7,z) =    cc1(goodPos(li));
                trackingMatrixZT(epos,8*t-6,z) =    cc2(goodPos(li));
                trackingMatrixZT(epos,8*t-5,z) =    cc5(goodPos(li));
                trackingMatrixZT(epos,8*t-4,z) =    cc6(goodPos(li));
                trackingMatrixZT(epos,8*t-3,z) =    cc7(goodPos(li));
                trackingMatrixZT(epos,8*t-2,z) =    cc8(goodPos(li));
                trackingMatrixZT(epos,8*t-1,z) =    cc3(goodPos(li));
                trackingMatrixZT(epos,8*t-0,z) =    cc4(goodPos(li));
                
            end % of for li
        end % of if-else    
        
        fprintf('\b\b\b\b\b\b\b\b\b\b\b');
                        
    end % of for t
    
    fprintf('\b\b\b\b\b\b\b\b\b\b');
    
end % of for z

fprintf('\n');

end % of function



