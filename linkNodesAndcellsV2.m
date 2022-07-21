function [nodeCellMat,nodeNodeMat,cellCellMat] = linkNodesAndcellsV2(mpm_nodes,mpm_cells,currentImage);
%CellsNearNode identifies the cells around all nodes and writes the matric PCCMatrix that houses the position of the node, class of node and cells surrounding the node.
%
%	INPUT:  mpm_nodes: mpm with node positions%
%           mpm_cells: mpm with cell positions
%           currentImage: segmented image
%
%	OUTPUT: nodeCellMat: 2-d matrix with node(row) to cell(columns)
%               connectivity
%           nodeNodeMat: 2-d matrix with node to node connectivity
%           cellCellMat: 2-d matrix with cell to cell connectivity

% Created by: Antonio Nava Jr.
% Version 1.0 - Last Modified 8/31/10
% modified 09/17/2010 DL
% modified 11/15/2016 TV  (line 140 changed if nn2<num_cells_use to  "<=")


%---Main---
[sx,sy] = size(currentImage);

[num_nodes,ny] = size(mpm_nodes);
[num_cells,ny2] = size(mpm_cells);
% Note: in node analysis of previously tracked cells, many entries in the
% mpm-cells matrix may be nan, because the positions refer to cells that
% appeared in earlier image frames
fpos_cells = find(isfinite(mpm_cells(:,1)));
num_cells_use = length(fpos_cells);

% identify cells by numbering them in the segmented image - unless the
% numbers are already pre-determined (e.g. in a previous tracking step)
% this is done directly with BWlabel
intMin = min(currentImage(:));
intMax = max(currentImage(:));
if intMax==1
    bwimage = bwlabel(currentImage,8);
    % NOTE: any objects that touch the edge of the image are set to value
    % -1 to mrak them as background - this will ensure in the subsequent
    % analysis that these cells are not used for direct characterization
    edgepixels1 = currentImage(:,1);
    edgepixels2 = currentImage(1,:);
    edgepixels3 = currentImage(:,sy);
    edgepixels4 = currentImage(sx,:);
    edgeAreas = unique( [edgepixels1(:),edgepixels2(:),edgepixels3(:),edgepixels4(:)] );
    for e=1:length(edgeAreas)
        cvalue = edgeAreas(e);
        if cvalue>0
            bwimage(bwimage==cvalue) = -1;
        end
    end   
else
    bwimage = currentImage;
end

% size of area for ring positions
wsize = 2;
% positions of the pixels in the ring around a square with size wsize
ringposvec = [1:wsize,wsize*[2:wsize-1],[wsize*wsize:-1:(wsize-1)*wsize+1],1+wsize*[(wsize-2):-1:1]];

nodeCellMat = zeros(num_nodes,num_cells+1);

% figure; compareInfo2(currentImage, mpm_cells);

   
for nn = 1:num_nodes
    
    % local area
    cpos = round(mpm_nodes(nn,:));
    xstart = max(1,cpos(1)-wsize); 
    xend = min(sx,cpos(1)+wsize);
    ystart = max(1,cpos(2)-wsize); 
    yend = min(sy,cpos(2)+wsize);   % Changed sx to sy on May 5, 2021, T.V.
            
    carea = bwimage(xstart:xend,ystart:yend);

    % touching cells
    touchingCells = unique( carea(:) );

    % set background to numcells +1
    touchingCells_rev = touchingCells;
    touchingCells_rev(touchingCells==-1) = num_cells+1;
    touchingCells_rev(touchingCells==0) = [];

    % enter into matrix
    nodeCellMat(nn,touchingCells_rev) = 1;

    
end


% find connected nodes and interfacing cells
nodeNodeMat = zeros(num_nodes,num_nodes);
cellCellMat = zeros(num_cells,num_cells+1);

for nn1=1:num_cells_use-1
    
    fprintf(' cell %05d',nn1);
    
    % nodes associated with this cell
    % fnodes1 = find(nodeCellMat(:,nn1)>0);
       
    for nn2 = nn1+1:num_cells_use
        % nodes associated with this cell
        % fnodes2 = find(nodeCellMat(:,nn2)>0);
        
        % intersection
        %intersectionNodes = intersect(fnodes1,fnodes2);
        
        % direct line readout
        % intersectionNodes = find(sum(nodeCellMat(:,[nn1,nn2]),2)==2);
        
        % readout position
        readpos1 = fpos_cells(nn1);
        readpos2 = fpos_cells(nn2);
        
        intersectionVector = ( nodeCellMat(:,readpos1) & nodeCellMat(:,readpos2) );
        intersectionLength = sum(intersectionVector);
        
%         intersectionVector = ( nodeCellMat(:,nn1) & nodeCellMat(:,nn2) );
%         intersectionLength = sum(intersectionVector);
        
        % if the cells have one node in common, they count as having
        % connection status 2; if they have 2 nodes in common, then they
        % count as having connection status 1 (interface), and the nodes
        % count as being directly connected (IF none of the cells is the
        % background)
       
        %if length(intersectionNodes)>0
        if intersectionLength>0
            
            intersectionNodes = find(intersectionVector);
            
            cellCellMat(readpos1,readpos2) = 2;
            cellCellMat(readpos2,readpos1) = 2;
            
            if length(intersectionNodes)>1
                cellCellMat(readpos1,readpos2) = 1;
                cellCellMat(readpos2,readpos1) = 1;
                
                if nn2<=num_cells_use   % < switched to <= on Nov 15/2016 from Roopa's investigation.
                    nodeNodeMat(intersectionNodes(1),intersectionNodes(2)) = 1;
                    nodeNodeMat(intersectionNodes(2),intersectionNodes(1)) = 1;
                end
            end
                
        end
                
    end % of for nn2
    
    fprintf('\b\b\b\b\b\b\b\b\b\b\b');
    
end % of for nn1
            
nodeNodeMat = sparse(nodeNodeMat);
cellCellMat = sparse(cellCellMat);
nodeCellMat = sparse(nodeCellMat);


end% of function


function []= compareInfo2(ImageBWlabel, mpm_cells);


np = size(mpm_cells,1);

for kk=1:np
    px = round(mpm_cells(kk,1)); py = round(mpm_cells(kk,2));
    if isfinite(px)
        compareVec(kk) = ImageBWlabel(py,px);
    end
end

plot(compareVec,'r.-')
    


end % of subfunction