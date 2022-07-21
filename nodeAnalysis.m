function [dstruct_nodes,dstruct_nodeNodeMat,dstruct_cellCellMat,dstruct_nodeCellMat] = nodeAnalysis(ImageMatrixSegment,cellCentroidStruct,ImageMatrixBWlabel);
% Input:    ImageMatrixSegment: 3-dim matrix with segmented images in
%               subsequent image layers
%           cellCentroidStruct: structure with cell centroid positions 
%           ImageMatrixBWlabel: 3-dim matrix with segmented images and
%               numbered regions in subsequent image layers 
%
% Output:   dstruct_nodes: structure with node positions in .position field
%           dstruct_nodeNodeMat: structure with 2-dim matrices that
%               indicate connectivity of nodes with each other
%           dstruct_cellCellMat: structure with 2-dim matrices that
%               indicate connectivity of cells with each other
%           dstruct_nodeCellMat: structure with 2-dim matrices that
%               indicate connectivity of nodes (rows) with cell (columns)
% 
% author: Julian Hagemeister
% modified 09/17/2010 DL

% ImageMatrixSegment can also be ImageMatrixBWlabel
[sx,sy,sz] = size(ImageMatrixSegment);

% loop over number of specified frames
for z=1:sz
    
    
    %fprintf(' section %03d',z);
    
    % extract current image
    currentImage = ImageMatrixSegment(:,:,z);
    if nargin>2
        currentImageBWlabel = ImageMatrixBWlabel(:,:,z);
    end
        
    %if image is not blank, continue with analysis
    fvec=find(currentImage(:)==0);
    
    if length(fvec)>0
        
        % find nodes in image
        [mpm_nodes] = detectNodes_ver2(currentImage);
        % add current node positions to nodeStruct
        dstruct_nodes(z).positions = mpm_nodes;
        
        mpm_cells = cellCentroidStruct(z).positions;
        
        % based on cell positions and cluster positions, determine which
        % cells meet at which nodes (nodeCellMat), which nodes are directly
        % connected to each other (nodeNodeMat), and which cells share an
        % interface (cellCellMat)
        if nargin>2
            [nodeCellMat,nodeNodeMat,cellCellMat] = linkNodesAndcellsV2(mpm_nodes,mpm_cells,currentImageBWlabel);
        else
            [nodeCellMat,nodeNodeMat,cellCellMat] = linkNodesAndcellsV2(mpm_nodes,mpm_cells,currentImage);
        end
        dstruct_nodeNodeMat(z).matrix = nodeNodeMat;
        dstruct_cellCellMat(z).matrix = cellCellMat;
        dstruct_nodeCellMat(z).matrix = nodeCellMat;
       
        %fprintf('\b\b\b\b\b\b\b\b\b\b\b\b');
    else
        %fprintf('\n');
        %display(['blank image at position ', num2str(z)]);
    end %if loop
      
end % of loop

% save results into the current folder (The following lines were commented
% because the main function handles the saving, fixed 1/7/13).
% save('dstruct_nodeNodeMat', 'dstruct_nodeNodeMat');
% save('dstruct_cellCellMat', 'dstruct_cellCellMat');
% save('dstruct_nodeCellMat', 'dstruct_nodeCellMat');
% save('dstruct_nodes', 'dstruct_nodes');

end % of function