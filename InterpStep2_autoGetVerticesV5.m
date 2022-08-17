function [X,Y,Z,Vinterp,fitobject,Verts] = InterpStep2_autoGetVerticesV5(V) % ,filtSz
%InterpStep2_autoGetVerticesV4 is an automated method for finding tissue
% surface positions (vertices) from the light sheet microscopy data. The
% only input is the 3D image volume for the cell membrane channel (if one
% exists). This function loops over 2D cross-sections (Y-Z sections)
% and detects the most apical signal to use as surface vertices

% first get the size of the image volume
[Nrows,Ncols,Nlayers] = size(V);

% initialize BW: a mask of where the epithelium is
BW = false(Nrows,Ncols,Nlayers);


% we fit the tissue surface with a 2D polygon of order 2 (parabola) in both
% dimensions, so the fittype is 'Poly22'.
fittype = 'Poly22';


% filter the image for better edge detection 
Vf = imgaussfilt3(V,[1 1 5]);
% find edges in 3D
BWedges = edge3(Vf,'approxcanny',0.15);

% remove all the non-vertical (non-axial) lines
BWedges = imopen(BWedges,cat(3,1,1,1));



% findSkewMask finds the region of the deskewed image that is all zeros, i.e.
% outside of the field of view and also takes a bit of the boundary of the
% image where there are edge effects due to deconvolution.
[mask_blackRegion] = findSkewMask(V,20);

% set the black region mask to zero to ensure none of our edges are
% boundary edges.
BWedges(mask_blackRegion) = false;


% for ii=1:10:Ncols
%     
%     Vslice = squeeze(Vf(:,ii,:));
%     BWslice = squeeze(BWedges(:,ii,:));
%     
%     BWslice = bwconvhull(BWslice);
%     BW(:,ii,:) = BWslice;
%     
%     imoverlay(imadjust(Vslice)',bwperim(BWslice)');
%     title(num2str(ii))
%     pause;
% end

for ii=1:10:Nrows
    
    Vslice = squeeze(Vf(ii,:,:));
    BWslice = squeeze(BWedges(ii,:,:));
    
    BWslice = bwconvhull(BWslice);
    BW(ii,:,:) = BWslice;
    
    imoverlay(imadjust(Vslice)',bwperim(BWslice)');
    title(num2str(ii))
    pause;
end



 

% we will use the convex hull to remove any concavities in the XY
% cross-sections
for ii=1:Nlayers
    BW(:,:,ii) = bwconvhull(BW(:,:,ii));
end

% findSurfaceVerts find the most Apical pixels that are present in the BW
% volume.
[Verts,Surfmask] = findSurfaceVerts(BW);


% Perform a fit of the surface, Z
[X,Y] = ndgrid(1:Nrows,1:Ncols);
fitobject = fit([Verts(:,1),Verts(:,2)],Verts(:,3),fittype);
Z = feval(fitobject,X,Y);


for ii= 1:Nlayers
    Vq = interpn(im2double(V),X,Y,Z+ii-1,'nearest');
    Vinterp(:,:,ii) = Vq;
end



for ii=1:50:Nrows
    
    BWslice = squeeze(Surfmask(ii,:,:))';
    Islice = squeeze(V(ii,:,:))';
    
    imoverlay(imadjust(Islice),logical(BWslice));
    
    %imshow(squeeze(V(ii,:,:))',[]);
    %Islice = I(ii,:,:);
    Yslice = Y(ii,:);
    Zslice = Z(ii,:);
    %Z2slice = Z2(ii,:);
    hold on
    plot(Yslice,Zslice,'LineWidth',2,'Color',[0 1 0]);
    %plot(Yslice,Z2slice,'LineWidth',2,'Color',[0 0 1]);hold off
    title(num2str(ii));
    pause;
    clf
    
end


% for ii=1:10:Ncols
%     
% 
%     I = squeeze(V(:,ii,:))';
%     BWslice = squeeze(BW(:,ii,:))';
% 
% %     [BW2,threshOut,Gv,Gh] = edge3(I,'Sobel');
%     I = imadjust(I);
% 
%     
%     figure(1)
%     out=imoverlay(I,BWslice);
%     imshow(out)
%     
%   
%     
% 
%     pause;
% end


end


function [maskBkgrnd] = findSkewMask(V,erosionVal)

% get image size and initialize the mask
[sx,sy,sz] = size(V);
maskImage = false(sx,sy,sz);

% since the deconvolved image has a bright edge effect before the black
% region we will essentially "erode" the image field by increasing the
% background mask

% loop over each z-layer
for z=1:sz
    
    colvecImage = ~all(V(:,:,z)==0); % logical vector of which columns the image is in
    
    
    startColImage = find(colvecImage,1,'first');
    stopColImage  = find(colvecImage,1,'last');
    colvec = (startColImage+erosionVal):(stopColImage-erosionVal);
    maskImage(:,colvec,z) = true;

end

maskBkgrnd = ~maskImage;
end


function [Verts,mask] = findSurfaceVerts(BW)

[sx,sy,sz] = size(BW);
[X,Y,Z] = ndgrid(1:sx,1:sy,1:sz);

BW(:,:,1:5) = false;

mask = BW & cumsum(BW,3) == 1;


% the following 3 lines are to remove parts of the mask near edges which are
% unreliable. The X and Y values that I chose were based on one movie. In
% the Y direction (left to right, or column direction) after 550 the mask
% matched the deskew boundary. In the X direction, edge detection didn't
% perform as well near both ends.
mask(Y>550) = false; % 550
mask(X<100)  = false; %75
mask(X>800) = false; %820




Verts = [X(mask),Y(mask),Z(mask)];

end

