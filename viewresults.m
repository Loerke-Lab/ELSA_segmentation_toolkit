function [] = viewresults(data,imageTypeStr,framedelay,tvec)
% viewresults(data,pspeed,tvec,type)
% VIEWRESULTS allows visualization of segmentation results and tracking data

% INPUT:    data: structure that contains lists of image files (as created by the
%                   function ELSA_1loadImageList) and the source (or path) to the
%                   images for each movie. The image file list should be
%                   stored in data.ImageFileList and the path to the images
%                   should be stored in data.Source.
%           framedelay: time interval between frames of playback
%           tvec (optional): time vector specifying the frames to be analyzed
%                   (e.g. [1:20]). If empty or left out tvec will be full
%                   movie.
%           imageTypeStr: string which specifies which type of results to view:
%                   'raw'       - view raw images.
%                   'tracking'  - view tracking labels.
%                   'seg'       - view bw image of segmentation
%                   'overlay'   - view segmentation lines overlaid on images.
%                   'seeds'     - view seeds and mask
% OUTPUT:   no function output - results are displayed on image for all
%           time points 
%
%
% T. Vanderleest 6/29/14

% record original directory (to return to it at the end)
od = cd;

% imshow magnification value
magval = 200;

% close all other figures
close all

% get image file list
list = data.ImageFileList;

% DEFAULT time vector are set to min and max available,
% unless different input is specified
tstart = 1;
tend = length(list);
tvecDefault = tstart:tend;
if nargin<4 || isempty(tvec)
    tvec = tvecDefault;
end



% loop over timepoints
for ti = 1:length(tvec)
    
    % specify t inside the loop
    t = tvec(ti);
    
    
    switch imageTypeStr
        case 'raw'
            % load image from imageFileList and convert to double
            image = imread(list{t});
            image = mat2gray(image);
            imshow(image,'InitialMagnification',magval);
            
        case 'tracking'
            % load segmentation results for current time point
            cd(data.Source);
            cd('SegmentationData');
            cframefoldername = sprintf('frame%04d',t);
            cd(cframefoldername);
            
            loadseg = load('ImageBWlabel_trackT.mat');
            BWlabel = loadseg.ImageBWlabel_trackT;
            BWlabel(BWlabel==-1) = 0;
            rgprops = regionprops(BWlabel,'Centroid');

            imshow(BWlabel,'InitialMagnification',magval)
            labels = unique(BWlabel(:));
            labels(labels == 0) = [];
            for i=1:length(labels)
                li = labels(i);
                C =rgprops(li).Centroid;
                text(C(1),C(2),num2str(li))
            end
            
        case 'seg'
            % load segmentation results for current time point
            cd(data.Source);
            cd('SegmentationData');
            cframefoldername = sprintf('frame%04d',t);
            cd(cframefoldername);
            
            loadseg = load('ImageSegment.mat');
            seg = loadseg.ImageSegment;
            imshow(seg,'InitialMagnification',magval)
            
        case 'overlay'
            % load segmentation results for current time point
            cd(data.Source);
            cd('SegmentationData');
            cframefoldername = sprintf('frame%04d',t);
            cd(cframefoldername);
            
            image = imread(list{t});
            image = mat2gray(image);
            if exist('ImageSegment.mat') == 2
                loadseg = load('ImageSegment');
                seg = loadseg.ImageSegment;
                imoverlay(image,seg==0);
            else
                 imshow(image,'InitialMagnification',magval);
            end
            
        case 'seeds'
             % load segmentation results for current time point
            cd(data.Source);
            cd('SegmentationData');
            cframefoldername = sprintf('frame%04d',t);
            cd(cframefoldername);
            
            image = imread(list{t});
            image = mat2gray(image);
            if exist('seeds.mat') == 2
                seeds = load('seeds.mat').seeds;
                mask = load('mask.mat').mask;
                out=imoverlay(image,seeds,[0 1 0]);
                imoverlay(out,mask,[0 0 1]);
            else
                 imshow(image,'InitialMagnification',magval);
            end
                
            
        otherwise
            % throw a warning if the type of results to view is not correct
            warning('Unexpected type');
            
            
    end
    
        

    title(['Time frame: ',num2str(t)],'FontSize',16)
    
    % pause between time points
    if ~isempty(framedelay)
        pause(framedelay);
    else
        pause;
    end
    
end

% return to original directory
cd(od)

end % function



function out = imoverlay(in, mask, color)
%IMOVERLAY Create a mask-based image overlay.
%   OUT = IMOVERLAY(IN, MASK, COLOR) takes an input image, IN, and a binary
%   image, MASK, and produces an output image whose pixels in the MASK
%   locations have the specified COLOR.
%
%   IN should be a grayscale or an RGB image of class uint8, uint16, int16,
%   logical, double, or single.  If IN is double or single, it should be in
%   the range [0, 1].  If it is not in that range, you might want to use
%   mat2gray to scale it into that range.
%
%   MASK should be a two-dimensional logical matrix.
%
%   COLOR should be a 1-by-3 vector of values in the range [0, 1].  [0 0 0]
%   is black, and [1 1 1] is white.
%
%   OUT is a uint8 RGB image.
%
%   Examples
%   --------
%   Overlay edge detection result in green over the original image.
%       
%       I = imread('cameraman.tif');
%       bw = edge(I, 'canny');
%       rgb = imoverlay(I, bw, [0 1 0]);
%       imshow(rgb)
%
%   Treating the output of peaks as an image, overlay the values greater than
%   7 in red.  The output of peaks is not in the usual grayscale image range
%   of [0, 1], so use mat2gray to scale it.
%
%       I = peaks;
%       mask = I > 7;
%       rgb = imoverlay(mat2gray(I), mask, [1 0 0]);
%       imshow(rgb, 'InitialMagnification', 'fit')

%   Steven L. Eddins, The MathWorks, Inc.
%   $Revision: 1.2 $  $Date: 2007/08/15 13:18:08 $

% If the user doesn't specify the color, use red.
DEFAULT_COLOR = [1 0 0];
if nargin < 3
    color = DEFAULT_COLOR;
end

% Make the uint8 the working data class.  The output is also uint8.
in_uint8 = im2uint8(in);
color_uint8 = im2uint8(color);

% Initialize the red, green, and blue output channels.
if ndims(in_uint8) == 2
    % Input is grayscale.  Initialize all output channels the same.
    out_red   = in_uint8;
    out_green = in_uint8;
    out_blue  = in_uint8;
else
    % Input is RGB truecolor.
    out_red   = in_uint8(:,:,1);
    out_green = in_uint8(:,:,2);
    out_blue  = in_uint8(:,:,3);
end

% Replace output channel values in the mask locations with the appropriate
% color value.
out_red(mask)   = color_uint8(1);
out_green(mask) = color_uint8(2);
out_blue(mask)  = color_uint8(3);

% Form an RGB truecolor image by concatenating the channel matrices along
% the third dimension.
out = cat(3, out_red, out_green, out_blue);

% Modification: Following lines added by Tim Vanderleest on 11/15/12
if nargout == 0
    imshow(out,'InitialMagnification',200)
end



end



