function[threshVector]=thresholdVector(reference,data,binvector,ierror)
%thresholdVector thresholds the data in the data vector according to the
%specified bins in the binvector
%
% SYNOPSIS [threshVector]=thresholdVector(reference,data,binvector);
%
% INPUT:    reference: reference data for the binvector (if data is e.g.
%                   intensity versus distance, and has to be thresholded
%                   for specified distance bins, then reference contains
%                   the original distances corresponding to all the
%                   intensity points)
%           data: vector with the data to be averaged for all the bins (in the
%                   above example, the vector containing the intensity
%                   data); this vector needs to have the same length as the
%                   reference vector!!
%           binvector: the vector containing the (distance) bins for the
%                   thresholding; the data in the bin are considered the
%                   limits of the bins, not the center points!!
%           ierror: include error (optional)
%                   0 (default) no change
%                   1 include error and n as second and third row
% OUTPUT:   threshVector: thresholded vector; the size of this vector is
%                   2xbn-1, the number of bins is 1 less than the number of
%                   points in binvector
%
% c: 04/27/2006 Dinah Loerke

% set default error, if error was inputted then change it
err = 0;
if nargin>3
    if ierror==1
        err=1;
    end
end

% convert reference and data to double class
reference = double(reference(:));
data = double(data(:));

% initialize output vector
threshVector = zeros(2,length(binvector)-1);
threshVector(:) = nan;
threshVector(1,:) = binvector(1:length(binvector)-1);

% throw an error if reference and data vector are different lengths
if (length(reference)~=length(data))
    error('reference and data vector of unequal length');
end

if length(reference)>0

    % check input vectors for orientation; they should be (nx1)
    % throw errors if the vector format is incorrect
    [checksizx,checksizy]=size(reference);

    if (min(checksizx,checksizy)>1)
        error('thresholdVector doesn''t know how to handle this vector format');
    end

    if ( checksizx < checksizy )
        referenceS = reference';
    else
        referenceS = reference;
    end

    [checksizx,checksizy]=size(data);

    if (min(checksizx,checksizy)>1)
        error('thresholdVector doesn''t know how to handle this vector format');
    end

    if ( checksizx < checksizy )
        dataS = data';
    else
        dataS = data;
    end
    
    usepos = find( isfinite(referenceS) & isfinite(dataS) );

    % sort data data according to reference data
    sortMat = sortrows([referenceS(usepos) dataS(usepos)],1);
    referenceSortvec = sortMat(:,1);
    dataSortvec = sortMat(:,2);

    % identify cutoff point
    maxpos = max(find(referenceSortvec<max(binvector)));
    
    %make provision for the case that maxpos = [], else there can be an
    %error message when the code is run
    if isempty(maxpos)
        return
    else
        
        binPosVector = 1:length(binvector);
        binPosVector(length(binvector)) = maxpos;

        %loop over bins, from end to start
        for b = (length(binvector)-1):-1:1
        
        % find position in sorted reference vector where the value is smaller
        % than the relevant bin value
        % NOTE: a previous version of this function - embedded into
        % correlationfromMPM - always had one [0 0] value in the reference and
        % data vector; since this function is supposed to be as general as
        % possible, and should be usable for negative data as well, some
        % changes have to be made for the case that NO data fulfill the
        % criterion, yielding an empty find vector    
            findVec = find( referenceSortvec(1:binPosVector(b+1)) < binvector(b) );
            if (length(findVec)>0)
                binPosVector(b)=max( findVec );
            else
                binPosVector(b)=0;
            end
    
        % if a smaller value than the previous exists (i.e. if there are any
        % values falling inside this bin), then the data are averaged
            if ( binPosVector(b)<binPosVector(b+1) )
                %threshVector(2,b) = nanmean( dataSortvec(binPosVector(b)+1:binPosVector(b+1)) );
                threshVector(2,b) = mean( dataSortvec(binPosVector(b)+1:binPosVector(b+1)),'omitnan' );
                if err==1
                    %threshVector(3,b) = nanstd( dataSortvec(binPosVector(b)+1:binPosVector(b+1)) );
                    threshVector(3,b) = std( dataSortvec(binPosVector(b)+1:binPosVector(b+1)),[],'omitnan' );
                    threshVector(4,b) = length( [binPosVector(b)+1:binPosVector(b+1)] );
                end

                % if the value has just dropped to zero, then there are no more
                % points left for the smaller bins, so the function can stop
                % calculating now
                if (binPosVector(b)==0)
                    return
                end

            end % of if

        end  %of for b (bins)
        
    end % of if maxpos is defined
    
end % of if reference vector isn't empty


end % of function