function [POI_coords,ordered_POI_coords] = SPIE_JMI_hand_min_intangle(edgecoords,POI_no)
    % find the POI_no points of max curvature given the ordered set of
    % pixels tracing the contour of the image
    
    no_edgecoords = size(edgecoords,1);
    
    % angle_window measurement size tells us how many pixels before and
    % after the current pixel shhould be used to measure the interior angle
    anglewindow = 30;
    
    % windowsize denotes the number of pixels to the left and right of the
    % current POI to be taken out of the running to be the next POI
    windowsize = ceil(no_edgecoords/(4*POI_no)); % (max(8,POI_no)*POI_no)
    
    POI_indlist = [];
    edgecoords_iplus1 = circshift(edgecoords,[-anglewindow 0]);
    edgecoords_iminus1 = circshift(edgecoords,[anglewindow 0]);
    
    % compute interior angle at each point
    diffvec_iplus1 = normr(edgecoords_iplus1-edgecoords); diffvec_iminus1 = normr(edgecoords_iminus1-edgecoords);
    phi = dot(diffvec_iplus1,diffvec_iminus1,2);
    phi = acos(phi);
 
    % find the point of max curvature. next, remove it from the list, along
    % with the points that are closest to it that may also be candidates
    % for max_POIs but should not be considered (because they are within
    % the same neighborhood of the current max)
    for i = 1:POI_no
        % identify the point along the edge with the ith highest curvature
        [m,maxcurv_ind] = min(phi); % [m,maxcurv_ind] = max(phi);
        POI_indlist = [POI_indlist;maxcurv_ind];
        
        % remove the points in the windowsize-neighborhood of the current
        % POI from the list possible POI candidates
        remove_start = maxcurv_ind-windowsize; 
        remove_end = maxcurv_ind+windowsize;  
        
        % reduce the curvature at the pixels to the left of current POI to
        % remove them from list of possible POI candidates
        if remove_start<1
            remove_start = mod(remove_start,no_edgecoords);
            phi(remove_start:end) = 100; 1e-10*phi(remove_start:end);
            phi(1:maxcurv_ind) = 100; 1e-10*phi(1:maxcurv_ind);
        elseif remove_start>1
            phi(remove_start:maxcurv_ind) = 100; 1e-10*phi(remove_start:maxcurv_ind);
        end
        
        % reduce the curvature at the pixels to the right of current POI to
        % remove them from list of possible POI candidates
        if remove_end<1
            remove_end = mod(remove_end,no_edgecoords);
            phi(maxcurv_ind:end) = 100; 1e-10*phi(maxcurv_ind:end);
            phi(1:remove_end) = 100; 1e-10*phi(1:remove_end);
        elseif remove_end>1
            phi(maxcurv_ind:remove_end) = 100; 1e-10*phi(maxcurv_ind:remove_end);
        end

    end
        
    ordered_POI_indlist = sort(POI_indlist);
    POI_coords = edgecoords(POI_indlist,:);
    ordered_POI_coords = edgecoords(ordered_POI_indlist,:);
end