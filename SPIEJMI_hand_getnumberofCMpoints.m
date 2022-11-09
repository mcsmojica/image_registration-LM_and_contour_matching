function CM_num = SPIEJMI_hand_getnumberofCMpoints(ordered_POIs,first_arcRxi)
    % determine a proportional number of contour-approximating landmarks given
    % the desired number of contour-approximating landmarks in the first arc
    % (the arc flanked by the first two major landmarks)
    %
    % inputs: 
    % POIs = ordered list of exact landmarks, whose consecutive
    % entries' distances will determine the number of landmarks in each
    % succeeding arc
    % first_arcRxi = array containing info about the desired number of
    % contour-matching points in the first arc of each experiment (i.e.,
    % the number of rows of this array should correspond to the number of
    % experiments to be performed
    
    % output:
    % CMpoints_num = nxnum_POI array containing the number of CM points per
    % arc for n experiments; num_POI = number of arcs = number of exact LMs
    
    num_expts = size(first_arcRxi,1);
    num_POIs = size(ordered_POIs,1);
    
    % get the distance between consecutive landmarks
    R_ordered_POIs_looping = [ordered_POIs;ordered_POIs(1,:)];
    POIdist = sqrt(sum(diff(R_ordered_POIs_looping).^2,2));

    % get the floor ratio of each POI distance with the distance between the
    % first two major POIs
    arclength_ratiotofirstarc = (1/POIdist(1))*POIdist';
    arclength_ratiotofirstarc = repmat(arclength_ratiotofirstarc,[num_expts,1]);
    
    CM_num = round(bsxfun(@times,arclength_ratiotofirstarc,first_arcRxi));
end