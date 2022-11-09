function [dist_TR,min_rowval,mindist_ind] = SPIEJMI_computedistfromRBFcenters_3d(points,RBFcenters)
    % compute distances between every points from every RBFcenterpoint
    % (usually the landmarks)
    num_RBFcenters = size(RBFcenters,1); num_points = size(points,1); 
    
    points_xrep = repelem(points(:,1),num_RBFcenters); points_yrep = repelem(points(:,2),num_RBFcenters); points_zrep = repelem(points(:,3),num_RBFcenters);
    Rxstacked = repmat(RBFcenters(:,1),[num_points 1]); 
    Rystacked = repmat(RBFcenters(:,2),[num_points 1]);
    Rzstacked = repmat(RBFcenters(:,3),[num_points 1]);
    
    dist_TR = sqrt((points_xrep - Rxstacked).^2 + (points_yrep - Rystacked).^2 + (points_zrep - Rzstacked).^2); 
    dist_TR = reshape(dist_TR,num_RBFcenters,num_points)';
    [min_rowval,mindist_ind] = min(dist_TR,[],2);   
end