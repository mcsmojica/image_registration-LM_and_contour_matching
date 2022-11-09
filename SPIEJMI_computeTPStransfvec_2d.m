function [Xt, Yt] = SPIEJMI_computeTPStransfvec_2d(points,landmarks,v,m)
    % transforms grid (points) using the parameters vec v
    % k = number of landmarks
    % points means the grid coordinates [X(:) Y(:) Z(:)] 
    % v = [c1^1 c2^1 ... ck^1 w0^1 w1^1 w2^1 w3^1 ...
    % ...  c1^2 c2^2 ... ck^2 w0^2 w1^2 w2^2 w3^2 ...
    % ...  c1^3 c2^3 ... ck^3 w0^3 w1^3 w2^3 w3^3]'
    % landmarks -> the radial basis functions are 'centered' at the
    % landmarks
    
    num_LM = size(landmarks,1); num_dimpar = num_LM + 3;
    num_gridpoints = size(points,1);
    firstdimpar = v(1:num_dimpar)'; seconddimpar = v(num_dimpar+1:2*num_dimpar)';
    [dist_arr,min_rowval,mindist_ind] = SPIEJMI_computedistfromRBFcenters_2d(points,landmarks);
    zero_loc = find(dist_arr==0); logdist = log(dist_arr); logdist(zero_loc) = 0;
    base_mat = [(dist_arr.^2).*logdist ones(num_gridpoints,1) points];
    
    firstdimpar_rep = repmat(firstdimpar,[num_gridpoints 1]); seconddimpar_rep = repmat(seconddimpar,[num_gridpoints 1]);
    y1 = sum(firstdimpar_rep.*base_mat,2); y2 = sum(seconddimpar_rep.*base_mat,2); 
    
    Xt = reshape(y1,m); Yt = reshape(y2,m);
end