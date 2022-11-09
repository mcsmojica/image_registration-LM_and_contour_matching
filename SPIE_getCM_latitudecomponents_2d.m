function [DP,partialv,mag_transfTvec_arr,Runitdiff_arr,transfTdiff_arr] = SPIE_getCM_latitudecomponents_2d(Tpoints,R_clusterpoints_slice,T_clusterpoints_slice,transfT_clusterpoints_slice,DP,partialv,mag_transfTvec_arr,Runitdiff_arr,transfTdiff_arr,looping)
    if looping
        number_clusterpoints = size(R_clusterpoints_slice,1);
        
        % get the "looping" list of unit vectors
        R_h_unitvec_cluster = normr(diff([R_clusterpoints_slice;R_clusterpoints_slice(1,:)])); T_h_unitvec_cluster = normr(diff([T_clusterpoints_slice;T_clusterpoints_slice(1,:)]));
        transfT_h_unitvec_cluster = normr(diff([transfT_clusterpoints_slice;transfT_clusterpoints_slice(1,:)]));
        
        % get the looping list of non-unit vectors
        R_h_vec_cluster = diff([R_clusterpoints_slice;R_clusterpoints_slice(1,:)]); T_h_vec_cluster = diff([T_clusterpoints_slice;T_clusterpoints_slice(1,:)]);
        transfT_h_vec_cluster = diff([transfT_clusterpoints_slice;transfT_clusterpoints_slice(1,:)]);
        
    else
        number_clusterpoints = size(R_clusterpoints_slice,1) - 1;
        % get the "non-looping" list of unit vectors
        R_h_unitvec_cluster = normr(diff(R_clusterpoints_slice)); T_h_unitvec_cluster = normr(diff(T_clusterpoints_slice));
        transfT_h_unitvec_cluster = normr(diff(transfT_clusterpoints_slice));
        
        % get the non-looping list of non-unit vectors
        R_h_vec_cluster = diff(R_clusterpoints_slice); T_h_vec_cluster = diff(T_clusterpoints_slice);
        transfT_h_vec_cluster = diff(transfT_clusterpoints_slice);
    end
    
    % get the dot product of corresponding unit vectors and add the DP
    % to the list of DPs
    DP = [DP;dot(R_h_unitvec_cluster,transfT_h_unitvec_cluster,2)];

    % Compute the distance arrays from the xis (points along the contour of the template image) of every segment from the Tpoints
    % [||x1-t1||  ... ||x1-tk||]
    % [||x2-t1||  ... ||x2-tk||]
    % ...
    % [||xn-t1||  ... ||xn-tk||]
    if looping
        latitude_Txi_looping = [T_clusterpoints_slice;T_clusterpoints_slice(1,:)];
    else
        latitude_Txi_looping = T_clusterpoints_slice;
    end
    
    % construct partialbasematrix/kernel matrix for CM term
    if size(Tpoints,2)==2
        [distfromTpoints,M1,M2] = SPIEJMI_computedistfromRBFcenters_2d(latitude_Txi_looping,Tpoints);
        zero_loc = find(distfromTpoints==0); logdist = log(distfromTpoints); logdist(zero_loc) = 0; TPSrho = (distfromTpoints.^2).*logdist;
        partialv = [TPSrho(2:end,:)-TPSrho(1:end-1,:) zeros(number_clusterpoints,1) T_h_vec_cluster];
    elseif size(Tpoints,2)==3
        [distfromTpoints,M1,M2] = SPIEJMI_computedistfromRBFcenters_3d(latitude_Txi_looping,Tpoints);
        partialv = [partialv;distfromTpoints(2:end,:)-distfromTpoints(1:end-1,:) zeros(number_clusterpoints,1) T_h_vec_cluster];
    end

    % Squared Magnitude vector of NON-unit vectors ||r_{i+1}-r_i||^2 AND ||x_{i+1}-x_i||^2
    mag_transfTvec = sum(transfT_h_vec_cluster.^2,2);
    % Add to array of magnitudes of non-unit transfTvecs
    mag_transfTvec_arr = [mag_transfTvec_arr;mag_transfTvec];
    % Add to array of "looping" list of unit Rvecs
    Runitdiff_arr = [Runitdiff_arr;R_h_unitvec_cluster];
    % Add to array of "looping" list of non-unit transfTvecs
    transfTdiff_arr = [transfTdiff_arr;transfT_h_vec_cluster];
end