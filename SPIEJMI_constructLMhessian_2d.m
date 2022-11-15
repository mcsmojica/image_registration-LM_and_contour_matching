function H = SPIEJMI_constructLMhessian_2d(Tpoints,alpha)
    num_LM = size(Tpoints,1);
    [dist_arr,min_rowval,mindist_ind] = SPIEJMI_computedistfromRBFcenters_2d(Tpoints,Tpoints);
    zero_loc = find(dist_arr==0); logdist = log(dist_arr); logdist(zero_loc) = 0;
    basemat_TPS = [(dist_arr.^2).*logdist ones(num_LM,1) Tpoints];
    % basemat_rho = basemat_TPS(:,1:num_LM);
    H_base = zeros(num_LM+3);
    H_base0 = H_base;
    
    % build the lower triangular part first
    for col = 1:num_LM+3
        for row = 1:col
          % elements of the same column j have a common multiplier center:
          % ||t_j-t_col|| for j=1,...,num_LM+3
          H_base(row,col) = 2*dot(basemat_TPS(:,col),basemat_TPS(:,row));
        end
    end
    H_base = (H_base+H_base') - eye(size(H_base,1)).*diag(H_base);
    H = [H_base H_base0;H_base0 H_base];
end