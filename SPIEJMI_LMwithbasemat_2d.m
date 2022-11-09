function der = SPIEJMI_LMwithbasemat_2d(transfTpoints,Tpoints,Rpoints,base_mat,alpha)
    % computes the derivative wrt TPS transf parameters of the LM distance
    % term sum_j ||theta_i(tj)-rj_i||, i = 1,2; j = 1,2,...,no.landmarks
    if alpha == 0
        der = 0;
    else
        num_LM = size(Tpoints,1);

        xdiff = 2*(transfTpoints(:,1) - Rpoints(:,1)); ydiff = 2*(transfTpoints(:,2) - Rpoints(:,2));

        partial1stdim = bsxfun(@times,base_mat,xdiff); partial2nddim = bsxfun(@times,base_mat,ydiff);
        der = [partial1stdim partial2nddim];
        der = sum(der,1)';
        der = alpha*der;
    end
end