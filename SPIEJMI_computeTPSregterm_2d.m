function [TPS,TPSder] = SPIEJMI_computeTPSregterm_2d(points,Tpoints,v,alpha)
    % computes the TPS regularizer and partial derivatives with respect to
    % each RBF transformation parameter
    
    if alpha == 0
        TPS = 0; TPSder = 0;
    else
        num_LM = size(Tpoints,1);
        [sq_dist,X_diff,Y_diff] = dist_arraytoarray(points,Tpoints);
        zero_loc = find(sq_dist==0);
        first_add = 1/log(10);
        second_add_denom = log(10)*sq_dist;
        second_add_num_x = 2*(X_diff).^2; second_add_num_y = 2*(Y_diff).^2;
        third_add = 2*log10(sq_dist.^(0.5));
        
        oneone = first_add + second_add_num_x./second_add_denom + third_add;
        twotwo = first_add + second_add_num_y./second_add_denom + third_add;
        onetwo = 2*(X_diff.*Y_diff)./second_add_denom;
        oneone(zero_loc) = 0; twotwo(zero_loc) = 0; onetwo(zero_loc) = 0;
        
        num_dimpars = num_LM + 3;
        firstdimLMpar = v(1:num_LM)'; seconddimLMpar = v(num_dimpars+1:num_dimpars+num_LM)';
        
        c1rep = repmat(firstdimLMpar,[size(points,1),1]);
        c2rep = repmat(seconddimLMpar,[size(points,1),1]);
        
        oneone_1 = c1rep.*oneone; twotwo_1 = c1rep.*twotwo; onetwo_1 = c1rep.*onetwo;
        oneone_2 = c2rep.*oneone; twotwo_2 = c2rep.*twotwo; onetwo_2 = c2rep.*onetwo;
        
        oneone_1 = sum(oneone_1,2); twotwo_1 = sum(twotwo_1,2); onetwo_1 = sum(onetwo_1,2); % build the entries of the symmetric hessian matrices for each pixel in dim 1
        oneone_2 = sum(oneone_2,2); twotwo_2 = sum(twotwo_2,2); onetwo_2 = sum(onetwo_2,2);
        
        hess_1 = oneone_1.^2 + 2*onetwo_1.^2 + twotwo_1.^2; % compute the inner prod of hessian with itself for each pixel in first dimension
        hess_2 = oneone_2.^2 + 2*onetwo_2.^2 + twotwo_2.^2; % compute the inner prod of hessian with itself for each pixel in second dimension
        
        hess_1 = sum(hess_1); % summing over all pixels
        hess_2 = sum(hess_2); % summing over all pixels
        
        TPS = alpha*(hess_1 + hess_2);
        
        %% compute the partial derivatives with respect to every parameter of
        % the RBF transformation
        TPSder_1stadd_1 = 2*repmat(oneone_1,[1,num_LM]); TPSder_1stadd_1 = TPSder_1stadd_1.*oneone;
        TPSder_2ndadd_1 = 4*repmat(onetwo_1,[1,num_LM]); TPSder_2ndadd_1 = TPSder_2ndadd_1.*onetwo;
        TPSder_3rdadd_1 = 2*repmat(twotwo_1,[1,num_LM]); TPSder_3rdadd_1 = TPSder_3rdadd_1.*twotwo;
        
        TPSder_1stadd_2 = 2*repmat(oneone_2,[1,num_LM]); TPSder_1stadd_2 = TPSder_1stadd_2.*oneone;
        TPSder_2ndadd_2 = 4*repmat(onetwo_2,[1,num_LM]); TPSder_2ndadd_2 = TPSder_2ndadd_2.*onetwo;
        TPSder_3rdadd_2 = 2*repmat(twotwo_2,[1,num_LM]); TPSder_3rdadd_2 = TPSder_3rdadd_2.*twotwo;
        
        partials_c1 = sum(TPSder_1stadd_1 + TPSder_2ndadd_1 + TPSder_3rdadd_1);
        partials_c2 = sum(TPSder_1stadd_2 + TPSder_2ndadd_2 + TPSder_3rdadd_2);
        
        TPSder = [partials_c1 0 0 0 partials_c2 0 0 0]';
        normTPSder = norm(TPSder);
        if normTPSder>0
            TPSder = alpha*(TPSder/normTPSder);
        end
    end
end

function [sq_dist,X_diff,Y_diff] = dist_arraytoarray(Tpoints,Rpoints)
    % compute distances between Tpoints and every Rpoint
    num_R = size(Rpoints,1); num_T = size(Tpoints,1); 
    
    Txrep = repelem(Tpoints(:,1),num_R); Tyrep = repelem(Tpoints(:,2),num_R); 
    Rxstacked = repmat(Rpoints(:,1),[num_T 1]); Rystacked = repmat(Rpoints(:,2),[num_T 1]);
    
    sq_dist = (Txrep - Rxstacked).^2 +(Tyrep - Rystacked).^2; sq_dist = reshape(sq_dist,num_R,num_T)';
    
    X_diff= Txrep - Rxstacked; Y_diff = Tyrep - Rystacked; 
    X_diff = reshape(X_diff,num_R,num_T)';
    Y_diff = reshape(Y_diff,num_R,num_T)';
end