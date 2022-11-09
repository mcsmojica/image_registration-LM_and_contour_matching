function [CM_fcnval,CM_der,partialv,CM_hess] = SPIEJMI_CMvalderHess_2d(Txi,transfTxi,Rxi,Tpoints,alpha)
    looping = 1;
    
    num_LM = size(Tpoints,1);
    
    [DP,partialv,transfTxi_mag,R_unitvec,transfTxi_vec] = SPIE_getCM_latitudecomponents_2d(Tpoints,Rxi,Txi,transfTxi,[],[],[],[],[],looping);
    CM_fcnval = 0.5*alpha*sum(1 - DP.^2);
    
    partialzero = zeros(size(partialv));
    
    LDR1 = [bsxfun(@times,partialv,R_unitvec(:,1)) partialzero]; LDR2 = [partialzero bsxfun(@times,partialv,R_unitvec(:,2))];
    LDR = (transfTxi_mag.^(-0.5)).*(LDR1 + LDR2);
    
    Right = -dot(transfTxi_vec,R_unitvec,2).*(transfTxi_mag.^(-1.5));
    RDL1 = [bsxfun(@times,partialv,transfTxi_vec(:,1)) partialzero]; RDL2 = [partialzero bsxfun(@times,partialv,transfTxi_vec(:,2))];
    RDL = bsxfun(@times,RDL1+RDL2,Right);
    
    der = -bsxfun(@times,LDR+RDL,DP);
    der = sum(der)';
    CM_der = alpha*der; % alpha*(der/norm(der)); %
    
    %% compute CM Hessian
    % product of non-unit transfTxivec and unitRxivec
    R = dot(transfTxi_vec,R_unitvec,2); 
       
    % get the squared magnitude and the components
    magsquared_transfTxi_vec = sum(transfTxi_vec.^2,2);
    % first column/component of vector of transformed Txi vectors
    % (square1 in code)
    transfTxi_vec1st = transfTxi_vec(:,1);
    % second column/component of vector of transformed Txi vectors
    % (square2 in code)
    transfTxi_vec2nd = transfTxi_vec(:,2);
    
    H_base = zeros(num_LM+3);
    H_base0 = H_base;
    
    partialv_trans = partialv';
    Rsquaredrow = (R.^2)';
    magsquared_transfTxi_vec_row = magsquared_transfTxi_vec';
    
    %% H1    
    Rcol1overpartialv = bsxfun(@times,partialv,R_unitvec(:,1));
    Rcol2overpartialv = bsxfun(@times,partialv,R_unitvec(:,2));
    
    H1_left = []; H1_right = [];
    for i = 1:num_LM+3 % sweep through the columns of H1
        
        % upper left part of H1
        H1_upperleft_factor1_i = Rcol1overpartialv(:,i)'; H1_upperright_factor1_i = Rcol2overpartialv(:,i)';
        H1_upperleft_factor1_i = repmat(H1_upperleft_factor1_i,[num_LM+3,1]); H1_upperright_factor1_i = repmat(H1_upperright_factor1_i,[num_LM+3,1]);
        
        H1_upperleft_factor2a = -2*R'.*(magsquared_transfTxi_vec_row.^(-2)).*transfTxi_vec1st';
        H1_upperleft_factor2a = bsxfun(@times,partialv_trans,H1_upperleft_factor2a); H1_upperright_factor2a = H1_upperleft_factor2a;
        
        H1_upperleft_factor2b = (magsquared_transfTxi_vec_row.^(-1));
        H1_upperleft_factor2b = bsxfun(@times,LDR1(:,1:num_LM+3)',H1_upperleft_factor2b); H1_upperright_factor2b = H1_upperleft_factor2b;
        
        H1_upperleft = H1_upperleft_factor1_i.*(H1_upperleft_factor2a + H1_upperleft_factor2b);
        H1_upperright = H1_upperright_factor1_i.*(H1_upperright_factor2a + H1_upperright_factor2b);
        
        % lowerleft half of H1
        H1_lowerleft_factor1_i = H1_upperleft_factor1_i; H1_lowerright_factor1_i = H1_upperright_factor1_i;
        
        H1_lowerleft_factor2a = -2*R'.*(magsquared_transfTxi_vec_row.^(-2)).*transfTxi_vec2nd';
        H1_lowerleft_factor2a = bsxfun(@times,partialv_trans,H1_lowerleft_factor2a); H1_lowerright_factor2a = H1_lowerleft_factor2a;
        
        H1_lowerleft_factor2b = (magsquared_transfTxi_vec_row.^(-1));
        H1_lowerleft_factor2b = bsxfun(@times,LDR2(:,num_LM+3+1:end)',H1_lowerleft_factor2b); H1_lowerright_factor2b = H1_lowerleft_factor2b;
        
        H1_lowerleft = H1_lowerleft_factor1_i.*(H1_lowerleft_factor2a + H1_lowerleft_factor2b); H1_lowerright = H1_lowerright_factor1_i.*(H1_lowerright_factor2a + H1_lowerright_factor2b);
        
        H1leftcol = [H1_upperleft;H1_lowerleft]; H1rightcol = [H1_upperright;H1_lowerright];
        H1_left = [H1_left sum(H1leftcol,2)]; H1_right = [H1_right sum(H1rightcol,2)];
    end
    
    H1 = [H1_left H1_right];
    
    
    %% H2
    H_zerohalf = zeros(2*(num_LM+3),num_LM+3);
    H2_left = [];
    
    for i = 1:num_LM+3
        % H2 upper half
        H2_upperleft_factor1_i = partialv(:,i)'; 
        H2_upperleft_factor1_i = -repmat(H2_upperleft_factor1_i,[num_LM+3,1]);
        
        H2_upperleft_factor2a = Rsquaredrow.*(magsquared_transfTxi_vec_row.^(-2));
        H2_upperleft_factor2a = bsxfun(@times,partialv_trans,H2_upperleft_factor2a);
        
        H2_upperleft_factor2b = -4*Rsquaredrow.*((transfTxi_vec1st.^2)').*(magsquared_transfTxi_vec_row.^(-3));
        H2_upperleft_factor2b = bsxfun(@times,partialv_trans,H2_upperleft_factor2b);
        
        H2_upperleft_factor2c = 2*(magsquared_transfTxi_vec_row.^(-2)).*(transfTxi_vec1st').*R';
        H2_upperleft_factor2c = bsxfun(@times,LDR1(:,1:num_LM+3)',H2_upperleft_factor2c);
        
        H2_upperleft = H2_upperleft_factor1_i.*(H2_upperleft_factor2a + H2_upperleft_factor2b + H2_upperleft_factor2c);
        
        
        % H2 lower half
        H2_lowerleft_factor1_i = H2_upperleft_factor1_i; 
        
        H2_lowerleft_factor2a = 0;
        
        H2_lowerleft_factor2b = -4*Rsquaredrow.*(transfTxi_vec1st').*(transfTxi_vec2nd').*(magsquared_transfTxi_vec_row.^(-3));
        H2_lowerleft_factor2b = bsxfun(@times,partialv_trans,H2_lowerleft_factor2b);
        
        H2_lowerleft_factor2c = 2*(magsquared_transfTxi_vec_row.^(-2)).*(transfTxi_vec1st').*R';
        H2_lowerleft_factor2c = bsxfun(@times,LDR2(:,num_LM+3+1:end)',H2_lowerleft_factor2c);
        
        H2_lowerleft = H2_lowerleft_factor1_i.*(H2_lowerleft_factor2a + H2_lowerleft_factor2b + H2_lowerleft_factor2c);     
        
        H2leftcol = [H2_upperleft;H2_lowerleft]; 
        H2_left = [H2_left sum(H2leftcol,2)]; 
    end
    
    H2 = [H2_left H_zerohalf];
    
    
    %% H3
    H3_right = [];
    
    for i = 1:num_LM+3
        % H3 upper half
        H3_upperright_factor1_i = partialv(:,i)'; 
        H3_upperright_factor1_i = -repmat(H3_upperright_factor1_i,[num_LM+3,1]);
        
        H3_upperright_factor2a = 0; 
        
        H3_upperright_factor2b = -4*Rsquaredrow.*(transfTxi_vec2nd').*(transfTxi_vec1st').*(magsquared_transfTxi_vec_row.^(-3));
        H3_upperright_factor2b = bsxfun(@times,partialv_trans,H3_upperright_factor2b);
        
        H3_upperright_factor2c = 2*(magsquared_transfTxi_vec_row.^(-2)).*(transfTxi_vec2nd').*R';
        H3_upperright_factor2c = bsxfun(@times,LDR1(:,1:num_LM+3)',H3_upperright_factor2c);
        
        H3_upperright = H3_upperright_factor1_i.*(H3_upperright_factor2a + H3_upperright_factor2b + H3_upperright_factor2c);
        
        
        % H3 lower half
        H3_lowerright_factor1_i = H3_upperright_factor1_i; 
        
        H3_lowerright_factor2a = Rsquaredrow.*(magsquared_transfTxi_vec_row.^(-2));
        H3_lowerright_factor2a = bsxfun(@times,partialv_trans,H3_lowerright_factor2a);
        
        H3_lowerright_factor2b = -4*Rsquaredrow.*((transfTxi_vec2nd.^2)').*(magsquared_transfTxi_vec_row.^(-3));
        H3_lowerright_factor2b = bsxfun(@times,partialv_trans,H3_lowerright_factor2b);
        
        H3_lowerright_factor2c = 2*(magsquared_transfTxi_vec_row.^(-2)).*(transfTxi_vec2nd').*R';
        H3_lowerright_factor2c = bsxfun(@times,LDR2(:,num_LM+3+1:end)',H3_lowerright_factor2c);
        
        H3_lowerright = H3_lowerright_factor1_i.*(H3_lowerright_factor2a + H3_lowerright_factor2b + H3_lowerright_factor2c);     
        
        H3rightcol = [H3_upperright;H3_lowerright]; 
        H3_right = [H3_right sum(H3rightcol,2)]; 
    end
    
    H3 = [H_zerohalf H3_right];
   
    
    CM_hess = -alpha*(H1+H2+H3);    
%     disp('1');
end

function norm_vecs = getnormed_diffvector(vec)
    norm_vecs = diff(vec); norm_vecs = normr(norm_vecs);
end 

function vecprod = columnbymatrix(col,matr)
    % multiply each column of a matrix by the corresponding element of a vector
    x = repmat(col,1,size(matr,2));
    vecprod = x.*matr;
end