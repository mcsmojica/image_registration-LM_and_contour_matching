function [v_small,f_LM_small,f_CM_small,Xt,Yt,transfTpoints,transfTxi,Ty] = SPIEJMI_LCM_2d(dataR,dataT,Rpoints,Tpoints,Rxi,Txi,LMalpha,CMalpha,titlestr,save_filename,save_file,save_images)
    TPSalpha=0;

    %% Set up grid.
    % Instensities are cell-centered; Landmarks are based on pixel values, so
    % translate the landmarks 0.5units down in order for the landmark locations
    % to coincide with some of the grid points and for each landmark to have an
    % associated intensity value
    [m,n] = size(dataR); num_pixels = m*n;
    x = 0.5:1:n-0.5; y = 0.5:1:m-0.5;
    [X,Y] = meshgrid(x,y); points = [X(:) Y(:)];

    num_LM = size(Rpoints,1); num_xi = size(Rxi,1);

    %% Initialization
    % Initialize stopping parameters and learning rate
    tol = 1e-8; rel_imp = 1;
    if CMalpha<1
        max_iter = 100000;
    else
        max_iter = 250000;
    end
    iter = 0;
    %     LMalpha = 1;
    vals = [];

    %% Specify initial guess for vector of TPS transf parameters
    v = (1e-8)*zeros(2*num_LM+6,1); v(num_LM+2) = 1; v(2*num_LM+3+3) = 1;
    [transfTpoints(:,1), transfTpoints(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints,Tpoints,v,size(Tpoints(:,1)));
    [transfTxi(:,1), transfTxi(:,2)] = SPIEJMI_computeTPStransfvec_2d(Txi,Tpoints,v,size(Txi(:,1)));

    %% Get the base matrix for LM term partial derivative and second partial der(does not change per iteration, except for the first multiplier that depends on the transfTpoints)
    [dist_arr,min_rowval,mindist_ind] = SPIEJMI_computedistfromRBFcenters_2d(Tpoints,Tpoints);
    zero_loc = find(dist_arr==0); logdist = log(dist_arr); logdist(zero_loc) = 0;
    partialbasemat_fLM = [(dist_arr.^2).*logdist ones(num_LM,1) Tpoints];
    H_LM = SPIEJMI_constructLMhessian_2d(Tpoints,LMalpha); % zeros(num_LM,num_LM);

    % Jacobian of the constraints
    J_base1 = zeros(3,num_LM+3); J_base0 = J_base1;
    J_base1(1,1:num_LM) = 1; J_base1(2:3,1:num_LM) = Tpoints';
    J = [J_base1 J_base0;J_base0 J_base1];
    Newt_mat0 = zeros(size(J,1),size(J,1));

    f_LM = SPIEJMI_computeLMfcnval_2d(transfTpoints,Rpoints,LMalpha);
    [f_CM,fprime_CMTPS,partialbasemat_fCM,H_CM] = SPIEJMI_CMvalderHess_2d(Txi,transfTxi,Rxi,Tpoints,CMalpha);

    f_new = f_LM + f_CM;

    f_small = f_new; transfTpoints_small = transfTpoints; transfTxi_small = transfTxi; v_small = v;
    f_LM_small = f_LM; f_CM_small = f_CM;
    no_imp = 0;
    while (rel_imp>tol & iter<max_iter & no_imp<=5000)
        fprime_LMTPS = SPIEJMI_LMwithbasemat_2d(transfTpoints,Tpoints,Rpoints,partialbasemat_fLM,LMalpha);
        if CMalpha>0
            [f_CM,fprime_CMTPS,partialbasemat_fCM,H_CM] = SPIEJMI_CMvalderHess_2d(Txi,transfTxi,Rxi,Tpoints,CMalpha);
        else
            f_CM = 0; fprime_CMTPS = 0;
        end

        fprime_old = fprime_LMTPS + fprime_CMTPS;
        f_old = f_new;

        H = LMalpha*H_LM + H_CM;
        h = SPIEJMI_computeTPSconstraint_2d(v,Tpoints);
        Newt_matA = [H J';J Newt_mat0]; Newt_matb = [-fprime_old;-h];
        Newt_matx = Newt_matA\Newt_matb;
        delta_v = Newt_matx(1:size(J,2),1);

        v = v + delta_v;
        [transfTpoints(:,1), transfTpoints(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints,Tpoints,v,size(Tpoints(:,1)));
        [transfTxi(:,1), transfTxi(:,2)] = SPIEJMI_computeTPStransfvec_2d(Txi,Tpoints,v,size(Txi(:,1)));
        fLM_old = f_LM; fCM_old = f_CM;

        % Compute new function values
        f_LM = SPIEJMI_computeLMfcnval_2d(transfTpoints,Rpoints,LMalpha);
        if CMalpha>0
            [f_CM,fprime_CMTPS,partialbasemat_fCM,H_CM] = SPIEJMI_CMvalderHess_2d(Txi,transfTxi,Rxi,Tpoints,CMalpha);
        else
            f_CM = 0; fprime_CMTPS = 0;
        end

        f_new = f_LM + f_CM;

        iter = iter + 1; rel_imp = (abs(f_old - f_new)/f_old);

        if f_new < f_small
            add_row = [iter f_new rel_imp]
            vals = [vals;add_row];
            decrease_fval = 0;
            if (f_LM<=fLM_old & f_CM<=fCM_old)% & f_TPS<=fTPS_old)
                disp('Every term is getting smaller...');
            else
                disp('NOT every term is getting smaller, but the total is decreasing...');
            end
            f_small = f_new;
            f_LM_small = f_LM; f_CM_small = f_CM;
            v_small = v;
        else
            iter
            no_imp = no_imp + 1;
        end
    end
    [Xt, Yt] = SPIEJMI_computeTPStransfvec_2d(points,Tpoints,v_small,size(dataT));
    [transfTpoints(:,1), transfTpoints(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints,Tpoints,v_small,size(Tpoints(:,1)));
    [transfTxi(:,1), transfTxi(:,2)] = SPIEJMI_computeTPStransfvec_2d(Txi,Tpoints,v_small,size(Txi(:,1)));
    Ty = griddata(Xt,Yt,dataT,X,Y,'cubic'); NaN_loc = find(isnan(Ty) == 1); Ty(NaN_loc) = mode(dataT(:));

    [prereg_dice,prereg_jaccard] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(dataT)); 
    [postreg_dice,postreg_jaccard] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(Ty));

    h = figure('units','normalized','outerposition',[0 0 1 1]); set(h, 'Visible', 'off');
    TitleStr = strcat('Pre-Reg Similarity:',num2str(prereg_dice,4), '; Post-Reg Similarity:', num2str(postreg_dice,4));

    if save_file==1
        save(strcat(save_filename,'.mat'));
    end
    if CMalpha==0
        sgtitle(strcat('TPS: ',titlestr, '--- LM',num2str(LMalpha),' CM',num2str(CMalpha),' TPS',num2str(TPSalpha), '---', TitleStr));
    else
        sgtitle(strcat('LMCM: ',titlestr, '--- LM',num2str(LMalpha),' CM',num2str(CMalpha),' TPS',num2str(TPSalpha), '---', TitleStr));
    end

    subplot(2,3,1); imshow(dataR,[]); title('R');
    subplot(2,3,2); imshow(dataT,[]); title('T');
    subplot(2,3,3); imshow(abs(dataR-dataT),[]); title('|R - T|');
    subplot(2,3,4); imshow(abs(dataR-Ty),[]); title('|R - T[y]|');
    rotatedRpoints = Rpoints; rotatedTpoints = Tpoints; rotatedtransfTpoints = transfTpoints;
    subplot(2,3,5); imshow(Ty,[]); title('T[y]'); hold on; plot(rotatedRpoints(:,1),rotatedRpoints(:,2),'m.','MarkerSize',25); hold on; plot(rotatedTpoints(:,1),rotatedTpoints(:,2),'b.','MarkerSize',25); plot(rotatedtransfTpoints(:,1), rotatedtransfTpoints(:,2), 'g.', 'MarkerSize',15);

    legend('Rpoints','Tpoints','Transf Tpoints','Location','Best');
    subplot(2,3,6); JanModersitzki_plotGrid([Xt(:);Yt(:)],[0 256 0 256],size(dataR),'spacing',[6,6]); title('Transformed Grid'); drawnow;
    axis([0 256 0 256]);
    if save_images==1
        saveas(gcf,strcat(save_filename,'-images.png'));
    end
end