% Compare TPS vs LCM using hands data and automatically detected landmarks
% with variations on number of missing landmarks

clear; clc; close  all;
SPIEJMI_hand_detectPOIs

adderror = 0; num_errors = 0;
save_file = 0; save_images = 0;

LMalpha_arr = [1]; 
CMalpha_arr = reshape(logspace(0,7,8),[2,4])';

% specify desired number of CM points in the first arc per experiment (each
% row should correspond to the number of desired number of CM points in the
% first arc)
first_arcRxi = [2;3;5;10;15];
Rxi_num = SPIEJMI_hand_getnumberofCMpoints(R_ordered_POIs,first_arcRxi);

% initialize sample points for dataR and dataT contours
Rxi = []; Txi = [];

Dice_arr = zeros(size(Rxi_num,1),6);
aveTREmissingLMs = Dice_arr;

if (adderror==1 & num_errors==1)
    % keep track of TREs for each possible landmark error + runtime
    TRE_TPS_arr = zeros(POI_no,length(LMalpha_arr),size(CMalpha_arr,1),size(Rxi_num,1));
    TRE_TPSall_arr = TRE_TPS_arr;
    TRE_Small_arr =  TRE_TPS_arr;
    TRE_Large_arr = TRE_TPS_arr;
    tictoc_TPS = zeros(POI_no,length(LMalpha_arr),size(CMalpha_arr,1),size(Rxi_num,1));   
else
    % keep track of runtime only
    tictoc_TPS = zeros(1,length(LMalpha_arr),size(CMalpha_arr,1),size(Rxi_num,1));   
end

tictoc_TPSall = tictoc_TPS; tictoc_Small = tictoc_TPS; tictoc_Large = tictoc_TPS;
 
if (adderror==1 & num_errors==2) | (adderror==0)
    % no need to run through all possible t_i errors since there is
    % only one possible location of incorrect landmarks (t1 and t5 if
    % numerrors is 2, and only one setup to consider if no errors)
    num_comb = 1;
elseif (adderror==1 & num_errors==1)
    num_comb = POI_no;
end

for POI_loc = 1:num_comb    
    for LMalpharow = 1:length(LMalpha_arr)
        LMalpha = LMalpha_arr(LMalpharow);
        
        for CMalpharow = 1:size(CMalpha_arr,1)
            LMvariation = 1;
            
            LMindex = 1:POI_no;
            
            for MajorLMnum = 1:size(Rxi_num,1)
                dataR = double(imread('hands-R.jpg'));
                dataT = double(imread('hands-T.jpg'));
                dataR = imresize(dataR,4);
                dataT = imresize(dataT,4);
                
                % initialize sample points for dataR and dataT contours
                Rxi = []; Txi = [];
                h1 = figure('Name','Binarized dataR'); set(h1, 'Visible', 'off');
                imshow(dataR,[]); hold on;
                for i = 1:POI_no
                    if i+1<=POI_no
                        Rcurr_POI = [R_ordered_POIs(i,:);R_ordered_POIs(i+1,:)];
                        Tcurr_POI = [T_ordered_POIs(i,:);T_ordered_POIs(i+1,:)];
                    else
                        Rcurr_POI = [R_ordered_POIs(i,:);R_ordered_POIs(mod(i+1,POI_no),:)];
                        Tcurr_POI = [T_ordered_POIs(i,:);T_ordered_POIs(mod(i+1,POI_no),:)];
                    end
                    
                    Rcurr_POI_contour_ind1 = find(R_edgecoords(:,1)==Rcurr_POI(1,1) & R_edgecoords(:,2)==Rcurr_POI(1,2));
                    Rcurr_subint = circshift(R_edgecoords,[-Rcurr_POI_contour_ind1+1 0]);
                    
                    Rcurr_POI_contour_ind1 = 1; % it's now at the start of the list
                    Rcurr_POI_contour_ind2 = find(Rcurr_subint(:,1)==Rcurr_POI(2,1) & Rcurr_subint(:,2)==Rcurr_POI(2,2));
                    plot(Rcurr_subint(Rcurr_POI_contour_ind1:Rcurr_POI_contour_ind2,1),Rcurr_subint(Rcurr_POI_contour_ind1:Rcurr_POI_contour_ind2,2),'g.','MarkerSize',10);
                    
                    Tcurr_POI_contour_ind1 = find(T_edgecoords(:,1)==Tcurr_POI(1,1) & T_edgecoords(:,2)==Tcurr_POI(1,2));
                    Tcurr_subint = circshift(T_edgecoords,[-Tcurr_POI_contour_ind1+1 0]);
                    Tcurr_POI_contour_ind1 = 1; % it's now at the start of the list
                    Tcurr_POI_contour_ind2 = find(Tcurr_subint(:,1)==Tcurr_POI(2,1) & Tcurr_subint(:,2)==Tcurr_POI(2,2));
                    
                    subint_minlength = min([Rcurr_POI_contour_ind2,Tcurr_POI_contour_ind2]);
                    subint_n = Rxi_num(MajorLMnum,i);
                    
                    if subint_n>subint_minlength
                        subint_n = subint_minlength;
                    end
                    
                    Rcurr_arc = Rcurr_subint(1:Rcurr_POI_contour_ind2,:);
                    Tcurr_arc = Tcurr_subint(1:Tcurr_POI_contour_ind2,:);
                    Raddthis = SPIEJMI_getXiStars(Rcurr_arc,subint_n);
                    Taddthis = SPIEJMI_getXiStars(Tcurr_arc,subint_n);
                    Rxi = [Rxi;Raddthis];
                    Txi = [Txi;Taddthis];
                    
                    plot(Raddthis(:,1),Raddthis(:,2),'b.','MarkerSize',20);
                    % pause;
                end
                
                Rpoints = R_ordered_POIs; Tpoints = T_ordered_POIs;
                Rpoints_TPSall = Rxi; Tpoints_TPSall = Txi;
                
                %% add error to specific Tpoints
                if adderror==1
                    errorvec = zeros(size(Tpoints));
                    if num_errors == 1
                        errorloc = [POI_loc];
                        wrongtloc_str = num2str(POI_loc);
                    elseif num_errors == 2
                        errorloc = [1 5];
                        wrongtloc_str = '1And5';
                    end
                    
                    errorvec(errorloc,:) = ones(length(errorloc),2); % add error to specific major LMs
                    randn('seed',0); r = 0.05*size(dataT,1) + 2.*randn(size(Tpoints)); randn('seed',0);
                    errorvec = errorvec.*r;
                    Tpoints_orig = Tpoints;
                    Tpoints = Tpoints + errorvec;
                    
                    % set up the Tpoints (with error) to be used for TPS + all LMs
                    % identify the indices of the Rpoints/Tpoints where errors are to be added
                    T_adderrortoCMpts_ind = find(ismember(Tpoints_TPSall,Tpoints_orig(errorloc,:),'rows')==1);
                    Tpoints_TPSall(T_adderrortoCMpts_ind,:) = Tpoints(errorloc,:);
                    Txi(T_adderrortoCMpts_ind,:) = [];
                    Rxi(T_adderrortoCMpts_ind,:) = [];
                end
                
                num_LM = size(Rpoints,1); num_xi = size(Rxi,1);
                
                [prereg_dice,prereg_jaccard] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(dataT));
                
                %% Use iterative TPS + Major LMs only (LCMTPS where CMalpha=0)
                CMalpha = 0; TPSalpha = 0;
                titlestr = strcat('TPS: ', num2str(num_LM),'LMs');
                save_filename = strcat('TPS-',num2str(num_LM),'LMs-TPSonly');
                tic;
                [v_TPS,f_LM_TPS,f_CM_TPS,Xq_TPS,Yq_TPS,transfTpoints_TPS,transfTxi_TPS,Ty_TPS] = SPIEJMI_LCM_2d(dataR,dataT,Rpoints,Tpoints,Rxi,Txi,LMalpha,CMalpha,titlestr,save_filename,save_file,save_images);
                TPStoc = toc; 
                [postreg_dice_TPS,postreg_jaccard_TPS] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(Ty_TPS));
                Dice_arr(MajorLMnum,2*(LMvariation-1)+1) = postreg_dice_TPS;
                
                pause(3);
                close all;
                
                %% Use iterative TPS + All Landmarks
                CMalpha = 0; TPSalpha = 0;
                titlestr = strcat('TPS: ', num2str(size(Rxi,1)),'LMs');
                save_filename = strcat('TPS-',num2str(num_LM),'LMs-TPSonly-AllLMs');
                tic;
                [v_TPSall,f_LM_TPSall,f_CM_TPSall,Xq_TPSall,Yq_TPSall,transfTpoints_TPSall,transfTxi_TPSall,Ty_TPSall] = SPIEJMI_LCM_2d(dataR,dataT,Rpoints_TPSall,Tpoints_TPSall,Rxi,Txi,LMalpha,CMalpha,titlestr,save_filename,save_file,save_images);
                TPSalltoc = toc;
                [postreg_dice_TPSall,postreg_jaccard_TPSall] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(Ty_TPSall));
                
                pause(3);
                close all;
                
                %% Use LCMTPS Small CMalpha
                CMalpha = CMalpha_arr(CMalpharow,1); CMsmall = CMalpha; TPSalpha = 0;
                titlestr = strcat(num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs-','CM',num2str(CMalpha),'-TPS',num2str(TPSalpha));
                save_filename = strcat('LCM-',num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs-','CM',num2str(CMalpha),'-TPS',num2str(TPSalpha));
                tic;
                [v_LCMsmall,f_LM_LCMsmall,f_CM_LCMsmall,Xq_LCMsmall,Yq_LCMsmall,transfTpoints_LCMsmall,transfTxi_LCMsmall,Ty_LCMsmall] = SPIEJMI_LCM_2d(dataR,dataT,Rpoints,Tpoints,Rxi,Txi,LMalpha,CMalpha,titlestr,save_filename,save_file,save_images);
                Smalltoc = toc; 
                [postreg_dice_LCMsmall,postreg_jaccard_LCMsmall] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(Ty_LCMsmall));
                
                pause(3);
                close all;
                
                %% Use LCMTPS Large CMalpha
                CMalpha = CMalpha_arr(CMalpharow,2); CMlarge = CMalpha; TPSalpha = 0;
                titlestr = strcat(num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs-','CM',num2str(CMalpha),'-TPS',num2str(TPSalpha));
                save_filename = strcat('LCM-',num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs-','CM',num2str(CMalpha),'-TPS',num2str(TPSalpha));
                tic;
                [v_LCMLarge,f_LM_LCMLarge,f_CM_LCMLarge,Xq_LCMLarge,Yq_LCMLarge,transfTpoints_LCMLarge,transfTxi_LCMLarge,Ty_LCMLarge] = SPIEJMI_LCM_2d(dataR,dataT,Rpoints,Tpoints,Rxi,Txi,LMalpha,CMalpha,titlestr,save_filename,save_file,save_images);
                Largetoc = toc;
                [postreg_dice_LCMLarge,postreg_jaccard_LCMLarge] = SPIEJMI_diceANDjaccardindex_2d(imbinarize(dataR), imbinarize(Ty_LCMLarge));
                Dice_arr(MajorLMnum,2*LMvariation) = postreg_dice_LCMLarge;
                
                pause(3);
                close all;
                
                pause(3);
                close all;
                
                %% Create image comparing TPS vs LCM results
                h = figure('units','normalized','outerposition',[0 0 1 1]); % set(h, 'Visible', 'off');
                
                sgtitle(strcat('LMalpha=',num2str(LMalpha),'LCM vs TPS: ', num2str(size(Rpoints,1)),'LM MajorLMs - ',num2str(size(Rxi,1)),'ApproxLMs - Pre-registration Dice: ',num2str(prereg_dice,2)));
                subplot(4,4,1); imshow(abs(dataR-dataT),[]); title('|R-T|, TPS Major LMs only'); hold on; plot(Rpoints(:,1),Rpoints(:,2),'m*','MarkerSize',10);
                subplot(4,4,5); imshow(abs(dataR-dataT),[]); title(strcat('|R-T|, TPS All LMs')); hold on; plot(Rxi(:,1),Rxi(:,2),'m*','MarkerSize',10);
                subplot(4,4,9); imshow(abs(dataR-dataT),[]); title(strcat('CM',num2str(CMsmall))); hold on; plot(Rpoints(:,1),Rpoints(:,2),'m*','MarkerSize',10); plot(Rxi(:,1),Rxi(:,2),'m.','MarkerSize',5);
                subplot(4,4,13); imshow(abs(dataR-dataT),[]); title(strcat('CM',num2str(CMlarge))); hold on; plot(Rpoints(:,1),Rpoints(:,2),'m*','MarkerSize',10); plot(Rxi(:,1),Rxi(:,2),'m.','MarkerSize',5);
                
                % TPS Using Major LMs Only
                subplot(4,4,2); imshow(Ty_TPS,[]); title(strcat('T[y]_{TPS}, Dice: ',num2str(postreg_dice_TPS,2)));
                subplot(4,4,3); imshow(abs(dataR-Ty_TPS),[]); title('|R - T[y]_{TPS}|');
                subplot(4,4,4); JanModersitzki_plotGrid([Xq_TPS(:);Yq_TPS(:)],[0 size(dataR,1) 0 size(dataR,2)],size(dataR),'spacing',[8,8]); title('TPS Transf. Grid'); drawnow; axis([0 size(dataR,1) 0 size(dataR,2)]);
                
                % TPS Using Major And Approx
                subplot(4,4,6); imshow(Ty_TPSall,[]); title(strcat('T[y]_{TPSall}, Dice: ',num2str(postreg_dice_TPSall,2)));
                subplot(4,4,7); imshow(abs(dataR-Ty_TPSall),[]); title('|R - T[y]_{TPSall}|');
                subplot(4,4,8); JanModersitzki_plotGrid([Xq_TPSall(:);Yq_TPSall(:)],[0 size(dataR,1) 0 size(dataR,2)],size(dataR),'spacing',[8,8]); title('TPSall Transf. Grid'); drawnow; axis([0 size(dataR,1) 0 size(dataR,2)]);
                
                % LCM_SmallCMAlpha
                subplot(4,4,10); imshow(Ty_LCMsmall,[]); title(strcat('T[y]_{LCMSmall}, Dice: ',num2str(postreg_dice_LCMsmall,2)));
                subplot(4,4,11); imshow(abs(dataR-Ty_LCMsmall),[]); title('|R - T[y]_{LCMSmall}|');
                subplot(4,4,12); JanModersitzki_plotGrid([Xq_LCMsmall(:);Yq_LCMsmall(:)],[0 size(dataR,1) 0 size(dataR,2)],size(dataR),'spacing',[8,8]); title(strcat('LCMSmall Transf. Grid - CMtoLMRatio=',num2str(f_CM_LCMsmall/f_LM_LCMsmall,2))); drawnow; axis([0 size(dataR,1) 0 size(dataR,2)]);
                
                % LCM_LargeCMAlpha
                subplot(4,4,14); imshow(Ty_LCMLarge,[]); title(strcat('T[y]_{LCMLarge}, Dice: ',num2str(postreg_dice_LCMLarge,2)));
                subplot(4,4,15); imshow(abs(dataR-Ty_LCMLarge),[]); title('|R - T[y]_{LCMLarge}|');
                subplot(4,4,16); JanModersitzki_plotGrid([Xq_LCMLarge(:);Yq_LCMLarge(:)],[0 size(dataR,1) 0 size(dataR,2)],size(dataR),'spacing',[8,8]); title(strcat('LCMLarge Transf. Grid - CMtoLMRatio=',num2str(f_CM_LCMLarge/f_LM_LCMLarge,2))); drawnow; axis([0 size(dataR,1) 0 size(dataR,2)]);
                
                if adderror==1
                    % compute TREs
                    [transfTpoints_correctTPS(:,1), transfTpoints_correctTPS(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints_orig,Tpoints,v_TPS,size(Tpoints_orig(:,1)));
                    [transfTpoints_correctTPSall(:,1), transfTpoints_correctTPSall(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints_orig,Tpoints_TPSall,v_TPSall,size(Tpoints_orig(:,1)));
                    [transfTpoints_correctLCMsmall(:,1), transfTpoints_correctLCMsmall(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints_orig,Tpoints,v_LCMsmall,size(Tpoints_orig(:,1)));
                    [transfTpoints_correctLCMLarge(:,1), transfTpoints_correctLCMLarge(:,2)] = SPIEJMI_computeTPStransfvec_2d(Tpoints_orig,Tpoints,v_LCMLarge,size(Tpoints_orig(:,1)));
                    
                    if num_errors==1
                        TRE_TPS = sqrt((transfTpoints_correctTPS(POI_loc,1)-Rpoints(POI_loc,1))^2 + (transfTpoints_correctTPS(POI_loc,2)-Rpoints(POI_loc,2))^2);
                        TRE_TPSall = sqrt((transfTpoints_correctTPSall(POI_loc,1)-Rpoints(POI_loc,1))^2 + (transfTpoints_correctTPSall(POI_loc,2)-Rpoints(POI_loc,2))^2);
                        TRE_Small = sqrt((transfTpoints_correctLCMsmall(POI_loc,1)-Rpoints(POI_loc,1))^2 + (transfTpoints_correctLCMsmall(POI_loc,2)-Rpoints(POI_loc,2))^2);
                        TRE_Large = sqrt((transfTpoints_correctLCMLarge(POI_loc,1)-Rpoints(POI_loc,1))^2 + (transfTpoints_correctLCMLarge(POI_loc,2)-Rpoints(POI_loc,2))^2);
                    elseif num_errors==2
                        TRE_TPS = 0.5*((sqrt((transfTpoints_correctTPS(errorloc(1),1)-Rpoints(errorloc(1),1))^2 + (transfTpoints_correctTPS(errorloc(1),2)-Rpoints(errorloc(1),2))^2)) + (sqrt((transfTpoints_correctTPS(errorloc(2),1)-Rpoints(errorloc(2),1))^2 + (transfTpoints_correctTPS(errorloc(2),2)-Rpoints(errorloc(2),2))^2)));
                        TRE_TPSall = 0.5*((sqrt((transfTpoints_correctTPSall(errorloc(1),1)-Rpoints(errorloc(1),1))^2 + (transfTpoints_correctTPSall(errorloc(1),2)-Rpoints(errorloc(1),2))^2)) + (sqrt((transfTpoints_correctTPSall(errorloc(2),1)-Rpoints(errorloc(2),1))^2 + (transfTpoints_correctTPSall(errorloc(2),2)-Rpoints(errorloc(2),2))^2)));
                        TRE_Small = 0.5*((sqrt((transfTpoints_correctLCMsmall(errorloc(1),1)-Rpoints(errorloc(1),1))^2 + (transfTpoints_correctLCMsmall(errorloc(1),2)-Rpoints(errorloc(1),2))^2)) + (sqrt((transfTpoints_correctLCMsmall(errorloc(2),1)-Rpoints(errorloc(2),1))^2 + (transfTpoints_correctLCMsmall(errorloc(2),2)-Rpoints(errorloc(2),2))^2)));
                        TRE_Large = 0.5*((sqrt((transfTpoints_correctLCMLarge(errorloc(1),1)-Rpoints(errorloc(1),1))^2 + (transfTpoints_correctLCMLarge(errorloc(1),2)-Rpoints(errorloc(1),2))^2)) + (sqrt((transfTpoints_correctLCMsmall(errorloc(2),1)-Rpoints(errorloc(2),1))^2 + (transfTpoints_correctLCMsmall(errorloc(2),2)-Rpoints(errorloc(2),2))^2)));
                    end
                end
                
                if num_errors==1
                    TRE_TPS_arr(POI_loc,LMalpharow,CMalpharow,MajorLMnum) = TRE_TPS;
                    TRE_TPSall_arr(POI_loc,LMalpharow,CMalpharow,MajorLMnum) = TRE_TPSall;
                    TRE_Small_arr(POI_loc,LMalpharow,CMalpharow,MajorLMnum) = TRE_Small;
                    TRE_Large_arr(POI_loc,LMalpharow,CMalpharow,MajorLMnum) = TRE_Large;
                    runtime_loc = POI_loc;
                else
                    runtime_loc = 1;
                end
                tictoc_TPS(runtime_loc,LMalpharow,CMalpharow,MajorLMnum) = TPStoc;
                tictoc_TPSall(runtime_loc,LMalpharow,CMalpharow,MajorLMnum) = TPSalltoc;
                tictoc_Small(runtime_loc,LMalpharow,CMalpharow,MajorLMnum) = Smalltoc;
                tictoc_Large(runtime_loc,LMalpharow,CMalpharow,MajorLMnum) = Largetoc;
                
                if adderror==0
                    save_filename = strcat(num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs_LMalpha',num2str(LMalpha),'_CMalphas',num2str(CMsmall),'-',num2str(CMlarge));
                elseif adderror==1 & num_errors==1
                    save_filename = strcat('With',num2str(num_errors),'errors-','-IncorrectT',wrongtloc_str,'-',num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs_LMalpha',num2str(LMalpha),'_CMalphas',num2str(CMsmall),'-',num2str(CMlarge));
                elseif adderror==1 & num_errors==2
                    save_filename = strcat('With',num2str(num_errors),'errors-','-IncorrectT',wrongtloc_str,'-',num2str(num_LM),'LMs-',num2str(num_xi),'ApproxLMs_LMalpha',num2str(LMalpha),'_CMalphas',num2str(CMsmall),'-',num2str(CMlarge));
                end
                
                if save_images
                    saveas(gcf,strcat(save_filename,'.png')); 
                end
                
                if save_file
                    save(strcat(save_filename,'-data.mat'));
                end
                pause(60);
            end
        end
    end
end