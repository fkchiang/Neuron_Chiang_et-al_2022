%% ridge regression for two-dimensional decoding
% This script is created mainly for 2D decoding.
% First at all, we would like to get the optimal k for ridge regression.
% Chiang,F-K, et al, NEURON 2021
% email: fengk-kuei.chiang@mssm.edu
clc;clear;close all;
load('EpochArray_zFR_XY_sac.mat','ArrayXY');
% ArrayXY.(SubName).STGON{1,sessions}{1,sig}(trial,[FR,X,Y,sac#])
%% subject: R(1) or Q(2)
for Subject = 1:2 
if Subject == 1     
    Dailysessions = 1:15;
    SubName = 'R';
elseif Subject == 2 
    Dailysessions = 1:10;
    SubName = 'Q';
end
EpochName = {'earlyFix','lateFix','Selection','HOLD','Reward'};
Epoch = {'earlyFixON','lateFixON','WinON','STGON','RewON'};
for sessions = Dailysessions
    for Ep =  1:5 % 5 for hold epoch
        dataEpoch = Epoch{Ep};
        FR = cellfun(@(x) x(:,1),ArrayXY.(SubName).(dataEpoch){1,sessions},'uniformoutput',false);
        FR_array = zscore(cell2mat(FR),0,1);yrange = 100;
        posX = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,2);
        posY = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,3);
        SacNO = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,4);
        num_obs = length(posX);
        num_k_fold = 10;
        %% Compute mean-squared error for regression using 10-fold cross validation
        % mean of square errors
        % (y - yhat)            % Errors
        % (y - yhat).^2         % Squared Error
        % mean((y - yhat).^2)   % Mean Squared Error
        X_set = FR_array;
        y1_set = posX;
        y2_set = posY;
        p = -8:1:8;
        Loc_cvMse_min_posX = zeros(100,num_k_fold);
        Loc_cvMse_min_posY = zeros(100,num_k_fold);
        for iter = 1:10
            %%%% for posX
            [TEST_set,TRAINING_set] = cvParti(num_obs,num_k_fold); % 10-fold cv 
            cvMse_posX = nan(length(p),num_k_fold);
            tic;
            for i = 1:length(p)
                k = 10^p(i);
                for kfold = 1:num_k_fold
                    xTRAIN = X_set(TRAINING_set{1,kfold},:);
                    xTEST  = X_set(TEST_set{1,kfold},:);
                    yTRAIN = y1_set(TRAINING_set{1,kfold},:);
                    y1real = y1_set(TEST_set{1,kfold},1);
                    betas  = ridge(yTRAIN,xTRAIN,k,0);
                    y1hat  = xTEST*betas(2:end)+betas(1);clear betas;
                    cvMse_posX(i,kfold) = mean(abs((y1real - y1hat)).^2);
                end
            end
            for kfold = 1:num_k_fold
                Loc_min = find(cvMse_posX(:,kfold) == min(cvMse_posX(:,kfold)));
                if ismember(1,Loc_min) || ismember(length(p),Loc_min)
                    Loc_cvMse_min_posX(iter,kfold) = nan;        % non-convergence
                else
                    Loc_cvMse_min_posX(iter,kfold) = Loc_min(1); % convergence
                    hpreX = plot(cvMse_posX(:,kfold),'-ob','MarkerFaceColor','b','MarkerEdgeColor','k');hold on;
                end
            end
            toc;
            %%%% for posY
            [TEST_set,TRAINING_set] = cvParti(num_obs,num_k_fold); % 10-fold cv
            cvMse_posY = nan(length(p),num_k_fold);
            tic;
            for i = 1:length(p)
                k = 10^p(i);
                for kfold = 1:num_k_fold
                    xTRAIN = X_set(TRAINING_set{1,kfold},:);
                    xTEST  = X_set(TEST_set{1,kfold},:);
                    yTRAIN = y2_set(TRAINING_set{1,kfold},:);
                    y2real = y2_set(TEST_set{1,kfold},1);
                    betas = ridge(yTRAIN,xTRAIN,k,0);
                    y2hat = xTEST*betas(2:end)+betas(1);clear betas;
                    cvMse_posY(i,kfold) = mean(abs((y2real - y2hat)).^2);
                end
            end
            for kfold = 1:num_k_fold
                Loc_min = find(cvMse_posY(:,kfold) == min(cvMse_posY(:,kfold)));
                if ismember(1,Loc_min) || ismember(length(p),Loc_min)
                    Loc_cvMse_min_posY(iter,kfold) = nan;        % non-convergence
                else
                    Loc_cvMse_min_posY(iter,kfold) = Loc_min(1); % convergence
                    hpreY = plot(cvMse_posY(:,kfold),'-or','MarkerFaceColor','r','MarkerEdgeColor','k');hold on;
                end
            end
        % plot section for Supplemental Figure S5A
        figure(iter);
        [~,nx] = find(cvMse_posX == min(cvMse_posX(:)));
        [~,ny] = find(cvMse_posY == min(cvMse_posY(:)));
        hpreX = plot(cvMse_posX(:,nx),'-b');hold on;
        hpreY = plot(cvMse_posY(:,ny),'-r');hold on;
        set(gca,'tickdir','out','xtick',1:1:17,'xticklabel',{'-8','-7','-6','-5','-4','-3','-2','-1','0',...
                                                              '1','2','3','4','5','6','7','8'});
        title(sprintf('Sbj%s-Ses%.2d-Epoch:%s',(SubName),sessions,(EpochName{Ep})));
        legend([hpreX,hpreY],'pred-X','pred-Y');
        xlabel('ridge k parameter (10^i)');
        ylabel('Mean squared error (MSE)');
        axis square;axis([0 17 20 70]);box off;        
        end
        [xN,xEDGES] = histcounts(Loc_cvMse_min_posX,'BinMethod','integers');
        Loc_optP_posX = [xN; p(xEDGES(1:end-1)+0.5)];
        
        [yN,yEDGES] = histcounts(Loc_cvMse_min_posY,'BinMethod','integers');
        Loc_optP_posY = [yN; p(yEDGES(1:end-1)+0.5)];
                
        optK_sum.(SubName).(dataEpoch).posX{1,sessions} = Loc_optP_posX;clear Loc_optP_posX;
        optK_sum.(SubName).(dataEpoch).posY{1,sessions} = Loc_optP_posY;clear Loc_optP_posY;
    end
end
end
% save('Broad_2D_decoder_intercept_optK.mat','Loc_cvMse_min_posX','Loc_cvMse_min_posY','optK_sum');