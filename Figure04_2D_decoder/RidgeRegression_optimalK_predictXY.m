%% ridge regression for two-dimensional decoding
% This script is created mainly for 2D decoding.
% Here, we do ridge regression with optimal K to predict target X or Y coordinates.
% Chiang,F-K, et al, NEURON 2021
% email: fengk-kuei.chiang@mssm.edu
clc;clear;close all;
load('EpochArray_zFR_XY_sac.mat','ArrayXY');
% ArrayXY.(SubName).STGON{1,sessions}{1,sig}(trial,[FR,X,Y,sac#])
load('Broad_2D_decoder_intercept_optK.mat','optK_sum');
% optK_sum.(SubName).STGON.posX{1,sessions}([freq,10^p])
%% subject: R(1) or Q(2)
for Subject = 1:2
if Subject == 1  
    Dailysessions = 15;
    SubName = 'R';
elseif Subject == 2
    Dailysessions = 10;
    SubName = 'Q';
end
EpochName = {'earlyFix','lateFix','Selection','HOLD','Reward'};
Epoch = {'earlyFixON','lateFixON','WinON','STGON','RewON'};
for sessions = 1:Dailysessions
    for Ep = 1:5     
        dataEpoch = Epoch{Ep};
        FR = cellfun(@(x) x(:,1),ArrayXY.(SubName).(dataEpoch){1,sessions},'uniformoutput',false);
        num_signal = length(FR);

        FR_array = zscore(cell2mat(FR),0,1);
        posX = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,2);
        posY = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,3);
        SacNO = ArrayXY.(SubName).(dataEpoch){1,sessions}{1,1}(:,4);
        num_obs = length(posX);
        num_k_fold = 10;
        %% Compute mean-squared error for regression using 10-fold cross validation
        % mean of square errors
        % (y - yhat)    % Errors
        % (y - yhat).^2   % Squared Error
        % mean((y - yhat).^2)   % Mean Squared Error
        X_set = FR_array;        
        y1_set = posX;
        y2_set = posY;

        pX_range = optK_sum.(SubName).(dataEpoch).posX{1,sessions};
        pY_range = optK_sum.(SubName).(dataEpoch).posY{1,sessions};
        [~,idx_maxX] = max(pX_range(1,:),[],2);
        px = pX_range(2,idx_maxX(1));
        [~,idx_maxY] = max(pY_range(1,:),[],2);
        py = pY_range(2,idx_maxY(1));        
        
        pforx = linspace(10^(px-2),10^(px+2));
        pfory = linspace(10^(py-2),10^(py+2));
        
        iterRun = 10;
        Loc_cvMse_min_posX = zeros(iterRun,num_k_fold);
        Loc_cvMse_min_posY = zeros(iterRun,num_k_fold);

        fprintf('Sbj%s-ses%.2d-Ep%.2d\n',SubName,sessions,Ep);
        for iter = 1:iterRun
        fprintf('iter-%.2d\n',iter);
        [TEST_set,TRAINING_set] = cvParti(num_obs,num_k_fold); % 10-fold cv
        cvMse_posX = nan(length(pforx),num_k_fold);
        % range from 10^-2 to 10^4 based on the result of SFWM_012
            tic;
            for i = 1:length(pforx)
                k = pforx(i);
                for kfold = 1:num_k_fold
                    xTRAIN = X_set(TRAINING_set{1,kfold},:);
                    xTEST  = X_set(TEST_set{1,kfold},:);
                    yTRAIN = y1_set(TRAINING_set{1,kfold},:);
                    y1real = y1_set(TEST_set{1,kfold},1);
                    betas = ridge(yTRAIN,xTRAIN,k,0);
                    y1hat = xTEST*betas(2:end)+betas(1);clear betas;
                    cvMse_posX(i,kfold) = mean(abs((y1real - y1hat)).^2);
                end
            end
            for kfold = 1:num_k_fold
                Loc_min = find(cvMse_posX(:,kfold) == min(cvMse_posX(:,kfold)));
                if ismember(1,Loc_min) || ismember(length(pforx),Loc_min) 
                    Loc_cvMse_min_posX(iter,kfold) = nan;        % non-convergence
                else
                    Loc_cvMse_min_posX(iter,kfold) = Loc_min(1); % convergence
                end
            end
            toc;

            [TEST_set,TRAINING_set] = cvParti(num_obs,num_k_fold); % 10-fold cv
            cvMse_posY = nan(length(pfory),num_k_fold);
            tic;
            for i = 1:length(pfory)
                k = pfory(i);
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
                if ismember(1,Loc_min) || ismember(length(pfory),Loc_min)
                    Loc_cvMse_min_posY(iter,kfold) = nan;        % non-convergence
                else
                    Loc_cvMse_min_posY(iter,kfold) = Loc_min(1); % convergence
                end
            end
            toc;
        end
        [xN,xEDGES] = histcounts(Loc_cvMse_min_posX,'BinMethod','integers');
        Loc_optPx = [xN;xEDGES(1:end-1)+0.5];
        [yN,yEDGES] = histcounts(Loc_cvMse_min_posY,'BinMethod','integers');
        Loc_optPy = [yN;yEDGES(1:end-1)+0.5];
        
        [~,idx_maxX] = max(Loc_optPx(1,:),[],2);
        opt_fineKx = pforx(Loc_optPx(2,idx_maxX(1)));
        [~,idx_maxY] = max(Loc_optPy(1,:),[],2);
        opt_fineKy = pfory(Loc_optPy(2,idx_maxY(1)));
        
        kx = opt_fineKx;
        ky = opt_fineKy;
        num_k_fold = 10;
        [loc_xTEST,loc_xTRAIN] = cvParti(num_obs,num_k_fold);
        [loc_yTEST,loc_yTRAIN] = cvParti(num_obs,num_k_fold);
        
        PredX = nan(ceil(num_obs/num_k_fold),num_k_fold);
        RealX = nan(ceil(num_obs/num_k_fold),num_k_fold);
        
        PredY = nan(ceil(num_obs/num_k_fold),num_k_fold);
        RealY = nan(ceil(num_obs/num_k_fold),num_k_fold);
        normX_dif = nan(1,num_k_fold);
        normY_dif = nan(1,num_k_fold);
        xBeta = nan(num_signal+1,num_k_fold); % include B0 (intercept)
        yBeta = nan(num_signal+1,num_k_fold); % include B0 (intercept)
        for i = 1:num_k_fold
            yTRAIN = y1_set(loc_xTRAIN{1,i});
            xTRAIN = X_set(loc_xTRAIN{1,i},:);
            xBeta(:,i) = ridge(yTRAIN,xTRAIN,kx,0);
            PredX(1:length(loc_xTEST{1,i}),i) = X_set(loc_xTEST{1,i},:)*xBeta(2:end,i) + xBeta(1,i);
            RealX(1:length(loc_xTEST{1,i}),i) = y1_set(loc_xTEST{1,i},1);
            normX_dif(1,i) = diff([norm(PredX(:,i)),norm(RealX(:,i))]); % diff(vector(2) - vector(1));

            yTRAIN = y2_set(loc_yTRAIN{1,i});
            xTRAIN = X_set(loc_yTRAIN{1,i},:);
            yBeta(:,i) = ridge(yTRAIN,xTRAIN,ky,0);
            PredY(1:length(loc_yTEST{1,i}),i) = X_set(loc_yTEST{1,i},:)*yBeta(2:end,i) + yBeta(1,i);
            RealY(1:length(loc_yTEST{1,i}),i) = y2_set(loc_yTEST{1,i},1);
            normY_dif(1,i) = diff([norm(PredY(:,i)),norm(RealY(:,i))]); % diff(vector(2) - vector(1));
        end
        optx_Beta = xBeta(:,normX_dif == min(normX_dif));
        opty_Beta = yBeta(:,normY_dif == min(normY_dif));

        snrxx = X_set*optx_Beta(2:end) + optx_Beta(1);
        prexx = awgn(snrxx,5,'measured');
        snryy = X_set*opty_Beta(2:end) + opty_Beta(1);
        preyy = awgn(snryy,5,'measured');
        
        optXY_betas.(SubName){sessions,Ep} = [optx_Beta,opty_Beta];
        predictionXY.(SubName){sessions,Ep} = [prexx,preyy];
        predictionXY_noAWGN.(SubName){sessions,Ep} = [snrxx,snryy];
    end
end
end
% save('EpochArray_zFR_XY_sac.mat','predictionXY','-append');