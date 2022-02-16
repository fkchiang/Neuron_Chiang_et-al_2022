%% Single unit decoder and best ensemble (BE) procedure
% This script is created mainly for target color decoder with BE procedure.
% Chiang,F-K, et al, NEURON 2021
% email: fengk-kuei.chiang@mssm.edu
clear;clc;close all;
figure('paperorientation','landscape');
load('fwmData.mat','Eyedata');
load('fwmSigSpikeCoInCo500TGON.mat', 'fwmSigSpikeCO500TGON');
%% subject: R(1) or Q(2)
for Subject = 1:2      
    if Subject == 1
        Dailysessions = 1:15;
        SubName = 'R';
    elseif Subject == 2
        Dailysessions = 1:10;
        SubName = 'Q';
    end
    EpochName = {'earlyFixON','lateFixON','WinON','STGON','RewON'};    
    for sessions = Dailysessions
        VarName = sprintf('%s%.3d',SubName,sessions);
        sumTable = Eyedata.(VarName).sumTable;
        NumTrials = length(sumTable);
        BlockNUM = sumTable(:,3);
        idx_CorTrials = ~isnan(sumTable(:,10)); % includes completed
        cor_blockNum = BlockNUM(idx_CorTrials);
        idx_cor_blockNum = unique(cor_blockNum,'stable');
        if Subject == 2 && sessions == 2
            blk_num = 1:5;
        else
            blk_num = 1:6;
        end
%% start with this secion
        for ep = 1:5
            FR = fwmSigSpikeCO500TGON.(SubName).meanFR_Signal.(EpochName{epp(ep)}){sessions};
            for blk = blk_num   % There are six different 6-taget configuration used in a session, so we calculate each one at a time.
                idx_corblock = cor_blockNum == idx_cor_blockNum(blk);
%% prepare FR and TG numbers
                sBlock = cellfun(@(x) x(:,idx_corblock),FR,'uniformoutput',false); % subblocks
                feaMat = cell2mat(cellfun(@(x) x(:),sBlock,'uniformoutput',false));% feature matrix
                idx_tgNum = sumTable(idx_CorTrials,2)';
                idx_rem_tgNum = rem(idx_tgNum+2,3)+1;
                tgNum = repmat(idx_rem_tgNum(idx_corblock),6,1);
                TGvector = tgNum(:);         % note: 6 correct saccades in a trial.
%% decoding section
                [nSacs,nSig] = size(feaMat); % total numbers of saccade in all trials. 
                trainingall = feaMat;        % all of your input data, in the form of [ntrials x features]    
                c = TGvector;                % categories (or labels) for data in trainingall
%%%% Part ONE: singel neuron decoder
                classes    = nan(nSacs,nSig,'single');
                err        = nan(nSacs,nSig,'single');
                posteriors = cell(nSig,1);
                logp       = nan(nSacs,nSig,'single');
                coef       = cell(nSig,1);
                decoMtx    = cell(nSig,1);
                tic;
                parfor j = 1:nSig
                    [sub_classes,sub_err,sub_posteriors,sub_logp,sub_coef] = parfor_Classfy_LOO(nSacs,trainingall,j,c); % Referred to the note section 

                    sub_decoMtx = zeros(3,3,'single');
                    for i = 1:nSacs
                        sub_decoMtx(c(i),sub_classes(i,1)) = sub_decoMtx(c(i),sub_classes(i,1)) + 1;
                    end
                    classes(:,j)     = sub_classes;    % clear sub_classes
                    err(:,j)         = sub_err;        % clear sub_err
                    posteriors{j,1}  = sub_posteriors; % clear sub_posteriors
                    logp(:,j)        = sub_logp;       % clear sub_logp
                    coef{j,1}        = sub_coef;       % clear sub_coef
                    decoMtx{j,1}     = sub_decoMtx;    % clear sub_decoMtx
                end
                ACCU_SIG = cellfun(@(x) sum(diag(x)) / sum(x(:)),decoMtx,'un',true);
                
                toc;                
%%%% Part TWO: Best Ensemble
                [~,Loc_informativeAVE] = sort(ACCU_SIG,'descend');               
                Loc_initial = Loc_informativeAVE(1); % initial index
                preBestSub = [];
                preRestSub = Loc_informativeAVE';
                [preBestSub, preRestSub] = func_idxBstEsb(Loc_initial, preBestSub, preRestSub);
                
                ACCU_AVE        = nan((nSig-1), (nSig-1));
                Loc_maxACCU_AVE = nan(       1, (nSig-1));
                
                while ~isempty(preRestSub)
                    trainingall = feaMat;
                    subSize = length(preBestSub);
                    nPairs  = length(preRestSub);
                    Loc_nPair      = cell(nPairs,1);
                    
                    classes_AVE    = nan(nSacs,nPairs,'single');
                    err_AVE        = nan(nSacs,nPairs,'single');
                    posteriors_AVE = cell(nPairs,1);
                    logp_AVE       = nan(nSacs,nPairs,'single');
                    coef_AVE       = cell(nPairs,1);
                    decoMtx_AVE    = cell(nPairs,1);                    
                    tic;
                    parfor j = 1:nPairs
                        [sub_classes,sub_err,sub_posteriors,sub_logp,sub_coef] = ...
                            parfor_Classfy_LOO(nSacs,trainingall,[preBestSub,preRestSub(j)],c); % Referred to the note section 

                        sub_decoMtx = zeros(3,3,'single');
                        for i = 1:nSacs
                            sub_decoMtx(c(i),sub_classes(i,1)) = sub_decoMtx(c(i),sub_classes(i,1)) + 1;
                        end
                        classes_AVE(:,j)     = sub_classes;    % clear sub_classes
                        err_AVE(:,j)         = sub_err;        % clear sub_err
                        posteriors_AVE{j,1}  = sub_posteriors; % clear sub_posteriors
                        logp_AVE(:,j)        = sub_logp;       % clear sub_logp
                        coef_AVE{j,1}        = sub_coef;       % clear sub_coef
                        decoMtx_AVE{j,1}     = sub_decoMtx;    % clear sub_decoMtx
                        Loc_nPair{j,1}       = [preBestSub,preRestSub(j)];
                    end
                    toc;
                    ACCU_AVE(1:nPairs,subSize) = cellfun(@(x) sum(diag(x)) / sum(x(:)),decoMtx_AVE,'un',true);
                    Loc_maxACCU_AVE(1,subSize) = find(ACCU_AVE(:,subSize) == max(ACCU_AVE(:,subSize)), 1, 'first');                                                          
                    [preBestSub, preRestSub] = func_idxBstEsb(preRestSub(Loc_maxACCU_AVE(1,subSize)), preBestSub, preRestSub);                    
                end
                
                ACCU_BE_AVE = nan(1,nSig-1);
                for ii = 1:(nSig-1)
                    ACCU_BE_AVE(1,ii) = ACCU_AVE(Loc_maxACCU_AVE(ii),ii);
                end                
%%%% plot section
                subplot(2,3,blk);
                hdesac = plot(1:nSig,ACCU_SIG(Loc_informativeAVE),'-k',(1:(nSig-1))+1,ACCU_BE_AVE,'-ok');
                if blk == 1
                    legend(hdesac,{'SU','BE'},'location','northwest');
                    title(sprintf('deColor-BE-%s%.3d-%s',SubName,sessions,EpochName{epp(ep)}));
                    ylabel('Accuracy(%)');
                    xlabel('ensemble size');
                end
                set(gca,'tickdir','out');
                hold on;axis square;axis([0.5 nSig+0.5 0 1.0]);box off;
                  
                decodeColor.(VarName).BE.(EpochName{epp(ep)}).classes{1,blk}  = classes;
                decodeColor.(VarName).BE.(EpochName{epp(ep)}).decoMtx{1,blk}  = decoMtx;
                decodeColor.(VarName).BE.(EpochName{epp(ep)}).ACCU_SIG{1,blk} = ACCU_SIG;

                decodeColor.(VarName).cumuAVE.(EpochName{epp(ep)}).Loc_informativeAVE{1,blk} = Loc_informativeAVE;
                decodeColor.(VarName).cumuAVE.(EpochName{epp(ep)}).classes_AVE{1,blk}        = classes_AVE;
                decodeColor.(VarName).cumuAVE.(EpochName{epp(ep)}).decoMtx_AVE{1,blk}        = decoMtx_AVE;
                decodeColor.(VarName).cumuAVE.(EpochName{epp(ep)}).ACCU_AVE{1,blk}           = ACCU_AVE;
                decodeColor.(VarName).cumuAVE.(EpochName{epp(ep)}).Loc_maxACCU_AVE{1,blk}    = Loc_maxACCU_AVE;
            end
            print(gcf,'-dpdf','-r300','-fillpage',sprintf('decodeColor-BE-%s%.3d-%s',SubName,sessions,EpochName{epp(ep)}));
            pause(2);
            clf;
        end
    end
end
% save('decodeColor_BE.mat','decodeColor');