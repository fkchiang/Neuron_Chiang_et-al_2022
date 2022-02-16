%% Single unit decoder and best single unit (BSU) procedure
% This script is created mainly for saccade order decoder with BSU procedure.
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
        if Subject == 1 && sessions == 11
            sumTable = sumTable(1:222,:);
        end
        NumTrials = length(sumTable);
        BlockNUM = sumTable(:,3);
        idx_CorTrials = ~isnan(sumTable(:,10));             % includes completed
        cor_blockNum = BlockNUM(idx_CorTrials);
        idx_cor_blockNum = unique(cor_blockNum,'stable');
        if Subject == 2 && sessions == 2
            blk_num = 1:5;
        else
            blk_num = 1:6;
        end
%% start with this secion
        for ep = 1:6
            FR = fwmSigSpikeCO500TGON.(SubName).meanFR_Signal.(EpochName{ep}){sessions};
            for blk = blk_num   % There are six different 6-taget configuration used in a session, so we calculate each one at a time.
                idx_corblock = cor_blockNum == idx_cor_blockNum(blk);
%% prepare FR and TG numbers
                sBlock = cellfun(@(x) x(:,idx_corblock),FR,'uniformoutput',false); % subblocks
                feaMat = cell2mat(cellfun(@(x) x(:),sBlock,'uniformoutput',false));% feature matrix
                idx_tgNum = sumTable(idx_CorTrials,5:10)';
                trials_Blk = sum(idx_corblock);
                tgNum = repmat((1:6)',1,trials_Blk);
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

                    sub_decoMtx = zeros(6,6,'single');
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
                toc;
%%%% Part TWO: best single ensemble
                [~,idx_informativeAVE] = sort(cellfun(@(x) sum(diag(x)) / sum(x(:)),decoMtx,'un',true),'descend');

                trainingAVE = feaMat(:,idx_informativeAVE);
                classes_AVE    = nan(nSacs,nSig,'single');
                err_AVE        = nan(nSacs,nSig,'single');
                posteriors_AVE = cell(nSig,1);
                logp_AVE       = nan(nSacs,nSig,'single');
                coef_AVE       = cell(nSig,1);
                decoMtx_AVE    = cell(nSig,1);
                tic;
                parfor j = 1:nSig      
                    [sub_classes,sub_err,sub_posteriors,sub_logp,sub_coef] = parfor_Classfy_LOO(nSacs,trainingAVE,1:1:j,c); % Referred to the note section 

                    sub_decoMtx = zeros(6,6,'single');
                    for i = 1:nSacs
                        sub_decoMtx(c(i),sub_classes(i,1)) = sub_decoMtx(c(i),sub_classes(i,1)) + 1;
                    end
                    classes_AVE(:,j)     = sub_classes;    % clear sub_classes
                    err_AVE(:,j)         = sub_err;        % clear sub_err
                    posteriors_AVE{j,1}  = sub_posteriors; % clear sub_posteriors
                    logp_AVE(:,j)        = sub_logp;       % clear sub_logp
                    coef_AVE{j,1}        = sub_coef;       % clear sub_coef
                    decoMtx_AVE{j,1}     = sub_decoMtx;    % clear sub_decoMtx
                end
                toc;
%%%% plot section  
                ACCU_SIG = cellfun(@(x) sum(diag(x)) / sum(x(:)),decoMtx,'un',true);
                ACCU_AVE = cellfun(@(x) sum(diag(x)) / sum(x(:)),decoMtx_AVE,'un',true);            
                subplot(2,3,blk);
                hdesac=plot(1:nSig,ACCU_SIG(idx_informativeAVE),'-k',1:nSig,ACCU_AVE,'-ok');
                if blk == 1
                    legend(hdesac,{'BSU','cumuBSU'},'location','northwest');
                    title(sprintf('deSac-BSU-%s%.3d-%s',SubName,sessions,EpochName{ep}));
                    ylabel('Accuracy(%)');
                    xlabel('ensemble size');
                end
                set(gca,'tickdir','out');
                hold on;axis square;axis([0.5 nSig+0.5 0 1.0]);box off;
                
                decodeSac.(VarName).BSU.(EpochName{ep}).classes{1,blk}  = classes;
                decodeSac.(VarName).BSU.(EpochName{ep}).decoMtx{1,blk}  = decoMtx;
                decodeSac.(VarName).BSU.(EpochName{ep}).ACCU_SIG{1,blk} = ACCU_SIG;

                decodeSac.(VarName).cumuAVE.(EpochName{ep}).idx_informativeAVE{1,blk} = idx_informativeAVE;
                decodeSac.(VarName).cumuAVE.(EpochName{ep}).classes_AVE{1,blk}        = classes_AVE;
                decodeSac.(VarName).cumuAVE.(EpochName{ep}).decoMtx_AVE{1,blk}        = decoMtx_AVE;
                decodeSac.(VarName).cumuAVE.(EpochName{ep}).ACCU_AVE{1,blk}           = ACCU_AVE;
            end
            print(gcf,'-dpdf','-r300','-fillpage',sprintf('decodeSac-BSU-%s%.3d-%s',SubName,sessions,EpochName{ep}));
            pause(2);
            clf;
        end
    end
end
% save('decodeSac_BSU.mat','decodeSac');