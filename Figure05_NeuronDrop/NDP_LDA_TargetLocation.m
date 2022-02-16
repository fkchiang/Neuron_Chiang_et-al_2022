%% Liner Discriment Analysis for target location with neuron drop procedure
% This script is created mainly for target location decoding.
% Chiang,F-K, et al, NEURON 2021
% email: fengk-kuei.chiang@mssm.edu
clear;clc;close all;
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
    disp(sessions);
    VarName = sprintf('%s%.3d',SubName,sessions);
    sumTable = Eyedata.(VarName).sumTable;
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
    for blk = blk_num
        loc_config = find(cor_blockNum == idx_cor_blockNum(blk));
        idx_corblock = cor_blockNum == idx_cor_blockNum(blk);
%% prepare FR and TG numbers
        for ep = 1:5
        FR = fwmSigSpikeCO500TGON.(SubName).meanFR_Signal.(EpochName{ep}){sessions};
        sBlock = cellfun(@(x) x(:,idx_corblock),FR,'uniformoutput',false);
        feaMat = cell2mat(cellfun(@(x) x(:),sBlock,'uniformoutput',false));
        idx_tgNum = sumTable(idx_CorTrials,5:10)';            
        tgNum = idx_tgNum(:,idx_corblock);
        TGvector = tgNum(:);
%% decoding section
        ntrials = length(TGvector);
        tic;
        for tdrop = 1:size(feaMat,2) % neuron dropping
            trainingall = feaMat;    % all of your input data, in the form of [ntrials x features]
            trainingall(:,tdrop) = [];
            c = TGvector;            % categories (or labels) for data in trainingall         
            classes    = nan(ntrials,1);
            err        = nan(ntrials,1);
            posteriors = nan(ntrials,6);
            logp       = nan(ntrials,1);
            coef       = cell(ntrials,1);
            for j = 1:ntrials
                training=trainingall;
                touse=std(training,0,1);
                training=training(:,touse~=0); % only use features w non-zero variance
                test=trainingall(j,:);
                test=test(:,touse~=0);
                training(j,:) = [];
                c2=c;
                c2(j)=[];                   
                [classes(j,1),err(j,1),posteriors(j,:),logp(j,1),coef{j,1}]=classify(test,training,c2);
            end
            CovMatrix = zeros(6,6);
            for i = 1:ntrials
                CovMatrix(TGvector(i),classes(i)) = CovMatrix(TGvector(i),classes(i)) + 1;
            end            
            decodeTG_CovMatrix.(VarName).(EpochName{ep}){tdrop,blk} = CovMatrix;
            toc;
        end
        end
    end
end
end
% save('decodeTG_NeuronDrop.mat','decodeTG_CovMatrix');clear decodeTG