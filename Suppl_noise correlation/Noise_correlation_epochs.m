%% Noise correlation
% This script is created mainly for noise correlation at each epoch.
% Chiang,F-K, et al, NEURON 2021
% email: fengk-kuei.chiang@mssm.edu
clear;clc;close all;
figure('paperorientation','landscape');
load('fwmData.mat','Eyedata');
load('fwmSigSpikeCoInCo500TGON.mat', 'fwmSigSpikeCO500TGON');
load('SFWMData.mat', 'celllist');
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
%% cell list for remove r-correlation pair recorded in the same electrode
        idx_ses = celllist.(SubName)(:,2) == sessions;
        Loc_num_ch_sig = celllist.(SubName)(idx_ses,3:5);
        Loc_ch = unique(Loc_num_ch_sig(:,2));
        num_ch = length(Loc_ch);
        elementNAN = [];
        for i = 1:num_ch
            if sum(Loc_num_ch_sig(:,2) == Loc_ch(i)) == 0
                error('check the cell list');
            elseif sum(Loc_num_ch_sig(:,2) == Loc_ch(i)) == 1
                continue
            elseif sum(Loc_num_ch_sig(:,2) == Loc_ch(i)) >= 2
                idx = Loc_num_ch_sig(:,2) == Loc_ch(i);
                list_signum = Loc_num_ch_sig(idx,1);
                elementNAN = [elementNAN; nchoosek(list_signum(1):list_signum(end),2)];
            end
        end
%% start with this secion
        for ep = 1:6
            FR = fwmSigSpikeCO500TGON.(SubName).meanFR_Signal.(EpochName{ep}){sessions};
            for blk = blk_num   % There are six different 6-taget configuration used in a session, so we calculate each one at a time.
                idx_corblock = cor_blockNum == idx_cor_blockNum(blk);
%% prepare FR and TG numbers
                sBlock = cellfun(@(x) x(:,idx_corblock),FR,'uniformoutput',false);  % subblocks
                feaMat = cell2mat(cellfun(@(x) x(:),sBlock,'uniformoutput',false)); % feature matrix
                idx_tgNum = sumTable(idx_CorTrials,5:10)';
                tgNum = idx_tgNum(:,idx_corblock);
                TGvector = tgNum(:); % note: 6 correct saccades in a trial.
%% coefficient of noise correlarion
                [coefMat_R,coefMat_pval] = corrcoef(feaMat);
                if any(elementNAN)
                    m = size(elementNAN,1);
                    for i = 1:m
                        coefMat_R(elementNAN(i,1),elementNAN(i,2)) = nan;
                    end
                end
                upT = triu(coefMat_R,1);  % get upper triangle part
                upT_nonzero = upT(upT~=0);
                upT_nonzero_nonnan = upT_nonzero(~isnan(upT_nonzero));
                upT_nonzero_nonnan_ztransformed = 0.5 .* log((1+upT_nonzero_nonnan)./(1-upT_nonzero_nonnan));

                upt_pval = triu(coefMat_pval,1);
                upt_pval_nonzero = upt_pval(upt_pval~=0);
                upT_pval_nonzero_nonnan = upt_pval_nonzero(~isnan(upT_nonzero));
                upT_pval_nonzero_nonnan_idx_sig = upT_pval_nonzero_nonnan < 0.05;
                
                sigTotalPosPairs = sum(upT_pval_nonzero_nonnan_idx_sig & (upT_nonzero_nonnan_ztransformed > 0));
                sigTotalNegPairs = sum(upT_pval_nonzero_nonnan_idx_sig & (upT_nonzero_nonnan_ztransformed < 0));

%% plot: visualization
                histogram(upT_nonzero_nonnan_ztransformed,'BinWidth',0.05,'BinLimits',[-1,1]);
                axis([-1.05 1.05 0 300]);
                axis square;
                pause(2);
                close all;
%% data summary
                numSigs = size(feaMat,2);
                totalPairs = length(upT_nonzero);
                removePairs = size(elementNAN,1);
                sigTotalPosPairs;
                sigTotalNegPairs;
%% file saved
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.elementNAN = elementNAN;
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.coefMat_R = triu(coefMat_R,1);
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.coefMat_pval = triu(coefMat_pval,1);
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.coefList_R_nonzero_nonnan = upT_nonzero_nonnan;
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.coefList_R_nonzero_nonnan_ztransformed = upT_nonzero_nonnan_ztransformed;
                NoiseCor.(SubName).(EpochName{ep}){sessions,1}{1,blk}.summaryNum = [numSigs, totalPairs, removePairs, sigTotalPosPairs, sigTotalNegPairs];
            end
        end
    end
end
% save('Suppl_NoiseCor.mat','NoiseCor');