%% High/low SI spatial tuning plots
% This script is created mainly for High/low SI spatial tuning plots.
% Chiang,F-K, et al, NEURON 2021
% email: feng-kuei.chiang@mssm.edu
clc;clear;close all;
load('spatialTuning_polar_dataset_reducedmodel.mat','mean_rad_LSI','mean_mag_LSI','mean_rad_HSI','mean_mag_HSI',...
     'idx_HSI_ftpOnlyXSIG','idx_HSI_ftpOnlyYSIG','idx_HSI_ftpXorYSIG','idx_HSI_ftpXandYSIG',...
     'idx_LSI_ftpOnlyXSIG','idx_LSI_ftpOnlyYSIG','idx_LSI_ftpXorYSIG','idx_LSI_ftpXandYSIG',...
     'idx_SI_Full_XorYSIG','sampleVector');
%% sample neuron: R015, signal#29
sdata = scaledata([sampleVector.weight_LSI;sampleVector.weight_HSI],0.01,1);
Loc_LSI = 1:length(sampleVector.weight_LSI);
Loc_HSI = (Loc_LSI(end)+1):length(sdata);

alpha_LSI  = sampleVector.alpha_LSI;
weight_LSI = sdata(Loc_LSI);
alpha_HSI  = sampleVector.alpha_HSI;
weight_HSI = sdata(Loc_HSI);
cirCode_LSI = [0.7 0.7 0.7]; % gray
cirCode_HSI = [0, 0, 0]; % 'r'
MarkAlpha = 1;
dim = 1;
%%
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
hpolarsctlowSI = polarscatter(alpha_LSI,weight_LSI,60,'ko','MarkerFaceColor',cirCode_LSI);
hpolarsctlowSI.MarkerFaceAlpha = MarkAlpha;
hold on;
hpolarscthighSI = polarscatter(alpha_HSI,weight_HSI,10,'ko','MarkerFaceColor',cirCode_HSI);
hpolarscthighSI.MarkerFaceAlpha = MarkAlpha;
hold on;
rlim([0 1]);
title('L-H SI (sample neuron)');
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig01');clf;
%%
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
[cossin_i_LSI, sample_mean_mag_LSI, sample_mean_rad_LSI] = circ_moment(alpha_LSI, weight_LSI, [], [], dim);
[cossin_i_HSI, sample_mean_mag_HSI, sample_mean_rad_HSI] = circ_moment(alpha_HSI, weight_HSI, [], [], dim);

hcLSI = polarscatter(sample_mean_rad_LSI, sample_mean_mag_LSI,100,'k^','MarkerFaceColor',cirCode_LSI);hold on;
hcHSI = polarscatter(sample_mean_rad_HSI, sample_mean_mag_HSI,100,'k^','MarkerFaceColor',cirCode_HSI);hold on;
rlim([0 0.07]);
hcHLdif = polarscatter(sample_mean_rad_LSI - sample_mean_rad_HSI, sample_mean_mag_LSI - sample_mean_mag_HSI,60,'k^','MarkerFaceColor','g');hold on;
title('mean resultant vector:L-H SI');hold off;
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig02');clf;
%% 
mean_rad_LSI_RQ = cell2mat(mean_rad_LSI);
mean_mag_LSI_RQ = cell2mat(mean_mag_LSI);
mean_rad_HSI_RQ = cell2mat(mean_rad_HSI);
mean_mag_HSI_RQ = cell2mat(mean_mag_HSI);

idx_SI_Full_XorYSIG_RQ   = [idx_SI_Full_XorYSIG.R;idx_SI_Full_XorYSIG.Q];

idx_HSI_reduce_XandY_RQ  = [idx_HSI_ftpXandYSIG.R;idx_HSI_ftpXandYSIG.Q];
idx_HSI_reduce_XorY_RQ   = [ idx_HSI_ftpXorYSIG.R; idx_HSI_ftpXorYSIG.Q];
idx_HSI_reduce_noneXY_RQ = ~idx_HSI_reduce_XorY_RQ;                      

idx_LSI_reduce_XandY_RQ  = [idx_LSI_ftpXandYSIG.R;idx_LSI_ftpXandYSIG.Q];
idx_LSI_reduce_XorY_RQ   = [ idx_LSI_ftpXorYSIG.R; idx_LSI_ftpXorYSIG.Q];
idx_LSI_reduce_noneXY_RQ = ~idx_LSI_reduce_XorY_RQ;                      

%% high or low SI blocks
dif_rad = mean_rad_LSI_RQ - mean_rad_HSI_RQ;
dif_mag = mean_mag_LSI_RQ - mean_mag_HSI_RQ;
dif_rad_XY_noneLHSI = dif_rad(idx_LSI_reduce_noneXY_RQ & idx_HSI_reduce_noneXY_RQ);
dif_mag_XY_noneLHSI = dif_mag(idx_LSI_reduce_noneXY_RQ & idx_HSI_reduce_noneXY_RQ);
dif_rad_XY_LSI_SIG = dif_rad(idx_LSI_reduce_XorY_RQ & ~idx_HSI_reduce_XorY_RQ);
dif_mag_XY_LSI_SIG = dif_mag(idx_LSI_reduce_XorY_RQ & ~idx_HSI_reduce_XorY_RQ);
dif_rad_XY_HSI_SIG = dif_rad(~idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
dif_mag_XY_HSI_SIG = dif_mag(~idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
dif_rad_twoXY_SIG = dif_rad(idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
dif_mag_twoXY_SIG = dif_mag(idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
%% histogram
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
hline3 = line([0 0],[0 120],'Color',[0.5 0.5 0.5],'LineStyle','--');hold on;
hmag_both = histogram(dif_mag_twoXY_SIG,'BinWidth',0.005,'BinLimits',[-0.08,0.08]);hold on;
hmag_both.FaceColor = [25, 137, 100] ./ 255;
hmag_both.FaceAlpha = 0.7;

hmag_none = histogram(dif_mag_XY_noneLHSI,'BinWidth',0.005,'BinLimits',[-0.08,0.08]);hold on;
hmag_none.FaceColor = [0.3 0.3 0.3];
hmag_none.FaceAlpha = 0.5;

hmag_LSI_only = histogram(dif_mag_XY_LSI_SIG,'BinWidth',0.005,'BinLimits',[-0.08,0.08]);hold on;
hmag_LSI_only.FaceColor = 'c';
hmag_LSI_only.FaceAlpha = 0.7;

hmag_HSI_only = histogram(dif_mag_XY_HSI_SIG,'BinWidth',0.005,'BinLimits',[-0.08,0.08]);hold on;
hmag_HSI_only.FaceColor = 'b';
hmag_HSI_only.FaceAlpha = 0.7;
%
xxboth = diff(hmag_both.BinEdges)/2 + hmag_both.BinEdges(1:(end-1));
[sigmaNew,muNew,Anew]=mygaussfit(xxboth,hmag_both.Values);
yfitbothSI=Anew*exp(-(xxboth-muNew).^2/(2*sigmaNew^2));
plot(xxboth,yfitbothSI,'-g','LineWidth',3,'Color',[25, 137, 100] ./ 255);hold on;
%
xxnone = diff(hmag_none.BinEdges)/2 + hmag_none.BinEdges(1:(end-1));
[sigmaNew,muNew,Anew]=mygaussfit(xxnone,hmag_none.Values);
yfitNoneSI=Anew*exp(-(xxnone-muNew).^2/(2*sigmaNew^2));
plot(xxnone,yfitNoneSI,'-g','LineWidth',3,'Color',[0.2 0.2 0.2]);hold on;
%
xxLSIonly = diff(hmag_LSI_only.BinEdges)/2 + hmag_LSI_only.BinEdges(1:(end-1));
[sigmaNew,muNew,Anew]=mygaussfit(xxLSIonly,hmag_LSI_only.Values);
yfitLSI=Anew*exp(-(xxLSIonly-muNew).^2/(2*sigmaNew^2));
plot(xxLSIonly,yfitLSI,'-g','LineWidth',3,'Color','c');hold on;
%
xxHSIonly = diff(hmag_HSI_only.BinEdges)/2 + hmag_HSI_only.BinEdges(1:(end-1));
[sigmaNew,muNew,Anew]=mygaussfit(xxHSIonly,hmag_HSI_only.Values);
yfitHSI=Anew*exp(-(xxHSIonly-muNew).^2/(2*sigmaNew^2));
plot(xxHSIonly,yfitHSI,'-g','LineWidth',3,'Color','b');hold on;
%
legend([hmag_both, hmag_none, hmag_LSI_only, hmag_HSI_only], {'Both', 'None', 'LSI only', 'HSI only'});
title('L-H SI (Magnitude)');
xlabel('Difference of magnitude between LSI and HSI');
set(gca,'tickdir','out');
axis([-0.1, 0.1, 0, 120]);axis square;box off;
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig03');clf;

%% raincloud plot
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
hline4 = line([0 0],[-60, 40],'Color',[0.5 0.5 0.5],'LineStyle','--');hold on;
plevel = 0.01;
[cb] = cbrewer('qual', 'Set3', 12, 'pchip');
h1 = raincloud_plot(dif_mag_twoXY_SIG, 'box_on', 1, 'color', [25, 137, 100] ./ 255, 'alpha', 1,...
     'box_dodge', 1, 'box_dodge_amount', .15, 'dot_dodge_amount', .15,'box_col_match', 0);
h2 = raincloud_plot(dif_mag_XY_noneLHSI, 'box_on', 1, 'color', [0.3 0.3 0.3], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .55, 'dot_dodge_amount', .55, 'box_col_match', 0);
h3 = raincloud_plot(dif_mag_XY_LSI_SIG, 'box_on', 1, 'color', [0 1 1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', .95, 'dot_dodge_amount', .95, 'box_col_match', 0);
h4 = raincloud_plot(dif_mag_XY_HSI_SIG, 'box_on', 1, 'color', [0 0 1], 'alpha', 0.5,...
     'box_dodge', 1, 'box_dodge_amount', 1.35, 'dot_dodge_amount', 1.35, 'box_col_match', 0);
legend([h1{1} h2{1} h3{1} h4{1}], {'Both', 'None', 'LSI only', 'HSI only'});
title('L-H SI (Magnitude)');
xlabel('Difference of magnitude between LSI and HSI');
ylabel('Number of neurons');
set(gca,'tickdir','out','ytick',0:10:50);
[~,P_both,CI_both,stats_both] = ttest(dif_mag_twoXY_SIG);
[~,P_none,CI_none,stats_none] = ttest(dif_mag_XY_noneLHSI);
[~,P_LSIonly,CI_LSIonly,stats_LSIonly] = ttest(dif_mag_XY_LSI_SIG);
[~,P_HSIonly,CI_HSIonly,stats_HSIonly] = ttest(dif_mag_XY_HSI_SIG);
text(-0.08,-8,sprintf('t(%.d)=%.2f, p=%.3f',stats_both.df,stats_both.tstat,P_both));
text(-0.08,-18,sprintf('t(%.d)=%.2f, p=%.3f',stats_none.df,stats_none.tstat,P_none));
text(-0.08,-28,sprintf('t(%.d)=%.2f, p=%.2e',stats_LSIonly.df,stats_LSIonly.tstat,P_LSIonly));
text(-0.08,-38,sprintf('t(%.d)=%.2f, p=%.2e',stats_HSIonly.df,stats_HSIonly.tstat,P_HSIonly));
axis([-0.1, 0.1, -60, 40]);box off;axis square;
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig04');clf;
%%
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
[   mean_twoXY, ~,    sem_twoXY]  = Bhvstats(   dif_mag_twoXY_SIG, 1);
[ mean_noneLHSI, ~, sem_noneLHSI] = Bhvstats( dif_mag_XY_noneLHSI, 1);
[  mean_onlyLSI, ~,  sem_onlyLSI] = Bhvstats(  dif_mag_XY_LSI_SIG, 1);
[  mean_onlyHSI, ~,  sem_onlyHSI] = Bhvstats(  dif_mag_XY_HSI_SIG, 1);

hbarSem = errorbar(1:4,[mean_twoXY,mean_onlyLSI,mean_onlyHSI,mean_noneLHSI],[sem_twoXY,sem_onlyLSI,sem_onlyHSI,sem_noneLHSI],'k.');hold on;
hbarMean = bar(1:4,[mean_twoXY,mean_onlyLSI,mean_onlyHSI,mean_noneLHSI]);hold on;
ylabel('Difference of magnitude between LSI and HSI');
set(gca,'tickdir','out','xtick',1:1:4,'xticklabel',{'Both','LSI only','HSI only','None'},...                        
                        'ytick',-0.01:0.005:0.01,'yticklabel',{'-0.01','-0.005','0','0.005','0.01'});
text(0.7,-0.01,sprintf('n=%d',length(dif_mag_twoXY_SIG)),'FontSize',12);
text(1.7,-0.01,sprintf('n=%d',length(dif_mag_XY_noneLHSI)),'fontSize',12);
text(2.7,-0.01,sprintf('n=%d',length(dif_mag_XY_LSI_SIG)),'fontSize',12);
text(3.7,-0.01,sprintf('n=%d',length(dif_mag_XY_HSI_SIG)),'fontSize',12);
title('L-H SI: Magnitude');
axis([0.5,4.5,-0.013,0.013]);axis square;box off;
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig05');clf;
%% conver theta to axial common scale
deg_dif_rad = rad2deg(dif_rad);
idx_n180 = deg_dif_rad < -180;
idx_p180 = deg_dif_rad > 180;
deg_dif_rad(idx_n180) = deg_dif_rad(idx_n180)+360;
deg_dif_rad(idx_p180) = deg_dif_rad(idx_p180)-360;
axial_dif_rad_XY_noneSI = deg_dif_rad(idx_LSI_reduce_noneXY_RQ & idx_HSI_reduce_noneXY_RQ);
axial_dif_rad_LSI_XY_SIG = deg_dif_rad(idx_LSI_reduce_XorY_RQ & ~idx_HSI_reduce_XorY_RQ);
axial_dif_rad_HSI_XY_SIG = deg_dif_rad(~idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
axial_dif_rad_allSIXY_SIG = deg_dif_rad(idx_LSI_reduce_XorY_RQ & idx_HSI_reduce_XorY_RQ);
% %
figure('paperorientation','landscape','units','normalized','position',[0.05 0.05 0.9 0.85]);
hallSI = histogram(axial_dif_rad_allSIXY_SIG,-180:15:180);hold on;
hallSI.FaceColor = [25, 137, 100] ./ 255; % [255,105,180]./255; % pink rgb([255,192,203]./255; [255,105,180]./255
hallSI.FaceAlpha = 0.7;
hnoneSI= histogram(axial_dif_rad_XY_noneSI,-180:15:180);hold on;
hnoneSI.FaceColor = 'k';
hnoneSI.FaceAlpha = 0.3;
hLSI = histogram(axial_dif_rad_LSI_XY_SIG,-180:15:180);hold on;
hLSI.FaceColor = 'b';
hLSI.FaceAlpha = 0.7;
hHSI = histogram(axial_dif_rad_HSI_XY_SIG,-180:15:180);hold on;
hHSI.FaceColor = 'c';
hHSI.FaceAlpha = 0.7;

title('L-H SI (spatial tuning)');
xlabel(['Low SI - High SI' ' ('  '''d' ')']);
ylabel('Number of neurons');
set(gca,'xtick',-180:60:180,'tickdir','out');
axis([-180 180 0 70]);
axis square;box off;
% print(gcf,'-dpdf','-r300','-fillpage','SFWM_50_spatial tuning_reducedModel_fig06');
%% fitting curve
x = (-170:15:180)-5;
yNoneSI = hnoneSI.Values;
yLSI = hLSI.Values;      
yHSI = hHSI.Values;      
yallSI = hallSI.Values;  
[sigmaNew,muNew,Anew]=mygaussfit(x,yallSI);
yfitallSI=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
% 
[sigmaNew,muNew,Anew]=mygaussfit(x,yLSI);
yfitLSI=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
% 
[sigmaNew,muNew,Anew]=mygaussfit(x,yHSI);
yfitHSI=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
% 
[sigmaNew,muNew,Anew]=mygaussfit(x,yNoneSI);
yfitNoneSI=Anew*exp(-(x-muNew).^2/(2*sigmaNew^2));
hold on; 
plot(x,yfitallSI,'-g','LineWidth',3,'Color',[25, 137, 100] ./ 255);hold on;
plot(x,yfitNoneSI,'-k','LineWidth',3,'Color',[0.2 0.2 0.2]);hold on;
plot(x,yfitLSI,'-b','LineWidth',3);hold on;
plot(x,yfitHSI,'-c','LineWidth',3);hold on;
legend('Both*** (n = 548)','n.s. (n = 281)','LSI only (n = 123)','HSI only (n = 125)');
axis([-180 180 0 90]);
axis square;box off;
% print(gcf,'-dpdf','-r300','-fillpage','spatial tuning_reducedModel_fig06');clf;close all;