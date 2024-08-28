% Load relevant variables
addpath('A:\eeglab_OldLSL_DataAna04072023')
addpath(genpath('Scripts'))
eeglab
close all;
clearvars;
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization');
% High low cons for 243 subjects
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\Pre_TwoTap_src_HighConsis243.mat')
source_all_high = Pre_TT_srcmap_HighConsis;
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\Pre_TwoTap_src_LowSlowConsis243.mat')
source_all_low = Pre_TT_srcmap_LowSlowConsis;
% High low cons for 81 subs
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\Pre_TwoTap_src_HighConsis.mat')
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\Pre_TwoTap_src_LowSlowConsis.mat')
% Resting data for subs
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\Pre_Rest_TT_March2023.mat')
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\tt_rest_outlierrm.mat')
load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\adl_rest_outlierrm.mat')
% tmp=xlsread('/Users/Dhak/Documents/2Tapfinaldata/Final_src/neuro324.xlsx');
% age=tmp(:,2);

savePath = 'A:\TwoTap\Manuscript\Figure';

%%
source_all_high2=cat(4,source_all_high,Pre_TT_srcmap_HighConsis);
source_all_low2=cat(4,source_all_low, Pre_TT_srcmap_LowSlowConsis);
rest=cat(4,tt_rest,adl_rest,Pre_Rest_srcmap);
%%
% missingrest=[257,258,272,287,289,291]
%source_all_high2(:,:,:,missingrest)=[];
%source_all_low2(:,:,:,missingrest)=[];
size(source_all_high2)
size(source_all_low2)
size(rest)
% size(age)

%% get rid of outliers
tmp=squeeze(nanmean(nanmean(nanmean(source_all_high2,1),2),3));
tmp_std=5*std(tmp);
tmp2=squeeze(nanmean(nanmean(nanmean(source_all_low2,1),2),3));
tmp2_std=5*std(tmp2);

tmp3=squeeze(nanmean(nanmean(nanmean(rest,1),2),3));
tmp3_std=5*std(tmp3);

outl1=find(tmp2>tmp2_std);
outl2=find(tmp>tmp_std);
outl3=find(tmp3>tmp3_std);
%%
outl=[outl1;outl2;outl3];

source_all_low=source_all_low2;
source_all_high=source_all_high2;

source_all_low(:,:,:,[outl])=[];
source_all_high(:,:,:,[outl])=[];
rest(:,:,:,[outl])=[];
% age(outl)=[];

size(source_all_low)
size(source_all_high)
size(rest)
% size(age)

%% For network-based analyses 
netwrk(1).name = 'FPN';
netwrk(1).roi = [5,6,55,56,59,60];
netwrk(2).name = 'CON';
netwrk(2).roi = [3,4,19,20,37,38,39,40,41,42,57,58,67,68];
netwrk(3).name = 'aDMN';
netwrk(3).roi = [11,12,25,26,29,30,53,54];
netwrk(4).name = 'pDMN';
netwrk(4).roi = [15,16,21,22,51,52];
netwrk(5).name = 'MTL-DMN';
netwrk(5).roi = [9,10,17,18,31,32,35,36,61,62,65,66];

%% organizing resting data into temporal epochs for analysis.

trest=-3996:4:0;
fr=[];
for a = 1:3%2
    tw_prer=[find(trest>-4000 & trest<0000)];
    tw_prer=tw_prer(1,[1,end]);
    f(a).r=squeeze(nanmean(rest(a,:,tw_prer(1):tw_prer(2),:),3));
   
    [~,p]=ttest(f(a).r');
    f(a).rp=p;
    f(a).rpc=fdr(p);
          
    for b = 1:5
        f(a).netr(:,b)=nanmean(f(a).r(netwrk(b).roi,:),1)';
    end
end

%%

% fdr_correction of all brain diff
tmpp1=[f(1).rp,f(2).rp,f(3).rp];
p_all1=fdr(tmpp1);
f(1).rpc2=p_all1(1:68);
f(2).rpc2=p_all1(69:136);
f(3).rpc2=p_all1(137:end);
close all;
clc;

%% organizing task data into temporal epochs for analysis.

t=-3996:4:3996;

for a = 1:3
    tw_pre=[find(t>-4000 & t<000)];
    tw_pre=tw_pre(1,[1,end]);
    f(a).high_pre=squeeze(nanmean(source_all_high(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).low_pre=squeeze(nanmean(source_all_low(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).df_pre=(f(a).low_pre-f(a).high_pre);
    f(a).highRest_pre = f(a).high_pre-f(a).r;
    
    [~,p]=ttest(f(a).df_pre');
    f(a).dfp=p;
    f(a).dfpc=fdr(p);
    
    [h,p]=ttest(f(a).high_pre');
    f(a).highp=p;
    f(a).highpc=fdr(p);
    
    [h,p]=ttest(f(a).low_pre');
    f(a).lowp=p;
    f(a).lowpc=fdr(p);

    [h,p]=ttest(f(a).highRest_pre');
    f(a).highRestp=p;
    f(a).highRestpc=fdr(p);

    
    for b = 1:5
        f(a).netlotime(:,:,b)=squeeze(nanmean(source_all_low(a,netwrk(b).roi,:,:),2)); 
        f(a).nethitime(:,:,b)=squeeze(nanmean(source_all_high(a,netwrk(b).roi,:,:),2)); 
        f(a).netdftime(:,:,b)=f(a).netlotime(:,:,b)-f(a).nethitime(:,:,b);
         
        f(a).netlotimeav(:,b)=squeeze(nanmean(squeeze(nanmean(source_all_low(a,netwrk(b).roi,:,:),2)),2)); 
        f(a).nethitimeav(:,b)=squeeze(nanmean(squeeze(nanmean(source_all_high(a,netwrk(b).roi,:,:),2)),2)); 
        f(a).netdftimeav(:,b)=squeeze(nanmean(f(a).netdftime(:,:,b),2)); 

        f(a).netlotimesem(:,b)=squeeze(nanstd(squeeze(nanmean(source_all_low(a,netwrk(b).roi,:,:),2)),[],2))/sqrt(size(source_all_low,3));  
        f(a).nethitimesem(:,b)=squeeze(nanstd(squeeze(nanmean(source_all_high(a,netwrk(b).roi,:,:),2)),[],2))/sqrt(size(source_all_low,3));  
        f(a).netdftimesem(:,b)=squeeze(nanstd(f(a).netdftime(:,:,b),[],2))/sqrt(size(source_all_low,3));
        
        f(a).netlo(:,b)=nanmean(f(a).low_pre(netwrk(b).roi,:),1)';
        f(a).nethi(:,b)=nanmean(f(a).high_pre(netwrk(b).roi,:),1)';
        f(a).netdf(:,b)=nanmean(f(a).df_pre(netwrk(b).roi,:),1)';
    end
end

% fdr_correction of all brain diff
tmpp1=[f(2).dfp];
p_all1=fdr(tmpp1);
%f(1).dfpc2=p_all1(1:68);
f(2).dfpc2=p_all1(1:68);
%f(3).dfpc2=p_all1(137:end);

tmpp2=[f(2).highp];
p_all2=fdr(tmpp2);
%f(1).highpc2=p_all2(1:68);
f(2).highpc2=p_all2(1:68);
%f(3).highpc2=p_all2(137:end);

tmpp3=[f(2).lowp];
p_all3=fdr(tmpp3);
%f(1).lowpc2=p_all3(1:68);
f(2).lowpc2=p_all3(1:68);
%f(3).lowpc2=p_all3(137:end);
close all;
clc;
%% Save Variables

T=table([f(2).netlo,f(2).nethi, f(2).netdf f(2).netr]);
writetable(T,[savePath '\tempData4Sec.xlsx'],'Sheet','networkData');

T=table([nanmean(f(2).r,2)]);
writetable(T,[savePath '\tempData4Sec.xlsx'],'Sheet','RestData');


%% Plot network data 

data = [nanmean(f(2).netlo,1) ;nanmean(f(2).nethi,1)];
sem = [nanstd(f(2).netlo) ;nanstd(f(2).nethi)]/sqrt(324);
close all
figure

subplot(121)
hold on
% plotErrBar(data',sem',{'b','r'})
allLow = reshape(f(2).netlo,[],1);
allHigh = reshape(f(2).nethi,[],1);
allRest = reshape(f(2).netr,[],1);

netOrder =(repmat(reshape(repmat([1:5],[324],1),[],1),3,1));
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
bbb=boxchart(netOrder, [allLow;allHigh; allRest],'GroupByColor',[ones(1,324*5), 2*ones(1,324*5), 3*ones(1,324*5)]','MarkerStyle','none')

x = repmat(1:5,324,1)-.33;
x=reshape(x,[],1);
swarmchart(x,[allLow],20,'.','XJitterWidth',0.1,'MarkerEdgeColor',[0 0.4470 0.7410]);
x = repmat(1:5,324,1);
x=reshape(x,[],1);
swarmchart(x,[allHigh],20,'.','XJitterWidth',0.1,'MarkerEdgeColor',[0.8500 0.3250 0.0980]);
x = repmat(1:5,324,1)+0.33;
x=reshape(x,[],1);
swarmchart(x,[allRest],20,'.','XJitterWidth',0.1,'MarkerEdgeColor',[0.9290 0.6940 0.1250]);


ylim([-0.0001 0.018])

set(gca,'XTickLabel',{netwrk.name})

legend('Low Cons','High Cons','Rest','box','off','location','northwest','fontsize',12,'box','off')
% set(gca,'XGrid','off','YGrid','on','Position',[0.1300 0.1100 0.7750 0.8121])
% set(gcf,'position',[811 609 774 344.6667])
% saveas(gcf,'A:\TwoTap\Manuscript\Figure\fig3_network.tiff')

% Plot Network Difference
data = [nanmean(f(2).netdf,1) ];
sem = [nanstd(f(2).netdf)]/sqrt(324);
% close all
% figure


%%% Difference box plots
subplot(122)
hold on
% plotErrBar(data',sem',{'k'})

x = repmat(1:5,324,1);
boxchart(f(2).netdf,'MarkerStyle','none')
swarmchart(x,f(2).netdf,10,'k.','XJitterWidth',0.2)

set(gca,'XTickLabel',{netwrk.name})
ylim([-4E-3 12E-3])
% xlim([-.0 5.9])
% legend('Dist - Atten','box','off','location','northwest','fontweight','bold','fontsize',12)
% set(gca,'XGrid','off','YGrid','on','Position',[0.1300 0.1100 0.7750 0.8121])
% set(gcf,'position',[811 609 774 344.6667])
% saveas(gcf,'A:\TwoTap\Manuscript\Figure\fig3_networkdiff.tiff')

%% Fig 2A: brain activity for different contrasts/conditions
close all;
for i =2
hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = nanmean(f(i).high_pre');
var(f(i).highpc>0.05)=0;
plot68roi(hm, T'*var',[-.001 .001],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\2Battended.tiff']);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = nanmean(f(i).low_pre');
var(f(i).lowpc>0.05)=0;
plot68roi(hm, T'*var',[-.001 .001],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\2Bdistracted.tiff']);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = nanmean(f(i).r');
var(f(i).rpc>0.05)=0;
plot68roi(hm, T'*var',[-0.001 .001],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\2Brest.tiff']);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = nanmean(f(i).df_pre');
var(f(i).dfpc>0.05)=0;
plot68roi(hm, T'*var',[-0.001 .001],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\2Ddiff.tiff']);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = nanmean(f(i).highRest_pre');
var(f(i).highRestpc>0.05)=0;
plot68roi(hm, T'*var',[-0.001 .001],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\2Dhighrestdf.tiff']);

end
%% Save Roi Data
sfile=[savePath '\tempData4Sec.xlsx'];

T=table();
T.attended = nanmean(f(2).high_pre,2)';
T.distracted = nanmean(f(2).low_pre,2)';
T.df = nanmean(f(2).df_pre,2)';
T.rest = nanmean(f(2).r,2)';

writetable(T,sfile,'Sheet','Subject_data');

%% Fig 2C CV across the brain.
clear CV*
for i = 1
    CV_lo(i,:)=nanstd(f(2).low_pre,[],2)./nanmean(f(2).low_pre,2);
    CV_hi(i,:)=nanstd(f(2).high_pre,[],2)./nanmean(f(2).high_pre,2);
    CV_rest(i,:)=nanstd(f(2).r,[],2)./nanmean(f(2).r,2);
end

CV_lom=nanmean(CV_lo');
CV_him=nanmean(CV_hi');
CV_restm=nanmean(CV_rest');

CV_losem=nanstd(CV_lo')/sqrt(68);
CV_hisem=nanstd(CV_hi')/sqrt(68);
CV_restsem=nanstd(CV_rest')/sqrt(68);

CVmall=[CV_lom;CV_him;CV_restm];
CVsemall=[CV_losem;CV_hisem;CV_restsem];

CVmall=CVmall';
CVsemall=CVsemall';

%% Plot CV
% figure;
% labels={'distracted','attended','rest'};
% xticklabels(labels);

close all
figure
hold on
CVplot = [CV_hi ;CV_lo; CV_rest]';
x = repmat(1:3,68,1);
boxchart(CVplot,'MarkerStyle','none')
swarmchart(x,CVplot,50,'k.')
ylabel('CV')
set(gca,'xticklabels',{'Attended', 'Distracted','Rest'},'fontweight','bold','fontsize',16,'xticklabelrotation',45)
set(gcf,'position',[1.1883e+03 773 371.6667 565])
saveas(gcf,[savePath '\2CCV.tiff']);

%% Save CV data
sfile=[savePath '\tempData4Sec.xlsx'];
T=table(CV_lo',CV_hi',CV_rest','VariableNames',{'CV_LowSlow','CV_attended','CV_rest'});
writetable(T,sfile,'Sheet','CV');

%% Fig 2C resting state and task activation correlations
wb=0; % set to 1 if you want to do whole-brain plots

close all
a=2;
        x1=f(a).r;
        x2=f(a).low_pre;
        x3=f(a).high_pre;
%        
    for i = 1:68
   
        [lowrr(i),lowrp(i)]=corr(x1(i,:)',x2(i,:)','Type','Spearman','Rows','pairwise');
        [highrr(i),highrp(i)]=corr(x1(i,:)',x3(i,:)','Type','Spearman','Rows','pairwise');
        [highlowr(i),highlowp(i)]=corr(x2(i,:)',x3(i,:)','Type','Spearman','Rows','pairwise');
        
        [lowrr2(i)]=partialcorr(x1(i,:)',x2(i,:)',age,'Type','Spearman','Rows','pairwise');
        [highrr2(i)]=partialcorr(x1(i,:)',x3(i,:)',age,'Type','Spearman','Rows','pairwise');
        [highlowr2(i)]=partialcorr(x2(i,:)',x3(i,:)',age,'Type','Spearman','Rows','pairwise');
        %
        lowrpc=fdr(lowrp);
        highrpc=fdr(highrp);
        highlowpc=fdr(highlowp);
    end

lowrr_alpha=squeeze(atanh(lowrr))';
highrr_alpha=squeeze(atanh(highrr))';
highlowr_alpha=squeeze(atanh(highlowr))';

lowrr_alpha2=squeeze(atanh(lowrr2))';
highrr_alpha2=squeeze(atanh(highrr2))';
highlowr_alpha2=squeeze(atanh(highlowr2))';




%% PLotting ISC bars
close all
figure
hold on
ISCplot = [highrr_alpha lowrr_alpha highlowr_alpha];
x = repmat(1:3,68,1);
boxchart(ISCplot)
swarmchart(x,ISCplot,50,'k.')
ylabel('Fisher Z')
set(gca,'xticklabels',{'Atten-Rest', 'Dist-Rest','Atten-Dist'},'xticklabelrotation',45)
set(gcf,'position',[1.1883e+03 953 328.7000 385])
saveas(gcf,[savePath '\2CISC.tiff']);
    
%%
T=table(lowrr_alpha2,highrr_alpha2,highlowr_alpha2);
sfile=[savePath '\tempData4Sec.xlsx'];
writetable(T,sfile,'Sheet','ISC_age');
%% Plotting roi by task

var1 = nanmean(f(2).high_pre');
var2 = nanmean(f(2).low_pre');
var3 = nanmean(f(2).r');

close all
figure
hold on
ROIsplot = [var1' var2' var3'];
ROIsplot=filloutliers(ROIsplot,nan,1);
x = repmat(1:3,68,1);
boxchart(ROIsplot,'markerstyle','none')
swarmchart(x,ROIsplot,50,'k.')
ylabel('Activity')
set(gca,'xticklabels',{'Atten', 'Dist','rest'},'xticklabelrotation',45)
set(gcf,'position',[1.1883e+03 953 328.7000 385])
saveas(gcf,[savePath '\2Droi_outRMV.tiff']);

%% Supp Fig 1: whole-brain correlation maps (data not used).
%  if wb==1
% hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
% T = hm.indices4Structure(hm.atlas.label);
% T = double(T)';
% P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
% %set(0,'DefaultFigureWindowStyle','docked');
% labelnames = {'low - rest corr'};
% var = lowrr';
% %var(p5d>0.05)=0;
% plotMultipanel_v2_jet_indBrainMaps_modifcorr(hm, T'*var',[-.6 .6],labelnames);
% set(gcf,'Position',[1,1,267,476]); 
% 
% hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
% T = hm.indices4Structure(hm.atlas.label);
% T = double(T)';
% P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
% %set(0,'DefaultFigureWindowStyle','docked');
% labelnames = {'high - rest corr'};
% var = highrr';
% %var(p5d>0.05)=0;
% plotMultipanel_v2_jet_indBrainMaps_modifcorr(hm, T'*var',[-.6 .6],labelnames);
% set(gcf,'Position',[1,1,267,476]); 
% 
% hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
% T = hm.indices4Structure(hm.atlas.label);
% T = double(T)';
% P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
% %set(0,'DefaultFigureWindowStyle','docked');
% labelnames = {'low - high corr'};
% var = highlowr';
% %var(p5d>0.05)=0;
% plotMultipanel_v2_jet_indBrainMaps_modifcorr(hm, T'*var',0,labelnames);
% set(gcf,'Position',[1,1,267,476]); 
%  end
% end
% 
% lowrr_alpha=squeeze(atanh(lowrr(2,:)))';
% highrr_alpha=squeeze(atanh(highrr(2,:)))';
% highlowr_alpha=squeeze(atanh(highlowr(2,:)))';
% T=table(lowrr_alpha,highrr_alpha,highlowr_alpha);
% sfile=['/Users/Dhak/Documents/2Tapfinaldata/June2023/figdata/Fig2C.xlsx'];
% writetable(T,sfile,'Sheet','ISC');
% 
% %% Network Correlations
% close all;
% 
% for a = 1:3
%     for i = 1:5
%         x1=f(a).netr(:,i);
%         x2=f(a).netlo(:,i);
%         x3=f(a).nethi(:,i);
%         
%     [netlorr(a,i),netlorp(a,i)]=corr(x1,x2,'Rows','pairwise','Type','Spearman');
%     [nethirr(a,i),nethirp(a,i)]=corr(x1,x3,'Rows','pairwise','Type','Spearman');
%     [nethilor(a,i),nethilop(a,i)]=corr(x2,x3,'Rows','pairwise','Type','Spearman');
%     end
% end
% 
% netlorpc=fdr(netlorp);
% nethirpc=fdr(nethirp);
% nethilopc=fdr(nethilop);

%% Supp Fig 2B Network NB Correlations

% labels2={'FPN','CON','aDMN','pDMN','mtDMN'};
% 
% for a = 1:3
%         netrestcorr{a}(1,:)=netlorr(a,:);
%         netrestcorr{a}(2,:)=nethirr(a,:);
%         netrestcorr{a}(3,:)=nethilor(a,:);
% end
% subplot(131)
% bar(netrestcorr{1}');
% axis([-inf inf 0 1]);
% title('theta');
% xticklabels(labels2);
% %xtickangle(90);
% box off;
% 
% subplot(132)
% bar(netrestcorr{2}');
% axis([-inf inf 0 1]);
% title('alpha');
% xticklabels(labels2);
% %xtickangle(90);
% box off;
% 
% subplot(133)
% bar(netrestcorr{1}');
% axis([-inf inf 0 1]);
% title('beta');
% xticklabels(labels2);
% %xtickangle(90);
% box off;
% 
% legend('rest - lo corr','rest - hi corr','lo - high corr');
% 
