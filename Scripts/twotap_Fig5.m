% Load relevant variables

% clearvars;
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

source_all_high2=cat(4,source_all_high,Pre_TT_srcmap_HighConsis);
source_all_low2=cat(4,source_all_low, Pre_TT_srcmap_LowSlowConsis);
rest=cat(4,tt_rest,adl_rest,Pre_Rest_srcmap);
%% get rid of outliers
tmp=squeeze(nanmean(nanmean(nanmean(source_all_high2,1),2),3));
tmp_std=5*std(tmp);
tmp2=squeeze(nanmean(nanmean(nanmean(source_all_low2,1),2),3));
tmp2_std=5*std(tmp2);


outl1=find(tmp2>tmp2_std);
outl2=find(tmp>tmp_std);
outl=[outl1;outl2];

source_all_low=source_all_low2;
source_all_high=source_all_high2;

source_all_low(:,:,:,[outl])=[];
source_all_high(:,:,:,[outl])=[];

size(source_all_low)
size(source_all_high)

%%
behdata=readtable("A:\TwoTap\Manuscript\TwoTap_ms_finaldata_May29.xlsx");
%%
% dep=table2array(behdata(:,4));
dep = behdata.Depression; 
% anx=table2array(behdata(:,3));
% mind=table2array(behdata(:,5));
mind = behdata.Mindfulness;
% beh=table2array(behdata(:,12));
% statemind=table2array(behdata(:,35));
statemind = behdata.state_mindfulness;
% age=table2array(behdata(:,2));
age = behdata.age;

% dep(outl)=[];
% anx(outl)=[];
% mind(outl)=[];
% beh(outl)=[];
% statemind(outl)=[];
% yall=[anx,dep,mind,beh];
% age(outl)=[];

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
netwrk(6).name = 'Visual';
netwrk(6).roi =[7 8 13 14 23 24 27 28 43 44];

%% Fig 4A Mindfulness Correlation
% organizing data into temporal epochs for analysis.

t=-3996:4:3996;
f=[];
for a = 1:3
    tw_pre=[find(t>-4000 & t<000)];
    tw_pre=tw_pre(1,[1,end]);
    f(a).high_pre=squeeze(nanmean(source_all_high(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).low_pre=squeeze(nanmean(source_all_low(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).df_pre=(f(a).low_pre-f(a).high_pre);

    f(a).high_prets=squeeze((source_all_high(a,:,tw_pre(1):tw_pre(2),:)));
    f(a).low_prets=squeeze((source_all_low(a,:,tw_pre(1):tw_pre(2),:)));

    [~,p]=ttest(f(a).df_pre');
    f(a).dfp=p;
    f(a).dfpc=fdr(p);
    
    [h,p]=ttest(f(a).high_pre');
    f(a).highp=p;
    f(a).highpc=fdr(p);
    
    [h,p]=ttest(f(a).low_pre');
    f(a).lowp=p;
    f(a).lowpc=fdr(p);
    
    for b = 1:5
        f(a).netlo(:,b)=nanmean(f(a).low_pre(netwrk(b).roi,:),1)';
        f(a).nethi(:,b)=nanmean(f(a).high_pre(netwrk(b).roi,:),1)';
        f(a).netdf(:,b)=nanmean(f(a).df_pre(netwrk(b).roi,:),1)';
        
        f(a).netlots(:,:,b)=squeeze(nanmean(f(a).low_prets(netwrk(b).roi,:,:),1));
        f(a).nethits(:,:,b)=squeeze(nanmean(f(a).high_prets(netwrk(b).roi,:,:),1));
    end
end
%% Plotting mindfulness corr
close all;
a = 2;
% MH = behdata.Consistency;
MH = mind;
% MH(age>60)=nan;

tmp5h=f(a).high_pre;
tmp5l=f(a).low_pre;
% tmp5h = filloutliers(tmp5h, nan);
% tmp5l = filloutliers(tmp5l, nan);
tmp5df=tmp5l-tmp5h;
% tmp5l=f(a).netlo';
% tmp5h=f(a).nethi';
% tmp5df=f(a).netdf';

% tmp5l=[];tmp5h=[];tmp5df=[];
% for n=1:5 
%     tmph=f(2).high_pre;
%     tmpl=f(2).low_pre;
%     tmph = filloutliers(tmph, nan);
%     tmpl = filloutliers(tmpl, nan);
% 
%     tmp5l(n,:) = nanmean(tmpl(netwrk(n).roi,:),1);
%     tmp5h(n,:) = nanmean(tmph(netwrk(n).roi,:),1);
%     tmp5df(n,:)= nanmean(tmpl(netwrk(n).roi,:),1)-nanmean(tmph(netwrk(n).roi,:),1);
% end

[mindr1,mindp1]=corr(MH,tmp5l','Type','Spearman','Rows','pairwise');
[mindr2,mindp2]=corr(MH,tmp5h','Type','Spearman','Rows','pairwise');
[mindr3,mindp3]=corr(MH,tmp5df','Type','Spearman','Rows','pairwise');
mindr=[mindr1' mindr2', mindr3'];
mindp=[mindp1' mindp2', mindp3'];
mindpc=[];
% mindpc=fdr(mindp);
% mindpc(:,1) = fdr(mindp(:,1));
% mindpc(:,2) = fdr(mindp(:,2));
% mindpc(:,3) = fdr(mindp(:,3));
mindpc=mindp;

%%% FDR Within network only
mindpc=ones(68,3);
for n=1:5
    mindpc(netwrk(n).roi,1) = fdr(mindp(netwrk(n).roi,1));
    mindpc(netwrk(n).roi,2) = fdr(mindp(netwrk(n).roi,2));
    mindpc(netwrk(n).roi,3) = fdr(mindp(netwrk(n).roi,3));
end

%%% reverse network
% mindtmp = zeros(68,3); 
% mindptmp= ones(68,3);
% for n=1:5
%     mindtmp(netwrk(n).roi,:) = repmat(mindr(n,:),length(netwrk(n).roi),1);
%     mindptmp(netwrk(n).roi,:) = repmat(mindpc(n,:),length(netwrk(n).roi),1);
% end
% mindr = mindtmp;
% mindpc = mindptmp;
%%%

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'dist','atten','diff'};
var = mindr;
var(mindpc >0.05)=0;
plot68roi(hm, T'*var,[],labelnames);
set(gcf,'Position',[965.6667 352.3333 820.6667 1.0087e+03]);
% saveas(gcf, savePath +"/Fig5aNB.tiff")

%% Network Connectivity corr with MH
MH = mind;
% ya - from twotap_Fig4 aDMN conn with distracted FPN distracted CON attended FPN attended CON 
% ym - from twotap_Fig4 mtDMN conn with distracted FPN distracted CON attended FPN attended CON 
% yp - from twotap_Fig4 pDMN conn with distracted FPN distracted CON attended FPN attended CON 
yp2 = [pDMNnetlocorrm(:,6:10) pDMNnethicorrm(:,6:10) pDMNnetrestm(:,6:10)]; % FPN-CON corrs

netcondata = yp;
% netcondata = netcondata(:,1:2)-netcondata(:,3:4);
clc
% [mindr1,mindp1]=corr(MH,netcondata,'Type','Spearman','Rows','pairwise')


% mindpc(:,1) = fdr(mindp(:,1));
% mindpc(:,2) = fdr(mindp(:,2));
% mindpc(:,3) = fdr(mindp(:,3));
%% Time chunk pDMN search for mindfulness
clc
MH = mind;
data =squeeze(source_all_high(a,:,1:1000,:));
pDMNData = squeeze(nanmean(data(netwrk(2).roi(6),:,:),1));

count=0;
winsize = 25;
ps=[];
for t=1:winsize:1000-winsize
    count=count+1;
    [rs(count),ps(count)] = corr(MH,nanmean(pDMNData(t:t+winsize,:),1)','Type','Spearman','Rows','pairwise');
end

find(fdr(ps')<0.05)



%% ANOVA with depression
a=2;
tbl = table;
tbl.dep3=1*(dep>4) + (dep>9);
tbl.dep2 = 1*(dep>4);
tbl.dep3(isnan(dep))=nan;
tbl.dep2(isnan(dep))=nan;
tbl.mind2 = 1*(mind>nanmedian(mind));
tbl.mind2(isnan(mind))=nan;

% tbl.atten=f(a).high_pre;
% tbl.dist=f(a).low_pre;
% tmp5h = filloutliers(tmp5h, nan);
% tmp5l = filloutliers(tmp5l, nan);
% tbl.df=tmp5l-tmp5h;

% tbl
% anova1(f(a).high_pre,tbl.dep3)
% Number of groups
numGroups = 3; % Assuming groups are 0, 1, and 2

% Initialize a cell array to store the results
% anovaResults = cell(1, 68);
data =  f(2).df_pre;
% data=yp';
data = filloutliers(data,nan,2);

% data = f(a).nethi';
groups = 1*tbl.mind2; 
% groups(age>60)=nan;
anovaResults=[];

for i = 1:size(data,1)
    % Extract the data for the current variable
    variableData = data(i,:);
    
    % Run the one-way ANOVA
    [p, at, stats] = anova1(variableData, groups, 'off');
    [~,~,stats]=glmfit(variableData, mind);

    % Store the results for the current variable
    anovaResults(i) = p;
    anovaTbl{i}=at;
end

clc
atlas.label(find((anovaResults)<0.05))
for n=1:6
    anovaResults(netwrk(n).roi) = fdr(anovaResults(netwrk(n).roi));
end
atlas.label(find((anovaResults)<0.05))

%% 
data = f(a).low_pre;
tbl.lpc = data(51,:)';
tbl.dep=dep;
tbl.mind=mind;

fitglm(tbl, 'dep3 ~ 1+ lpc + mind')

% [b,dev,stats]=glmfit(data(51,:)',dep,'poisson')
%% Precuneus anova
groups = 1*tbl.mind2; 
singleData  = f(2).df_pre(51,:);
% singleData = filloutliers(singleData,nan);

close all
[p, at, stats] = anova1(singleData, groups, 'on');
xticklabels({'Low Mindfulness','High Mindfulness'})
set(gca,'fontweight','bold','fontsize',12)
title('rPC Attended','fontweight','bold')

%% PLotting Anova as box plots
groups = 1*tbl.mind2; 
data =  f(2).high_pre(51,:);

goodSubs = ~isnan(groups) & ~isnan(data');
groups = groups(goodSubs);
data = data(goodSubs);
catData = [data(groups==0) data(groups==1)];
% close all

% subplot(4,4,[11,15]) % with line plots
subplot(2,2,4)
hold on
netOrder = ones(size(catData));
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
boxchart(netOrder', catData,'GroupByColor',[ones(1,sum(groups==0)), 2*ones(1,sum(groups))],'MarkerStyle','none')

% data=round(data,3);
x = repmat(1,(sum(groups==0)),1)-0.25;
x=reshape(x,[],1);
swarmchart(x,[data(groups==0)],50,'.','XJitterWidth',0.4,'MarkerEdgeColor',[0 0.4470 0.7410]);
x = repmat(1,(sum(groups)),1)+0.25;
x=reshape(x,[],1);
swarmchart(x,[data(groups==1)],50,'.','XJitterWidth',0.4,'MarkerEdgeColor',[0.8500 0.3250 0.0980]);

ylabel('lPC')
title('high Cons')
% legend({'Low Mindfulness','High Mindfulness'},'box','off')
% Precuneus over time
% clc
% MH = mind;
% tbl=table();
% tbl.mind2 = 1*(mind>nanmedian(mind));
% tbl.mind2(isnan(mind))=nan;
% group=1*tbl.mind2;
% lPC = squeeze(source_all_low(2,51,1:1000,:)-source_all_high(2,51,1:1000,:));
% 
% count=0;
% winsize = 10;
% p=[];betas=[];
% for t=1:winsize:1000-winsize
%     count=count+1;
%     variableData = nanmean(lPC(t:t+winsize,:),1)';
%     % Run the one-way ANOVA
% %     [p(count), at, stats] = anova1(variableData, groups, 'off');
% 
%     [b,~,stats]=glmfit(variableData, tbl.mind2);
%     p(count)=stats.p(2);
%     betas(count)=b(2);
%     % Store the results for the current variable
% %     anovaTbl{i}=at;
% end
% % close all
% subplot(4,4,[12,16])
% sigIdx = find(fdr(p)<0.05);
% tvect = linspace(-4000,0,length(betas)); 
% plot(tvect, betas,'linewidth',1.5)
% hold on
% plot(tvect(sigIdx), ones(size(sigIdx)), '*')
% % title('Distracted Trials','fontweight','bold')
% ylabel('glm Beta','fontweight','bold')
% xlabel('Time (sec)','fontweight','bold')

%% Overlap between mindfulness and depression 

close all;
a = 2;
MH = mind;
% MH(age>60)=nan;

tmp5h=f(a).high_pre;
tmp5l=f(a).low_pre;
% tmp5h = filloutliers(tmp5h, nan);
% tmp5l = filloutliers(tmp5l, nan);
tmp5df=tmp5l-tmp5h;

%%% Mindfulness Corrs
[mindr1,mindp1]=corr(MH,tmp5l','Type','Spearman','Rows','pairwise');
[mindr2,mindp2]=corr(MH,tmp5h','Type','Spearman','Rows','pairwise');
[mindr3,mindp3]=corr(MH,tmp5df','Type','Spearman','Rows','pairwise');
mindr=[mindr1' mindr2', mindr3'];
mindp=[mindp1' mindp2', mindp3'];
mindpc=mindp;

tbl = table;
tbl.dep3=(dep>5) + (dep>10);
tbl.dep2 = (dep>10);
groups = tbl.dep3;
anovaResults=[];

% low
for i = 1:68
    [p, at, stats] = anova1(tmp5l(i,:), groups, 'off');
    anovaResults(i,1) = p;
end

% high
for i = 1:68
    [p, at, stats] = anova1(tmp5h(i,:), groups, 'off');
    anovaResults(i,2) = p;
end

% diff
for i = 1:68
    [p, at, stats] = anova1(tmp5df(i,:), groups, 'off');
    anovaResults(i,3) = p;
end
clc
bothSigs = (mindpc <0.05 & anovaResults<0.05);
atlas.label(bothSigs(:,1))
atlas.label(bothSigs(:,2))
atlas.label(bothSigs(:,3))

%%
hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'dist','atten','diff'};
var = zeros(68,1);
var(51)=-1;
plot68roi(hm, T'*var,[],labelnames);
set(gcf,'Position',[965.6667 352.3333 820.6667 1.0087e+03]);
clc

%% plot scatter 
x1=mind;
y1=tmp5l(51,:);

% Create figure
figure1 = figure('Position',[1200 900 180 125]);
mdl = fitlm(x1,y1,'linear','RobustOpts','on')
plot1=plot(mdl,'MarkerSize',7,'Marker','.','MarkerEdgeColor',[0    0.4470    0.7410]);

set(plot1(2),'LineStyle','-','LineWidth',1);
set(plot1(3),'LineStyle',':','LineWidth',2);
set(plot1(4),'LineStyle',':','LineWidth',2);
box off
legend('off')
title ('')
xlabel('mindfulness');
ylabel('L Precunes activity');
%% Fig 4B Correlation with mindfulness and inc act of lo con trials.
t=-3996:4:3996;
f=[];
for a = 1:3
    tw_pre=[find(t>-4000 & t<4000)];
    tw_pre=tw_pre(1,[1,end]);
    f(a).high_pre=squeeze(nanmean(source_all_high(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).low_pre=squeeze(nanmean(source_all_low(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).df_pre=(f(a).low_pre-f(a).high_pre);

    f(a).high_prets=squeeze((source_all_high(a,:,tw_pre(1):tw_pre(2),:)));
    f(a).low_prets=squeeze((source_all_low(a,:,tw_pre(1):tw_pre(2),:)));

    [~,p]=ttest(f(a).df_pre');
    f(a).dfp=p;
    f(a).dfpc=fdr(p);
    
    [h,p]=ttest(f(a).high_pre');
    f(a).highp=p;
    f(a).highpc=fdr(p);
    
    [h,p]=ttest(f(a).low_pre');
    f(a).lowp=p;
    f(a).lowpc=fdr(p);
    
    for b = 1:5
        f(a).netlo(:,b)=nanmean(f(a).low_pre(netwrk(b).roi,:),1)';
        f(a).nethi(:,b)=nanmean(f(a).high_pre(netwrk(b).roi,:),1)';
        f(a).netdf(:,b)=nanmean(f(a).df_pre(netwrk(b).roi,:),1)';
        
        f(a).netlots(:,:,b)=squeeze(nanmean(f(a).low_prets(netwrk(b).roi,:,:),1));
        f(a).nethits(:,:,b)=squeeze(nanmean(f(a).high_prets(netwrk(b).roi,:,:),1));
        
    end
end



for a=2
    for j = 1:5
    %netlots2(j,:,:)=(squeeze((zscore(f(a).netlots(:,:,j),0,1))));
   % nethits2(j,:,:)=(squeeze((zscore(f(a).nethits(:,:,j),0,1))));
    netlots2(j,:,:)=(squeeze(f(a).netlots(:,:,j)));
   nethits2(j,:,:)=(squeeze(f(a).nethits(:,:,j)));
    
%    netlotsem(a,j,:)=(squeeze(nanstd(f(a).netlots(:,:,j),0,2)))/sqrt(ns);
 %   nethitsem(a,j,:)=(squeeze(nanstd(f(a).nethits(:,:,j),0,2)))/sqrt(ns);

    end
end
%%
t2=t(tw_pre(1):tw_pre(2));
close all
figure1 = figure('Position',[1200 900 350 120]);
networknames={'FPN','CON','aDMN','pDMN','mtDMN'};
for i = 1:5
    tmp=squeeze(netlots2(i,:,:)-nethits2(i,:,:));
    tmp=squeeze(nethits2(i,:,:));
    tmp=squeeze(netlots2(i,:,:));
    subplot(1,5,i)
    [mndcorr,p]=corr(mind,tmp','Type','Spearman','Rows','pairwise');
    plot(t2,mndcorr);     hold on; 
    box off
%    plot([-4000;4000],[.12;.12],'r--');
 %   plot([-4000;4000],[-.12;-.12],'r--');
    axis([-4000 0 -0.3 .3]);
    title(networknames{i});
    [~,~,~,pc]=fdr_bh(p(1:end));
    sig=find(pc<0.05);
    plot(t2(sig),ones(length(sig),1)*-.25,'.r');
%close all;
end  
    
