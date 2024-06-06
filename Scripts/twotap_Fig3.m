% Load relevant variables

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
% age(outl)=[];

size(source_all_low)
size(source_all_high)
% size(age)

%% For network-based analyses
netwrk(1).name = 'FPN'
netwrk(1).roi = [5,6,55,56,59,60];
netwrk(2).name = 'CON'
netwrk(2).roi = [3,4,19,20,37,38,39,40,41,42,57,58,67,68];
netwrk(3).name = 'aDMN'
netwrk(3).roi = [11,12,25,26,29,30,53,54];
netwrk(4).name = 'pDMN'
netwrk(4).roi = [15,16,21,22,51,52];
netwrk(5).name = 'MTL-DMN'
netwrk(5).roi = [9,10,17,18,31,32,35,36,61,62,65,66];

%% organizing data into temporal epochs for analysis.

t=-3996:4:3996;
f=[];
for a = 1:3
    tw_pre=[find(t>-2000 & t<000)];
    tw_pre=tw_pre(1,[1,end]);
    f(a).high_pre=squeeze(nanmean(source_all_high(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).low_pre=squeeze(nanmean(source_all_low(a,:,tw_pre(1):tw_pre(2),:),3));
    f(a).df_pre=(f(a).low_pre-f(a).high_pre);
    
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
    end
end

% fdr_correction of all brain diff
tmpp1=[f(1).dfp,f(2).dfp,f(3).dfp];
p_all1=fdr(tmpp1);
f(1).dfpc2=p_all1(1:68);
f(2).dfpc2=p_all1(69:136);
f(3).dfpc2=p_all1(137:end);

tmpp2=[f(1).highp,f(2).highp,f(3).highp];
p_all2=fdr(tmpp2);
f(1).highpc2=p_all2(1:68);
f(2).highpc2=p_all2(69:136);
f(3).highpc2=p_all2(137:end);

tmpp3=[f(1).lowp,f(2).lowp,f(3).lowp];
p_all3=fdr(tmpp3);
f(1).lowpc2=p_all3(1:68);
f(2).lowpc2=p_all3(69:136);
f(3).lowpc2=p_all3(137:end);
close all;
clc;

% %% Figure 3A whole brain contrasts at each frequency.
% close all
% for i = 2%1:3
%     hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
%     T = hm.indices4Structure(hm.atlas.label);
%     T = double(T)';
%     P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%     labelnames = {''};
%     var = nanmean(f(i).df_pre');
%     var(f(i).dfpc2>0.05)=0;
%     plotMultipanel_v2_jet_indBrainMaps_modifcorr(hm, T'*var',[-0.0005 .0005],labelnames);
%     set(gcf,'Position',[1,1,267,476]);
% end


%% Figure 3A Images of networks.

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'FPN'};
var=zeros(1,68);
var(netwrk(1).roi)=1;
plotMultipanel_v2_jet_indBrainMaps_modifcorr_nolegend(hm, T'*var',0,labelnames);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'CON'};
var=zeros(1,68);
var(netwrk(2).roi)=1;
plotMultipanel_v2_jet_indBrainMaps_modifcorr_nolegend(hm, T'*var',0,labelnames);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'aDMN'};
var=zeros(1,68);
var(netwrk(3).roi)=1;
plotMultipanel_v2_jet_indBrainMaps_modifcorr_nolegend(hm, T'*var',0,labelnames);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'pDMN'};
var=zeros(1,68);
var(netwrk(4).roi)=1;
plotMultipanel_v2_jet_indBrainMaps_modifcorr_nolegend(hm, T'*var',0,labelnames);

hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
%set(0,'DefaultFigureWindowStyle','docked');
labelnames = {'mtDMN'};
var=zeros(1,68);
var(netwrk(5).roi)=1;
plotMultipanel_v2_jet_indBrainMaps_modifcorr_nolegend(hm, T'*var',0,labelnames);

%% Fig 3 data (fig / stats prism).
% freqs={'theta','alpha','beta'};
for i = 1:3    
%         Fig3Adata(i).raw=(f(i).df_pre')';
%         Fig3Adata(i).mean=nanmean(f(i).df_pre')';
%         Fig3Adata(i).sem=nanstd(f(i).df_pre')/sqrt(size((f(i).df_pre'),1))';
%         Fig3Adata(i).p=f(1).dfpc2';
%        Fig2Adata(i).freq=freqs{i};

        Fig3ABdata(i).lo=(f(i).netlo);
        Fig3ABdata(i).hi=(f(i).nethi);
        Fig3ABdata(i).diff=(f(i).netdf);        
        %Fig3Bdata(i).mean=nanmean(f(i).netdf,1);
        %Fig3Bdata(i).sem=nanstd(f(i).netdf)'/sqrt(size(f(i).netdf,1));        
 %       Fig2bBdata(i).freq=freqs{i};
end
 sdir=['/Users/Dhak/Documents/2Tapfinaldata/June2023/figdata/Fig3a_bdata.mat'];
 save (sdir,'Fig3ABdata','-v7.3');

