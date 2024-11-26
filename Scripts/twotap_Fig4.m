% Load relevant variables

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

size(source_all_low)
size(source_all_high)
size(rest)

%%
behdata=readtable("A:\Resting\Total_Data_324_for_Correlation.xlsx");
%%
dep=behdata.Depression;
anx=behdata.Anxiety;
mind=behdata.Mindfulness;
% beh=table2array(behdata(:,12));
statemind=behdata.state_mindfulness;
age=behdata.Age;

dep(outl)=[];
anx(outl)=[];
mind(outl)=[];
% beh(outl)=[];
statemind(outl)=[];
% yall=[anx,dep,mind,beh];
age(outl)=[];
size(age)

%%

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

%% organizing resting data into temporal epochs for analysis.

trest=-3996:4:0;

for a = 1:3%2
    tw_prer=[find(trest>-4000 & trest<4000)];
    tw_prer=tw_prer(1,[1,end]);
    f(a).r=squeeze(nanmean(rest(a,:,tw_prer(1):tw_prer(2),:),3));
    f(a).rts=squeeze((source_all_low(a,:,tw_prer(1):tw_prer(2),:)));

    [~,p]=ttest(f(a).r');
    f(a).rp=p;
    f(a).rpc=fdr(p);
          
    for b = 1:5
        f(a).netr(:,b)=nanmean(f(a).r(netwrk(b).roi,:,:),1)';
        f(a).netrts(:,:,b)=squeeze(nanmean(f(a).rts(netwrk(b).roi,:,:),1));

    end
end 

%% Freq x Freq correlations
win=250;
%[H,P,ADSTAT,CV] = adtest(f(2).netlots(:))


clear FPNnet* CONnet* aDMNnet* pDMNnet* mtDMNnet*
ns=size(f(2).netlots,2);

for f1=1:3
    %for f2=1
        for roi2=1:5%
            for j = 1:size(f(f1).netlots,2)% N subjects
                cnt=0;
                cnt1=0;
                for i = 1:25:1998-win
                    cnt=cnt+1;
%                     FPNnetlocorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,1)),zscore(f(f1).netlots(i:win+i,j,net2)),'Type','Pearson'));
%                     FPNnethicorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,1)),zscore(f(f1).nethits(i:win+i,j,net2)),'Type','Pearson'));
%                     CONnetlocorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,2)),zscore(f(f1).netlots(i:win+i,j,net2)),'Type','Pearson'));
%                     CONnethicorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,2)),zscore(f(f1).nethits(i:win+i,j,net2)),'Type','Pearson'));
%                     
%                     aDMNnetlocorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,3)),zscore(f(f1).netlots(i:win+i,j,net2)),'Type','Pearson'));
%                     aDMNnethicorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,3)),zscore(f(f1).nethits(i:win+i,j,net2)),'Type','Pearson'));
%                     pDMNnetlocorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,4)),zscore(f(f1).netlots(i:win+i,j,net2)),'Type','Pearson'));
%                     pDMNnethicorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,4)),zscore(f(f1).nethits(i:win+i,j,net2)),'Type','Pearson'));
%                     mtDMNnetlocorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,5)),zscore(f(f1).netlots(i:win+i,j,net2)),'Type','Pearson'));
%                     mtDMNnethicorr(f1,net2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,5)),zscore(f(f1).nethits(i:win+i,j,net2)),'Type','Pearson'));
% %                     
                    FPNnetlocorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,1)),zscore(f(f1).netlots(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    FPNnethicorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,1)),zscore(f(f1).nethits(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    CONnetlocorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,2)),zscore(f(f1).netlots(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    CONnethicorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,2)),zscore(f(f1).nethits(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    
                    aDMNnetlocorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,3)),zscore(f(f1).netlots(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    aDMNnethicorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,3)),zscore(f(f1).nethits(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    pDMNnetlocorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,4)),zscore(f(f1).netlots(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    pDMNnethicorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,4)),zscore(f(f1).nethits(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    mtDMNnetlocorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).netlots(i:win+i,j,5)),zscore(f(f1).netlots(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    mtDMNnethicorr(f1,roi2,j,cnt)=atanh(corr(zscore(f(f1).nethits(i:win+i,j,5)),zscore(f(f1).nethits(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                
                end
                for i = 1:25:1000-win
                    cnt1=cnt1+1;
                    FPNnetlocorr_r(f1,roi2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,1)),zscore(f(f1).netrts(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    CONnetlocorr_r(f1,roi2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,2)),zscore(f(f1).netrts(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));                    
                    aDMNnetlocorr_r(f1,roi2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,3)),zscore(f(f1).netrts(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    pDMNnetlocorr_r(f1,roi2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,4)),zscore(f(f1).netrts(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                    mtDMNnetlocorr_r(f1,roi2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,5)),zscore(f(f1).netrts(i:win+i,j,roi2)),'rows','pairwise','Type','Spearman'));
                  
%                     FPNnetlocorr_r(f1,net2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,1)),zscore(f(f1).netrts(i:win+i,j,net2)),'rows','pairwise','Type','Pearson'));
%                     CONnetlocorr_r(f1,net2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,2)),zscore(f(f1).netrts(i:win+i,j,net2)),'rows','pairwise','Type','Pearson'));                    
%                     aDMNnetlocorr_r(f1,net2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,3)),zscore(f(f1).netrts(i:win+i,j,net2)),'rows','pairwise','Type','Pearson'));
%                     pDMNnetlocorr_r(f1,net2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,4)),zscore(f(f1).netrts(i:win+i,j,net2)),'rows','pairwise','Type','Pearson'));
%                     mtDMNnetlocorr_r(f1,net2,j,cnt1)=atanh(corr(zscore(f(f1).netrts(i:win+i,j,5)),zscore(f(f1).netrts(i:win+i,j,net2)),'rows','pairwise','Type','Pearson'));
                end
            end
        %end
    end
end
%% Save activity network conn data
f1=2;
pDMNFPN_loConn = [];pDMNFPN_hiConn=[];pDMNFPN_lohiConn=[];
for r1=1:length(netwrk(4).roi)
    roi1=netwrk(4).roi(r1);
    for r2 =1:length(netwrk(4).roi)
        roi2=netwrk(1).roi(r2);
        for j = 1:324 % # subjects
            cnt=0;
            cnt1=0;
            for i = 1:25:1000-win
                cnt=cnt+1;
%                 pDMNFPN_loConn(r1,r2,j,cnt)=atanh(corr(zscore(squeeze(source_all_low(f1,roi1,i:win+i,j))),zscore(squeeze(source_all_low(f1,roi2,i:win+i,j))),'rows','pairwise','Type','Spearman'));
%                 pDMNFPN_hiConn(r1,r2,j,cnt)=atanh(corr(zscore(squeeze(source_all_high(f1,roi1,i:win+i,j))),zscore(squeeze(source_all_high(f1,roi2,i:win+i,j))),'rows','pairwise','Type','Spearman'));
                pDMNFPN_lohiConn(r1,r2,j,cnt)=atanh(corr(zscore(squeeze(source_all_low(f1,roi1,i:win+i,j))),zscore(squeeze(source_all_high(f1,roi2,i:win+i,j))),'rows','pairwise','Type','Spearman'));

%                 pROIFPN_hiConn(r1,1,j,cnt) = atanh(corr(zscore(squeeze(source_all_high(f1,roi1,i:win+i,j))),zscore(f(f1).nethits(i:win+i,j,1)),'rows','pairwise','Type','Spearman'));
%                 pROIFPN_loConn(r1,1,j,cnt) = atanh(corr(zscore(squeeze(source_all_low(f1,roi1,i:win+i,j))),zscore(f(f1).netlots(i:win+i,j,1)),'rows','pairwise','Type','Spearman'));

            end
        end
    end
end

%% roi x roi connectivity
MH=mind;
tbl.mind2 = 1*(mind>nanmedian(mind));
tbl.mind2(isnan(mind))=nan;
group=tbl.mind2;
% MH(age>60)=nan;

% PFlow = squeeze(nanmean(pDMNFPN_loConn,4));
% PFhigh = squeeze(nanmean(pDMNFPN_hiConn,4));
PFlohigh = squeeze(nanmean(pDMNFPN_lohiConn,4));

mp=[];
clc
for r1=1:6
    for r2=1:6
%         [mr(r1,r2),mp(r1,r2)]=corr(MH,squeeze(PFlohigh(r1,r2,:)),'Type','Spearman','Rows','pairwise');
        [mp(r1,r2), at, stats] = anova1(squeeze(PFlohigh(r1,r2,:)), groups, 'off');
    end
end
mp
%% single roi to all networks 
MH=mind;
MH(age>60)=nan;
PFlow = squeeze(nanmean(pROIFPN_loConn,4));
PFhigh = squeeze(nanmean(pROIFPN_hiConn,4));
mp=[];mr=[];
for r2=1:6
    [mr(r2),mp(r2)]=corr(MH,squeeze(PFhigh(r2,:))','Type','Spearman','Rows','pairwise');
    [mp(r2), at, stats] = anova1(squeeze(PFlow(r2,:))', groups, 'off');

end
mp



%% Al Freq Figure (3A)
Netnames={'FPN','CON','aDMN','pDMN','mtDMN'};
t3=[1:25:1998-win]+round(win/2);
t3=t(t3)/1000;

t3r=[1:25:1000-win]+round(win/2);
t3r=trest(t3r)/1000;

tmp1=squeeze(nanmean(FPNnetlocorr(2,:,:,:),3));
tmp2=squeeze(nanmean(FPNnethicorr(2,:,:,:),3));
tmp3=squeeze(nanmean(FPNnetlocorr_r(2,:,:,:),3));

tmp1std=squeeze(nanstd(FPNnetlocorr(2,:,:,:),[],3))/sqrt(ns);
tmp2std=squeeze(nanstd(FPNnethicorr(2,:,:,:),[],3))/sqrt(ns);
tmp3std=squeeze(nanstd(FPNnetlocorr_r(2,:,:,:),[],3))/sqrt(ns);
figure;

for i = 2:5
    subplot(5,5,i)
    shadedErrorBar(t3,tmp1(i,:),tmp1std(i,:),'lineProps','b'); hold on;
    shadedErrorBar(t3,tmp2(i,:),tmp2std(i,:),'lineProps','r')
%    shadedErrorBar(t3r,tmp3(i,:),tmp3std(i,:),'lineProps','g')
    title([Netnames{i}]);
    axis([-2.500 0 -inf inf]);
end

tmp1=squeeze(nanmean(CONnetlocorr(2,:,:,:),3));
tmp2=squeeze(nanmean(CONnethicorr(2,:,:,:),3));

tmp1std=squeeze(nanstd(CONnetlocorr(2,:,:,:),[],3))/sqrt(ns);
tmp2std=squeeze(nanstd(CONnethicorr(2,:,:,:),[],3))/sqrt(ns);

for i = 3:5
    subplot(5,5,i+5)
    shadedErrorBar(t3,tmp1(i,:),tmp2std(i,:),'lineProps','b'); hold on;
    shadedErrorBar(t3,tmp2(i,:),tmp2std(i,:),'lineProps','r')
    title([Netnames{i}]);
    axis([-2.500 1 -inf inf]);
end


tmp1=squeeze(nanmean(aDMNnetlocorr(2,:,:,:),3));
tmp2=squeeze(nanmean(aDMNnethicorr(2,:,:,:),3));
tmp1std=squeeze(nanstd(aDMNnetlocorr(2,:,:,:),[],3))/sqrt(ns);
tmp2std=squeeze(nanstd(aDMNnethicorr(2,:,:,:),[],3))/sqrt(ns);

for i = 4:5
    subplot(5,5,i+10)
    shadedErrorBar(t3,tmp1(i,:),tmp2std(i,:),'lineProps','b'); hold on;
    shadedErrorBar(t3,tmp2(i,:),tmp2std(i,:),'lineProps','r')
    title([Netnames{i}]);
    axis([-2.500 1 -inf inf]);
end

tmp1=squeeze(nanmean(pDMNnetlocorr(2,:,:,:),3));
tmp2=squeeze(nanmean(pDMNnethicorr(2,:,:,:),3));
tmp1std=squeeze(nanstd(pDMNnetlocorr(2,:,:,:),[],3))/sqrt(ns);
tmp2std=squeeze(nanstd(pDMNnethicorr(2,:,:,:),[],3))/sqrt(ns);

for i = 5:5
    subplot(5,5,i+15)
    shadedErrorBar(t3,tmp1(i,:),tmp2std(i,:),'lineProps','b'); hold on;
    shadedErrorBar(t3,tmp2(i,:),tmp2std(i,:),'lineProps','r')    
    title([Netnames{i}]);
    axis([-2.500 1 -inf inf]);
end


tmp1=squeeze(nanmean(mtDMNnetlocorr(2,:,:,:),3));
tmp2=squeeze(nanmean(mtDMNnethicorr(2,:,:,:),3));
tmp1std=squeeze(nanstd(mtDMNnetlocorr(2,:,:,:),[],3))/sqrt(ns);
tmp2std=squeeze(nanstd(mtDMNnethicorr(2,:,:,:),[],3))/sqrt(ns);


%% analysis within tw of interest
tmptw=find(t3>-4 & t3<0);

FPNnetlocorrm=reshape(permute(squeeze(nanmean(FPNnetlocorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
CONnetlocorrm=reshape(permute(squeeze(nanmean(CONnetlocorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
aDMNnetlocorrm=reshape(permute(squeeze(nanmean(aDMNnetlocorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
pDMNnetlocorrm=reshape(permute(squeeze(nanmean(pDMNnetlocorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
mtDMNnetlocorrm=reshape(permute(squeeze(nanmean(mtDMNnetlocorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';

FPNnethicorrm=reshape(permute(squeeze(nanmean(FPNnethicorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
CONnethicorrm=reshape(permute(squeeze(nanmean(CONnethicorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
aDMNnethicorrm=reshape(permute(squeeze(nanmean(aDMNnethicorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
pDMNnethicorrm=reshape(permute(squeeze(nanmean(pDMNnethicorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';
mtDMNnethicorrm=reshape(permute(squeeze(nanmean(mtDMNnethicorr(:,:,:,tmptw),4)),[2,1,3]),[15,ns])';

FPNnetrestm=reshape(permute(squeeze(nanmean(FPNnetlocorr_r(:,:,:,5:25),4)),[2,1,3]),[15,ns])';
CONnetrestm=reshape(permute(squeeze(nanmean(CONnetlocorr_r(:,:,:,5:25),4)),[2,1,3]),[15,ns])';
aDMNnetrestm=reshape(permute(squeeze(nanmean(aDMNnetlocorr_r(:,:,:,5:25),4)),[2,1,3]),[15,ns])';
pDMNnetrestm=reshape(permute(squeeze(nanmean(pDMNnetlocorr_r(:,:,:,5:25),4)),[2,1,3]),[15,ns])';
mtDMNnetlrestm=reshape(permute(squeeze(nanmean(mtDMNnetlocorr_r(:,:,:,5:25),4)),[2,1,3]),[15,ns])';


% FPNnetdiff=FPNnetlocorrm-FPNnethicorrm;
% CONnetdiff=CONnetlocorrm-CONnethicorrm;
% aDMNnetdiff=aDMNnetlocorrm-aDMNnethicorrm;
% pDMNnetdiff=pDMNnetlocorrm-pDMNnethicorrm;
% mtDMNnetdiff=mtDMNnetlocorrm-mtDMNnethicorrm;

%% for SPSS rmANOVA arrangement (data for Fig 3B-C)

%Theta_corr_spsslo=[FPNnetlocorrm(:,2:5),CONnetlocorrm(:,3:5),aDMNnetlocorrm(:,4:5),pDMNnetlocorrm(:,5)];
%Theta_corr_spsshi=[FPNnethicorrm(:,2:5),CONnethicorrm(:,3:5),aDMNnethicorrm(:,4:5),pDMNnethicorrm(:,5)];
%Theta_corr_diff=Theta_corr_spsslo-Theta_corr_spsshi;

Alpha_corr_spsslo=[FPNnetlocorrm(:,7:10),CONnetlocorrm(:,8:10),aDMNnetlocorrm(:,9:10),pDMNnetlocorrm(:,10)];
Alpha_corr_spsshi=[FPNnethicorrm(:,7:10),CONnethicorrm(:,8:10),aDMNnethicorrm(:,9:10),pDMNnethicorrm(:,10)];
Alpha_corr_rest=[FPNnetrestm(:,7:10),CONnetrestm(:,8:10),aDMNnetrestm(:,9:10),pDMNnetrestm(:,10)];

%Alpha_corr_diff=Alpha_corr_spsslo-Alpha_corr_spsshi;

%Beta_corr_spsslo=[FPNnetlocorrm(:,12:15),CONnetlocorrm(:,13:15),aDMNnetlocorrm(:,14:15),pDMNnetlocorrm(:,15)];
%Beta_corr_spsshi=[FPNnethicorrm(:,12:15),CONnethicorrm(:,13:15),aDMNnethicorrm(:,14:15),pDMNnethicorrm(:,15)];
%Beta_corr_diff=Beta_corr_spsslo-Beta_corr_spsshi;

%%
t=-3996:4:3996;
f=[];
for a = 1:3
    tw_pre=[find(t>-4000 & t<000)];
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

%% Fig 4 saving data
% aDMN
savePath = 'A:\TwoTap\Manuscript\Figure';
sfile=[savePath '\tempData.xlsx'];

ya = [aDMNnetlocorrm(:,6:7) aDMNnethicorrm(:,6:7) aDMNnetrestm(:,6:7)]; % FPN-CON corrs
T=array2table(ya, 'VariableNames',{'dist_FPN','dist_CON','atten_FPN','atten_CON','rest_FPN','rest_CON'});
% writetable(T, sfile, 'sheet','aDMN_netConn')
% mtlDMN
ym = [mtDMNnetlocorrm(:,6:7) mtDMNnethicorrm(:,6:7) mtDMNnetlrestm(:,6:7)]; % FPN-CON corrs
T=array2table(ym, 'VariableNames',{'dist_FPN','dist_CON','atten_FPN','atten_CON','rest_FPN','rest_CON'});
% writetable(T, sfile, 'sheet','mtlDMN_netConn')
% pDMN
yp = [pDMNnetlocorrm(:,6:7) pDMNnethicorrm(:,6:7) pDMNnetrestm(:,6:7)]; % FPN-CON corrs
T=array2table(yp, 'VariableNames',{'dist_FPN','dist_CON','atten_FPN','atten_CON','rest_FPN','rest_CON'});
% writetable(T, sfile, 'sheet','pDMN_netConn')

%% Fig 4 Grouped bar charts

close all
figure

% subplot(131)
% data = reshape(nanmean(ya,1),[2,3]);
% sem = reshape(nanstd(ya,[],1),[2,3])/sqrt(324);
% plotErrBar(data',sem',{'white','k'})
% set(gca,'xticklabels',{'Dist', 'Atten','Rest'},'fontweight','bold','fontsize',16,'xticklabelrotation',90)
% title('aDMN')
% ylim([0 0.6])
% % legend('Distracted Trials','Attended Trials','box','off','location','northwest')
% % set(gca,'XGrid','off','YGrid','on','Position',[0.1300 0.1100 0.7750 0.8121])
% 
% 
% subplot(132)
% data = reshape(nanmean(ym,1),[2,3]);
% sem = reshape(nanstd(ym,[],1),[2,3])/sqrt(324);
% plotErrBar(data',sem',{'white','k'})
% set(gca,'xticklabels',{'Dist', 'Atten','Rest'},'fontweight','bold','fontsize',16,'xticklabelrotation',90)
% title('mtDMN')
% ylim([0 0.6])


% subplot(133)
subplot(131)
hold on
data = reshape(nanmean(yp,1),[2,3]);
sem = reshape(nanstd(yp,[],1),[2,3])/sqrt(324);
% plotErrBar(data',sem',{'white','k'})
allFPN = reshape(yp(:,[1,3,5]),[],1);
allCON = reshape(yp(:,[2,4,6]),[],1);
netOrder =(repmat(reshape(repmat([1:3],[324],1),[],1),2,1));
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
bbb=boxchart(netOrder, [allFPN;allCON],'GroupByColor',[ones(1,324*3), 2*ones(1,324*3)]','MarkerStyle','none')

x = repmat(1:3,324,1)-0.25;
x=reshape(x,[],1);
swarmchart(x,[allFPN],20,'.','XJitterWidth',0.1,'MarkerEdgeColor',[0 0.4470 0.7410]);
x = repmat(1:3,324,1)+0.25;
x=reshape(x,[],1);
swarmchart(x,[allCON],20,'.','XJitterWidth',0.1,'MarkerEdgeColor',[0.8500 0.3250 0.0980]);
% ylim([-0.0001 0.018])
legend('Low Cons','High Cons','box','off','location','northwest','box','off')
set(gca,'xticklabels',{'Dist', 'Atten','Rest'},'xticklabelrotation',90)
% title('pDMN')
% ylim([0 0.6])
legend('FPN','CON')
set(gcf, 'position', [821.6667 567.6667 766.0000 298])
% saveas(gcf,[savePath '\Fig4a.tiff']);
%% Fig 4D Activity Conn correlation
%panel 1

% y = aDMNnetlocorrm(:,6:7);
% y1= aDMNnetlocorrm(:,7); % alpha CON
% x1 = f(2).netlo(:,3); % aDMN 
% y1=filloutliers(y1,nan);
% x1=filloutliers(x1,nan);
% [Act_conn_corr1,Act_conn_corrp1]=corr(x1,y1,'Type','Spearman','Rows','pairwise')
% close all

% Create figure
% subplot(131)
% % figure1 = figure('Position',[1200 900 180 125]);
% mdl = fitlm(y1,x1,'linear','RobustOpts','on')
% plot1=plot(mdl,'MarkerSize',7,'Marker','.','MarkerEdgeColor',[0    0.4470    0.7410]);
% 
% set(plot1(2),'LineStyle','-','LineWidth',1);
% set(plot1(3),'LineStyle',':','LineWidth',2);
% set(plot1(4),'LineStyle',':','LineWidth',2);
% % ylim([-0.02 0.02])
% box off
% legend('off')
% title ('')
% xlabel('aDMN-CON conn');
% ylabel('aDMN activity');
% title("r=" + string(round(Act_conn_corr1,2)) + "; p=" + string(round(Act_conn_corrp1,3)))

% panel 2 mtDMN
% % y = mtDMNnetlocorrm(:,6:7);
% y2 = mtDMNnetlocorrm(:,7);
% x2 = f(2).netlo(:,5);
% y2=filloutliers(y2,nan);
% x2=filloutliers(x2,nan);
% [Act_conn_corr3,Act_conn_corrp3]=corr(x2,y2,'Type','Spearman','Rows','pairwise');
% %[~,~,~,Act_conn_corrp3]=fdr_bh(Act_conn_corrp3)
% % figure1 = figure('Position',[1200 900 180 125]);
% subplot(132)
% mdl2 = fitlm(y2,x2,'linear','RobustOpts','on')
% plot1=plot(mdl2,'MarkerSize',7,'Marker','.','MarkerEdgeColor',[0    0.4470    0.7410]);
% 
% set(plot1(2),'LineStyle','-','LineWidth',1);
% set(plot1(3),'LineStyle',':','LineWidth',2);
% set(plot1(4),'LineStyle',':','LineWidth',2);
% box off
% legend('off')
% title ('')
% ylabel('mtDMN activity');
% xlabel('mtDMN-CON conn');
% title("r=" + string(round(Act_conn_corr3,2)) + "; p=" + string(round(Act_conn_corrp3,3)))

% panel 3 pDMN

% y = pDMNnetlocorrm(:,6:7);
x3= pDMNnetlocorrm(:,7);
y3 = f(2).netlo(:,4);
% y3=filloutliers(y3,nan);
% x3=filloutliers(x3,nan);
[Act_conn_corr2,Act_conn_corrp2m]=corr(x3,y3,'Type','Spearman','Rows','pairwise')

% Create figure
% figure1 = figure('Position',[1200 900 180 125]);
subplot(224)
mdl3 = fitlm(y3,x3,'linear','RobustOpts','on')
plot1=plot(mdl3,'MarkerSize',7,'Marker','.','MarkerEdgeColor',[0    0.4470    0.7410]);

set(plot1(2),'LineStyle','-','LineWidth',1);
set(plot1(3),'LineStyle',':','LineWidth',2);
set(plot1(4),'LineStyle',':','LineWidth',2);
box off
legend('off')
title ('')
xlabel('pDMN activity');
ylabel('pDMN-CON conn');
% title("r=" + string(round(Act_conn_corr2,2)) + "; p=" + string(round(Act_conn_corrp2,3)))


%% Fig 3D Activity Conn correlation
clc
y = aDMNnetlocorrm(:,6:7);
x1 = f(2).netlo(:,3);
% x1=filloutliers(x1,nan);
% y=filloutliers(y,nan);
na = sum(~isnan(y) & ~isnan(x1));
[Act_conn_corr1b,Act_conn_corrp1b]=corr(x1,y,'Type','Spearman','Rows','pairwise');

y = pDMNnetlocorrm(:,6:7);
x2 = f(2).netlo(:,4);
% x2=filloutliers(x2,nan);
% y=filloutliers(y,nan);
np = sum(~isnan(y) & ~isnan(x2));

[Act_conn_corr2b,Act_conn_corrp2b]=corr(x2,y,'Type','Spearman','Rows','pairwise');

y = mtDMNnetlocorrm(:,6:7);
x3 = f(2).netlo(:,5);
% x3=filloutliers(x3,nan);
% y=filloutliers(y,nan);
nm = sum(~isnan(y) & ~isnan(x3));

[Act_conn_corr3b,Act_conn_corrp3b]=corr(x3,y,'Type','Spearman','Rows','pairwise');

%d[~,~,~,allpc]=fdr_bh(allps)
savePath = 'A:\TwoTap\Manuscript\Figure';
sfile=[savePath '\tempData4Sec.xlsx'];

allRhoLo = [ Act_conn_corr1b;Act_conn_corr3b;Act_conn_corr2b];
allPLo = [Act_conn_corrp1b; Act_conn_corrp3b;Act_conn_corrp2b];
allNlo = [na; nm ;np];

% writetable(array2table(allRhoHi),sfile,'Sheet','ActNetCorr','WriteVariableNames',false, 'Range','B3')
% writetable(array2table(allPHi),sfile,'Sheet','ActNetCorr','WriteVariableNames',false, 'Range','B10')
%%
y = aDMNnetrestm(:,6:7);
x1 = f(2).netr(:,3);
% x1=filloutliers(x1,nan);
% y=filloutliers(y,nan);
na = sum(~isnan(y) & ~isnan(x1));

[Act_conn_corr1r,Act_conn_corrp1r]=corr(x1,y,'Type','Spearman','Rows','pairwise');

y = pDMNnetrestm(:,6:7);
x2 = f(2).netr(:,4);
% x2=filloutliers(x2,nan);
% y=filloutliers(y,nan);
np = sum(~isnan(y) & ~isnan(x2));

[Act_conn_corr2r,Act_conn_corrp2r]=corr(x2,y,'Type','Spearman','Rows','pairwise');

%
y = mtDMNnetlrestm(:,6:7);
x3 = f(2).netr(:,5);
% x3=filloutliers(x3,nan);
% y=filloutliers(y,nan);
nm = sum(~isnan(y) & ~isnan(x3));

[Act_conn_corr3r,Act_conn_corrp3r]=corr(x3,y,'Type','Spearman','Rows','pairwise');

allRhoRest = [ Act_conn_corr1r;Act_conn_corr3r;Act_conn_corr2r];
allPRest = [Act_conn_corrp1r; Act_conn_corrp3r;Act_conn_corrp2r];
allNRest = [na; nm; np];
% writetable(array2table(allRhoRest),sfile,'Sheet','ActNetCorr','WriteVariableNames',false, 'Range','J3')
% writetable(array2table(allPRest),sfile,'Sheet','ActNetCorr','WriteVariableNames',false, 'Range','J10')


%% Plot activity vs network Conn
% close all
bigData = cat(3,allRhoLo,allRhoHi,allRhoRest);
bigData = permute(bigData,[3,2,1]);

bigN = cat(3,allNlo,allNhi,allNRest);
bigN = permute(bigN,[3,2,1]);

% subplot(131)
% data=squeeze(bigData(:,:,1));
% sem=zeros(size(data));
% plotErrBar(data,sem,{'w';'k'})
% set(gca,'xticklabels',{'Dist', 'Atten','Rest'},'fontweight','bold','fontsize',16,'xticklabelrotation',90)
% title('aDMN')
% ylim([-0.35 0.15])

% subplot(132)
% data=squeeze(bigData(:,:,2));
% sem=zeros(size(data));
% plotErrBar(data,sem,{'w';'k'})
% set(gca,'xticklabels',{'Dist', 'Atten','Rest'},'fontweight','bold','fontsize',16,'xticklabelrotation',90)
% title('mtDMN')
% ylim([-0.35 0.15])

subplot(223)
data=squeeze(bigData(:,2,3));
sem=zeros(size(data));
plotErrBar(data,sem,{'k'})
set(gca,'Xtick',[1:3],'xticklabels',{'Dist', 'Atten','Rest'},'xticklabelrotation',90)
% title('pDMN')
ylim([-0.35 0.15])
xlim([-0.2 4.2])
title('CON-Connectivity')
% sgtitle('Network activity corr with network connectivity')
set(gcf, 'position', [663 313.6667 1.1987e+03 710.6667])
box off
% saveas(gcf,[savePath '\Fig4corrBar.tiff']);

%% Dunn and Clarkfs fisher z transform

% compare atten/dist to rest
%%% aDMN 
% distracted con vs dist Rest
% atten Con vs atten Rest
r1 = bigData(1,2,1); % distracted, CON, aDMN 
r2 = bigData(2,2,1); % attended, CON, aDMN 
rr = bigData(3,2,1); % rest, CON, aDMN 
n1 = bigN(1,2,1); % distracted, CON, aDMN 
n2 = bigN(2,2,1); % attended, CON, aDMN 
nr = bigN(3,2,1); % rest, CON, aDMN 
p1 = compare_correlation_coefficients(rr,r2,nr,n2);
p2 = compare_correlation_coefficients(r1,rr,n1,nr);
newP = fdr([p1 p2])

%%% mtDMN 
r1 = bigData(1,2,2); % distracted, CON, aDMN 
r2 = bigData(2,2,2); % attended, CON, aDMN 
rr = bigData(3,2,2); % rest, CON, aDMN 
n1 = bigN(1,2,2); % distracted, CON, aDMN 
n2 = bigN(2,2,2); % attended, CON, aDMN 
nr = bigN(3,2,2); % rest, CON, aDMN 
p1 = compare_correlation_coefficients(rr,r2,nr,n2);
p2 = compare_correlation_coefficients(r1,rr,n1,nr);
newP = fdr([p1 p2])

%% pDMN
clc
r1 = bigData(1,2,3); % distracted, FPN, pDMN 
r2 = bigData(2,2,3); % attended, FPN, pDMN 
rr = bigData(3,2,3); % rest, FPN, pDMN 
n1 = bigN(1,2,3); % distracted, FPN, pDMN 
n2 = bigN(2,2,3); % attended, FPN, pDMN 
nr = bigN(3,2,3); % rest, FPN, pDMN 
p1 = compare_correlation_coefficients(rr,r2,nr,n2);
p2 = compare_correlation_coefficients(r1,rr,n1,nr);
newP = fdr([p1 p2])

%% Compare between networks

% aDMN
r1 = bigData(1,1,1); % distracted, FPN, aDMN 
r2 = bigData(1,2,1); % Distracted, CON, aDMN 
r3 = bigData(2,1,1); % attended, FPN, aDMN 
r4 = bigData(2,2,1); % attended, CON, aDMN 

n1 = bigN(1,1,1); % distracted, FPN, aDMN 
n2 = bigN(1,2,1); % Distracted, CON, aDMN 
n3 = bigN(2,1,1); % attended, FPN, aDMN 
n4 = bigN(2,2,1); % attended, CON, aDMN 

p1 = compare_correlation_coefficients(r1,r2,n1,n2);
p2 = compare_correlation_coefficients(r3,r4,n3,n4);
[p1 p2]
% newP = fdr([p1 p2])

%%% mtDMN 
r1 = bigData(1,1,2); % distracted, FPN, mtDMN 
r2 = bigData(1,2,2); % Distracted, CON, mtDMN 
r3 = bigData(2,1,2); % attended, FPN, mtDMN 
r4 = bigData(2,2,2); % attended, CON, mtDMN 

n1 = bigN(1,1,2); % distracted, FPN, mtDMN 
n2 = bigN(1,2,2); % Distracted, CON, mtDMN 
n3 = bigN(2,1,2); % attended, FPN, mtDMN 
n4 = bigN(2,2,2); % attended, CON, mtDMN 

p1 = compare_correlation_coefficients(r1,r2,n1,n2);
p2 = compare_correlation_coefficients(r3,r4,n3,n4);
[p1 p2]
% newP = fdr([p1 p2])

%%% pDMN
r1 = bigData(1,1,3); % distracted, FPN, pDMN 
r2 = bigData(1,2,3); % Distracted, CON, pDMN 
r3 = bigData(2,1,3); % attended, FPN, pDMN 
r4 = bigData(2,2,3); % attended, CON, pDMN 

n1 = bigN(1,1,3); % distracted, FPN, pDMN 
n2 = bigN(1,2,3); % Distracted, CON, pDMN 
n3 = bigN(2,1,3); % attended, FPN, pDMN 
n4 = bigN(2,2,3); % attended, CON, pDMN 

p1 = compare_correlation_coefficients(r1,r2,n1,n2);
p2 = compare_correlation_coefficients(r3,r4,n3,n4);
[p1 p2]
% newP = fdr([p1 p2])


%% rest relationship
% Fig 3D Activity Conn correlation

y = aDMNnetrestm(:,6:7);
x1 = f(2).netr(:,3);
[Act_conn_corr1r,Act_conn_corrp1r]=corr(x1,y,'Type','Spearman','Rows','pairwise');

y = pDMNnetrestm(:,6:7);
x2 = f(2).netr(:,4);
[Act_conn_corr2r,Act_conn_corrp2r]=corr(x2,y,'Type','Spearman','Rows','pairwise');

%
y = mtDMNnetlrestm(:,6:7);
x3 = f(2).netr(:,5);
[Act_conn_corr3r,Act_conn_corrp3r]=corr(x3,y,'Type','Spearman','Rows','pairwise');

allps=[Act_conn_corrp1,Act_conn_corrp2,Act_conn_corrp3,Act_conn_corrp1b,Act_conn_corrp2b,Act_conn_corrp3b,Act_conn_corrp1r,Act_conn_corrp2r,Act_conn_corrp3r]

keyboard;

%% adding in age co-variate

y = aDMNnetlocorrm(:,6:7);
x1 = f(2).netlo(:,3);
[Act_conn_corr1,Act_conn_corrp1]=partialcorr(x1,y,age,'Type','Spearman','Rows','pairwise');

y = pDMNnetlocorrm(:,6:7);
x2 = f(2).netlo(:,4);
[Act_conn_corr2,Act_conn_corrp2]=partialcorr(x2,y,age,'Type','Spearman','Rows','pairwise');

%
y = mtDMNnetlocorrm(:,6:7);
x3 = f(2).netlo(:,5);
[Act_conn_corr3,Act_conn_corrp3]=partialcorr(x3,y,age,'Type','Spearman','Rows','pairwise');

y = aDMNnethicorrm(:,6:7);
x1 = f(2).nethi(:,3);
[Act_conn_corr1b,Act_conn_corrp1b]=partialcorr(x1,y,age,'Type','Spearman','Rows','pairwise');
%
y = pDMNnethicorrm(:,6:7);
x2 = f(2).nethi(:,4);
[Act_conn_corr2b,Act_conn_corrp2b]=partialcorr(x2,y,age,'Type','Spearman','Rows','pairwise');

%
y = mtDMNnethicorrm(:,6:7);
x3 = f(2).nethi(:,5);
[Act_conn_corr3b,Act_conn_corrp3b]=partialcorr(x3,y,age,'Type','Spearman','Rows','pairwise');
%[~,~,~,Act_conn_corrp3b]=fdr_bh(Act_conn_corrp3b)
%%
allps=[Act_conn_corrp1,Act_conn_corrp2,Act_conn_corrp3,Act_conn_corrp1b,Act_conn_corrp2b,Act_conn_corrp3b]


%% 
% for i = 1:5
%     tmp = pDMNnethicorr;%-FPNnethicorr;    
%     tmp=squeeze(tmp(2,i,:,:));
%     tmpcorr=corr(mind,tmp,'Type','Spearman','Rows','pairwise');
%     plot(t3, tmpcorr);  hold on; 
%     plot([-3.5;3.5],[.12;.12],'r--');
%     plot([-3.5;3.5],[-.12;-.12],'r--');
%     axis([-inf inf -0.2 .2]);
%     title(['Network ' num2str(i)]);
%     pause;
%     close all;
% end
 