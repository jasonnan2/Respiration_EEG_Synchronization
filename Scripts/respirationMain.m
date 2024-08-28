%% Repiration test for EEGO
%%%%%%%%%%%%%% ONLY USE FOR EEGO TWOTAP REVISIONS - NOT COMPATIBLE WITH REGULAR
%%%%%%%%%%%%%% BRAINE PROCESSING
%%%%%%%%%%%%%% MAJOR CHANGES INCLUDE EPOCHING AND DEFINING PROCESSVARS
% Jason Nan
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization')
addpath('A:\eeglab_OldLSL_DataAna04072023')
addpath(genpath('Scripts'))
eeglab

%% 
% id: 2 - didn't do task right 
% id: 7 - double pressed 
% id: 15 - error isn't normal behavior 
% id: 18 - error isn't normal behavior - user paused
subIds = [1,3,5,6,7,9,10,11,13,14,15,16,17,18,19,20];
count=0;
avgFast = []; avgSlow=[];consistency=[];
for id = subIds
    count=count+1;

    rawFile = "A:\TwoTap\rawData\twotap"+string(id)+"\twotap"+string(id)+".xdf";
    behFile = pickfiles({"A:\TwoTap\rawData\twotap"+string(id)}, {"trial_summary.csv"});
    codeFile = pickfiles({"A:\TwoTap\rawData\twotap"+string(id)}, {".csv"},{".csv"},{'trial_summary','fixed_output'});
    behTbl = readtable(behFile);
    codeTbl = readtable(codeFile);
    if id==18
        behTbl.ResponseTime(36) = 9450;
    end
    [behTbl, responseTbl]=fixTwoTapBehFile(behTbl,codeTbl); % Fixing behTbl
    % Behavior Variables
    med = median(behTbl.ResponseTime);
    singleBreath=med/2;
    madRT = mad(behTbl.ResponseTime,1);
    COV_m = mean(behTbl.ResponseTime(behTbl.ResponseTime> 0));
    COV_sd = std(behTbl.ResponseTime(behTbl.ResponseTime>0));
    consistency(count) = (1-(COV_sd / COV_m));
    respRate(count)=med;

    
    attendedPresses = behTbl.Accuracy==1 & behTbl.ResponseTime<med+madRT & behTbl.ResponseTime>med-madRT;
    slowPresses = behTbl.Accuracy==1 & behTbl.ResponseTime>med+madRT;
    
    % Calculating behavior trials
    streams = load_xdf(rawFile);
    for s=1:length(streams)
        if strcmp(streams{s}.info.type,'RB')
            idxRespiration = s;
        elseif strcmp(streams{s}.info.type,'EEG')
            idxEEG = s;
        elseif strcmp(streams{s}.info.type,'Markers')
            idxMarker = s;
        end
    end
    % reading in data
    respiration=double(streams{idxRespiration}.time_series); % respiration is sampled at 10Hz
    timeSeries=streams{idxRespiration}.time_stamps; 
    
    markers=streams{idxMarker}.time_series;
    markersTime=streams{idxMarker}.time_stamps;
    
    filterLen = 10*ceil(median(behTbl.ResponseTime/8000))-1;
    respiration = highpass(respiration,0.04,10); % cutoff at 0.04Hz - lets in breaths <25 seconds
    respiration = movmean(respiration,5); % Very simple low pass filter
    respiration = sgolayfilt(respiration, 2, filterLen);  % Smooth the signal
    
    dydx = gradient(respiration(:)) ./ gradient(timeSeries(:)); % calculate derivitive
    dydx = lowpass(dydx,0.9,10); % lets in breaths > 2 sec

    % Getting Button presses
    pressesLoc = find(markers==11); % all button presses
    pressedTime=markersTime(pressesLoc);
    slowTimes = pressedTime(slowPresses);
    attendedTimes = pressedTime(attendedPresses);
    [sm,sl]=min(abs(timeSeries-slowTimes')');
    [am,al]=min(abs(timeSeries-attendedTimes')');

    %%% Phase Stuff 
    respirationPhase = angle(hilbert(respiration)); % calculate phase 
    [~,loc]=min(abs(timeSeries'-pressedTime));
    slowAngles=respirationPhase(loc(slowPresses));
    fastAngles=respirationPhase(loc(attendedPresses));
    avgFast(count) = circ_mean(fastAngles'); % save phase data
    avgSlow(count) = circ_mean(slowAngles'); % save phase data
    %%%

    % Single trial
    attendedTrials = find(attendedPresses);
    distractedTrials = find(slowPresses);

    [discardedAttended(count), attendedI(count), attendedE(count), attenedCycle(count), attendedRate(count)]=calculateBreathCycles(pressedTime,timeSeries,streams,  idxRespiration, behTbl, attendedTrials,0);
    [discardedDistracted(count), distractedI(count), distractedE(count), distractedCycle(count),distractedRate(count)]=calculateBreathCycles(pressedTime,timeSeries,streams,  idxRespiration, behTbl, distractedTrials,0);

    totalAttended(count)=length(attendedTrials);
    totalDistracted(count)=length(distractedTrials);
end

clc
1
%%
keptDistracted = totalDistracted'-discardedDistracted';
% 100*(discardedDistracted./totalDistracted)'
keptAttended = totalAttended'-discardedAttended';
% 100*(discardedAttended./totalAttended)'

% removing any subs with <3 trials
attendedI(keptAttended<3)=nan;
attendedE(keptAttended<3)=nan;
attenedCycle(keptAttended<3)=nan;

distractedI(keptDistracted<3)=nan;
distractedE(keptDistracted<3)=nan;
distractedCycle(keptDistracted<3)=nan;

[~,pI]=ttest(attendedI,distractedI)
[~,pE]=ttest(attendedE,distractedE)
[~,pBC]=ttest(attenedCycle,distractedCycle)
% [~,pdiff]=ttest(attendedE-attendedI,distractedE - distractedI)
% [~,pdiff]=ttest(attendedE+attendedI,distractedE + distractedI)

breathTbl=table();
breathTbl.attendedI = attendedI'; 
breathTbl.attendedE = attendedE'; 
breathTbl.attendedCycle = attenedCycle'; 
breathTbl.distractedI = distractedI'; 
breathTbl.distractedE = distractedE'; 
breathTbl.distractedCycle = distractedCycle'; 

% 
% writetable(breathTbl,'A:\TwoTap\Manuscript\Figure\tempData4Sec.xlsx','sheet','breathData')
%% Plotting breath Response

breathTbl = readtable('A:\TwoTap\Manuscript\Figure\tempData4Sec.xlsx','sheet','breathData');
dataMat = table2array(breathTbl);
dataMat = dataMat(~any(isnan(dataMat),2),:);
% avgdata = nanmean(dataMat,1);
% data = [avgdata(4:6); avgdata(1:3)];
% sem = nanstd(dataMat,[],1)./sqrt([16,16,16,12,12,12]);
% sem = [sem(4:6);sem(1:3)];
% close all
% plotErrBar(data',sem',{'b','r'})

% plot(meanWeight,'-o')
subplot(2,3,[4,5])
hold on
allI = reshape(dataMat(:,[4,1]),[],1);
allE = reshape(dataMat(:,[5,2]),[],1);
netOrder =(repmat(reshape(repmat([1:2],[12],1),[],1),2,1));
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
bbb=boxchart(netOrder, [allI;allE],'GroupByColor',[ones(1,12*2), 2*ones(1,12*2)]','MarkerStyle','none')

x = [ones(12,1)-0.25 + 0.1*(rand(12,1)-0.5) ones(12,1)+0.25 + 0.1*(rand(12,1)-0.5)];
for tt=1:12
    plot(x(tt,:),dataMat(tt,[4,1]),'ko-')
end

x = [2*ones(12,1)-0.25 + 0.1*(rand(12,1)-0.5) 2*ones(12,1)+0.25 + 0.1*(rand(12,1)-0.5)];
for tt=1:12
    plot(x(tt,:),dataMat(tt,[5,2]),'ko-')
end
set(gca,'XTickLabel',{'avg Inhalation','avg Exhalation'})
ylabel('Seconds')


subplot(2,3,6)
hold on
allBC = reshape(dataMat(:,[6,3]),[],1);
netOrder =(repmat(reshape(repmat([1:1],[12],1),[],1),2,1));
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
bbb=boxchart(netOrder, [allBC],'GroupByColor',[ones(1,12*1),2*ones(1,12*1) ]','MarkerStyle','none')

x = [ones(12,1)-0.25 + 0.1*(rand(12,1)-0.5) ones(12,1)+0.25 + 0.1*(rand(12,1)-0.5)];
for tt=1:12
    plot(x(tt,:),dataMat(tt,[6,3]),'ko-')
end
ylabel('Breath Cycle')
legend({'Low Cons','High Cons'},'box','off')


%%
% close all
% [asdf, avgI, avgE, avgBreathCycle]=calculateBreathCycles(pressedTime,timeSeries,streams,  idxRespiration, behTbl, distractedTrials,0)
% calculateBreathCycles(pressedTime,timeSeries,streams,  idxRespiration, behTbl, attendedTrials)
%% Phase anlaysis
% figure

clc
phaseP= circ_wwtest(avgFast,avgSlow) % ttest between slow and fast angles
% close all
% Define your vector of angles (in radians)
% Set the radius of the circle
radius = 1; % You can adjust this value as needed
% Create a polar plot
% figure;
subplot(233)
polarplot(avgSlow, radius*ones(size(avgSlow)), 'b*'); % slow
hold on;
polarplot(avgFast, radius*ones(size(avgFast)), 'r*'); % fast
% title('Distracted and Attended Trials across Subjects');
grid on;
set(gca,'rticklabels','')
rlim([0 1.01])
% legend('Low Cons','High Cons')

%% Sample plots for Illustration


id=9;
rawFile = "A:\TwoTap\rawData\twotap"+string(id)+"\twotap"+string(id)+".xdf";
behFile = pickfiles({"A:\TwoTap\rawData\twotap"+string(id)}, {"trial_summary.csv"});
codeFile = pickfiles({"A:\TwoTap\rawData\twotap"+string(id)}, {".csv"},{".csv"},{'trial_summary','fixed_output'});
behTbl = readtable(behFile);
codeTbl = readtable(codeFile);
if id==18
    behTbl.ResponseTime(36) = 9450;
end
[behTbl, responseTbl]=fixTwoTapBehFile(behTbl,codeTbl); % Fixing behTbl
% Behavior Variables
med = median(behTbl.ResponseTime);
singleBreath=med/2;
madRT = mad(behTbl.ResponseTime,1);

attendedPresses = behTbl.Accuracy==1 & behTbl.ResponseTime<med+madRT & behTbl.ResponseTime>med-madRT;
slowPresses = behTbl.Accuracy==1 & behTbl.ResponseTime>med+madRT;

% Calculating behavior trials
streams = load_xdf(rawFile);
for s=1:length(streams)
    if strcmp(streams{s}.info.type,'RB')
        idxRespiration = s;
    elseif strcmp(streams{s}.info.type,'EEG')
        idxEEG = s;
    elseif strcmp(streams{s}.info.type,'Markers')
        idxMarker = s;
    end
end
% reading in data
respiration=double(streams{idxRespiration}.time_series); % respiration is sampled at 10Hz
rawRespiration =respiration; 
timeSeries=streams{idxRespiration}.time_stamps; 

markers=streams{idxMarker}.time_series;
markersTime=streams{idxMarker}.time_stamps;

filterLen = 10*ceil(median(behTbl.ResponseTime/8000))-1;
respiration = highpass(respiration,0.04,10); % cutoff at 0.04Hz - lets in breaths <25 seconds
respiration = movmean(respiration,5); % Very simple low pass filter
respiration = sgolayfilt(respiration, 2, filterLen);  % Smooth the signal

dydx = gradient(respiration(:)) ./ gradient(timeSeries(:)); % calculate derivitive
dydx = lowpass(dydx,0.9,10); % lets in breaths > 2 sec

% Getting Button presses
pressesLoc = find(markers==11); % all button presses
pressedTime=markersTime(pressesLoc);
slowTimes = pressedTime(slowPresses);
attendedTimes = pressedTime(attendedPresses);
[sm,sl]=min(abs(timeSeries-slowTimes')');
[am,al]=min(abs(timeSeries-attendedTimes')');

%%% Phase Stuff 
respirationPhase = angle(hilbert(respiration)); % calculate phase 
[~,loc]=min(abs(timeSeries'-pressedTime));
slowAngles=respirationPhase(loc(slowPresses));
fastAngles=respirationPhase(loc(attendedPresses));
avgFast(count) = circ_mean(fastAngles'); % save phase data
avgSlow(count) = circ_mean(slowAngles'); % save phase data
%%%

%%
close all
figure
subplot(2,3,[1,2])
hold on
plot(timeSeries,(rawRespiration),'k')
plot(slowTimes,rawRespiration(sl),'b*')
plot(attendedTimes,rawRespiration(al),'r*')
ylabel('Respiration Signal')
xlim([3247330 3247350])
ax=gca();
set(gca,'xtick',ax.XTick,'xticklabels',ax.XTick-ax.XTick(1))
legend({'','Low Cons','High Cons'},'box','off')
xlabel('Time (sec)')

% subplot(312)
% hold on
% plot(timeSeries,movmean(dydx,1),'k')
% plot(slowTimes,dydx(sl),'b*')
% plot(attendedTimes,dydx(al),'r*')
% plot(timeSeries,zeros(size(timeSeries)),'k--')
% plot(timeSeries(zeroCross),dydx(zeroCross),'color','#D95319','marker','+','linestyle','none')
% xlim([3.2473e+06 3.2473e+06+86])
% ylabel('Derivative of respiration Signal')
% 
% subplot(212)
% hold on
% plot(timeSeries,(respirationPhase),'k')
% plot(slowTimes,respirationPhase(sl),'b*')
% plot(attendedTimes,respirationPhase(al),'r*')
% ylabel('Phase (rad)')
% xlabel('Time (sec)')
% plot(timeSeries,zeros(size(timeSeries)),'k--')
% xlim([3.2473e+06 3.2473e+06+86])
% ax=gca();
% set(gca,'xticklabels',ax.XTick-ax.XTick(1),'fontweight','bold','fontsize',12)
% 



%% Testing bimodal distribution of data

gp1 = breathTbl.distractedCycle>2;
gp2 = breathTbl.distractedCycle<2;


[~,pBC1]=ttest(breathTbl.attendedCycle(gp1),breathTbl.distractedCycle(gp1))
[~,pBC2]=ttest(breathTbl.attendedCycle(gp2),breathTbl.distractedCycle(gp2))

[~,pBC1]=ttest(breathTbl.attendedCycle(gp1)-breathTbl.distractedCycle(gp1))
[~,pBC2]=ttest(breathTbl.attendedCycle(gp2)-breathTbl.distractedCycle(gp2))

[~,pBC2]=ttest(abs(breathTbl.attendedCycle-breathTbl.distractedCycle))

[pBC2,~,stats]=signrank(breathTbl.attendedCycle,breathTbl.distractedCycle)

clc
[~,pI,~,statsI]=ttest(breathTbl.attendedI(gp1),breathTbl.distractedI(gp1))
[~,pE,~,statsE]=ttest(breathTbl.attendedE(gp1),breathTbl.distractedE(gp1))





%% 64 vs 24 channel validations

load('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization\collated\SourceValidation.mat')
hm = headModel.loadFromFile('headModel_templateFile_newPEBplus.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';
labelnames = {''};
var = atanh(fisherZ);
plot68roi(hm, T'*var',[0 2],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\figS1_fisherZ_allScale.tiff']);
%%
var = alphas';
plot68roi(hm, T'*var',[0 1],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\figS1_zchronbach.tiff']);

var = alphas2;
plot68roi(hm, T'*var',[0 1],labelnames);
set(gcf,'Position',[985 611 390 523.3333]); 
saveas(gcf,[savePath '\figS1_chronbach2.tiff']);

%% Plotting distribution of corr and chronbach alphas
close all
figure
hold on
validationPlot = atanh(fisherZ);
x = repmat(1,68,1);
boxchart(validationPlot,'MarkerStyle','none')
swarmchart(x,validationPlot,50,'k.')
% ylim([0 1])
ylabel('Fisher Z')
set(gca,'xticklabels',{ 'Spearmans'},'xticklabelrotation',45)
set(gcf,'position',[1.1883e+03 1.1043e+03 140.0333 233.6667])

saveas(gcf,[savePath '\figS1_fisherZbox.tiff']);

%%
rMat; % size 16 x 68
close all
for n=[1,3,4,2,5]
    if n ==2 
        subplot(2,3,[4,4.5])
    elseif n ==5
        subplot(2,3,[5.5,6])
    elseif n==3
        subplot(2,3,2)
    elseif n==4
        subplot(2,3,3)
    else
        subplot(2,3,n)
    end
    hold on
    rois = netwrk(n).roi;
    nROI = length(rois);
    roiData = rMat(:,rois);
    roiName = atlas.label(rois);

    x = repmat(1:nROI,16,1);
    boxchart(roiData,'markerstyle','none')
    swarmchart(x,roiData,50,'k.')
    ylim([0.8,1])
    title(netwrk(n).name)
    set(gca,'xticklabels',roiName,'xticklabelrotation',45,'fontsize',12)

end

%% pDMN low vs High in 16 subs
close all
figure
hold on
ROIsplot = 1E13*[allLow'; allHigh' ];

netOrder = ones(1,32);
netOrder = categorical(arrayfun(@num2str, netOrder, 'UniformOutput', 0));
boxchart(netOrder', ROIsplot,'GroupByColor',[ones(1,16), 2*ones(1,16)],'MarkerStyle','none')

% data=round(data,3);
x = repmat(1,16,1)-0.25;
x=reshape(x,[],1);
swarmchart(x,1e13*allLow,50,'.','XJitterWidth',0.25,'MarkerEdgeColor',[0 0.4470 0.7410]);
x = repmat(1,16,1)+0.25;
x=reshape(x,[],1);
swarmchart(x,1e13*allHigh,50,'.','XJitterWidth',0.25,'MarkerEdgeColor',[0.8500 0.3250 0.0980]);

% ylabel('lPC')
title('pDMN 16 Subs low high')




