%% Repiration test for EEGO
%%%%%%%%%%%%%% ONLY USE FOR EEGO TWOTAP REVISIONS - NOT COMPATIBLE WITH REGULAR
%%%%%%%%%%%%%% BRAINE PROCESSING
%%%%%%%%%%%%%% MAJOR CHANGES INCLUDE EPOCHING AND DEFINING PROCESSVARS
% Jason Nan
cd('C:\Users\jason\OneDrive\Desktop\MS Beng\Neatlabs\Respiration_EEG_Synchronization')
addpath('A:\eeglab_OldLSL_DataAna04072023')
addpath('Scripts')
eeglab

%% 
% id: 2 - didn't do task right 
% id: 7 - double pressed 
% id: 15 - error isn't normal behavior 
% id: 18 - error isn't normal behavior - user paused
subIds = [1,3,5,6,7,9,10,11,13,14,15,16,17,18,19,20];
count=0;

id = 5

rawFile = "A:\TwoTap\twotap"+string(id)+"\twotap"+string(id)+".xdf";
behFile = pickfiles({"A:\TwoTap\twotap"+string(id)}, {"trial_summary.csv"});
codeFile = pickfiles({"A:\TwoTap\twotap"+string(id)}, {".csv"},{".csv"},{'trial_summary','fixed_output'});
behTbl = readtable(behFile);
codeTbl = readtable(codeFile);

if id==18
    behTbl.ResponseTime(36) = 9450;
end
% Fixing behTbl
responseTbl = codeTbl(codeTbl.code==11,:); % response generated from codes table
diffRows = abs(size(responseTbl,1)-size(behTbl,1));
if size(responseTbl,1)~=size(behTbl,1)
    count=0;
    for i=1:size(responseTbl,1)% go through all codes in responseTbl
        while abs(responseTbl.duration(i) - behTbl.ResponseTime(i)) >10 & abs(responseTbl.duration(i+1) - behTbl.ResponseTime(i+1)) >10 % buffer of 10mS
            % as long as the trial in behTbl is not the same as the
            % response and the next one aren't the same
            % remove that row. "i" stays the same so it will check the next row once the
            % previous one is removed
            behTbl.ResponseTime(i)=behTbl.ResponseTime(i)+behTbl.ResponseTime(i+1);
            behTbl(i+1,:)=[];
            count=count+1;
        end
    end
%     writetable(behTbl, replace(behFile,'trial_summary','fixed_output'))
end
% Checking if the number of RT differences is the same as number of
% different rows initially
if count ~= diffRows
    error('Check behTbl')
end


%%
med = median(behTbl.ResponseTime);
singleBreath=med/2;
madRT = mad(behTbl.ResponseTime);

attendedPresses = behTbl.Accuracy==1 & behTbl.ResponseTime<med+madRT & behTbl.ResponseTime>med-madRT;
slowPresses = behTbl.Accuracy==1 & behTbl.ResponseTime>med+madRT;

%% Calculating behavior trials
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
% respiration = lowpass(respiration,1,10); % lets in breaths > 2 sec
respiration = sgolayfilt(respiration, 2, filterLen);  % Smooth the signal

% Interpolation 
% respiration = interp1(timeSeries,respiration, timeSeries(1):0.02:timeSeries(end));
% timeSeries = timeSeries(1):0.02:timeSeries(end);


% median(behTbl.ResponseTime)/2000

dydx = gradient(respiration(:)) ./ gradient(timeSeries(:)); % calculate derivitive
dydx = lowpass(dydx,0.9,10); % lets in breaths > 2 sec
% dydx = sgolayfilt(dydx, 2, filterLen);  % Smooth the signal

[~,~,zeroCross]=zerocrossrate(dydx);

d2y = gradient(dydx(:)) ./ gradient(timeSeries(:)); % inflection point

respirationPhase = angle(hilbert(respiration)); % calculate phase 

% Plotting respiration data
pressesLoc = find(markers==11); % all button presses
pressedTime=markersTime(pressesLoc);

slowTimes = pressedTime(slowPresses);
attendedTimes = pressedTime(attendedPresses);
[sm,sl]=min(abs(timeSeries-slowTimes')');
[am,al]=min(abs(timeSeries-attendedTimes')');

% close all
figure
subplot(311)
hold on
plot(timeSeries,(respiration),'k')
plot(slowTimes,respiration(sl),'b*')
plot(attendedTimes,respiration(al),'r*')
ylabel('Respiration Signal')

subplot(312)
hold on
plot(timeSeries,movmean(dydx,1),'k')
plot(slowTimes,dydx(sl),'b*')
plot(attendedTimes,dydx(al),'r*')
plot(timeSeries,zeros(size(timeSeries)),'k--')
plot(timeSeries(zeroCross),dydx(zeroCross),'color','#D95319','marker','+','linestyle','none')

ylabel('Derivative of respiration Signal')

subplot(313)
hold on
plot(timeSeries,(respirationPhase),'k')
plot(slowTimes,respirationPhase(sl),'b*')
plot(attendedTimes,respirationPhase(al),'r*')
ylabel('Phase')
plot(timeSeries,zeros(size(timeSeries)),'k--')

%% Fitting perfect sinusoid

% Detect start and end points of breath cycles

% Segment breath cycles

% Curve fitting and breath width calculation



%% Single trial
attendedTrials = find(attendedPresses);

for t=1:length(attendedTrials)

    RT = behTbl.ResponseTime(attendedTrials(t)); % get response time from behTbl
    t1=pressedTime(attendedTrials(t)); % time stamp of button press
    t0 = t1-RT/1000;% get previous button press
    t2 = t1+behTbl.ResponseTime(attendedTrials(t)+1)/1000; % next Button press
    [~,pb]=min(abs(timeSeries'-t0)); % index of previous button press relative to all time
    [~,cb]=min(abs(timeSeries'-t1)); % index of current button press relative to all time
    [~,nb]=min(abs(timeSeries'-t2)); % index of current button press relative to all time

    nextBuffer = round(RT/200); % go back half a breath
    prevBuffer = round(RT/200); % go back half a breath
    loc0=pb-prevBuffer;
    loc1=cb+nextBuffer;
    timeSection = timeSeries(loc0:loc1);% get section of time stamps
    respSection = streams{idxRespiration}.time_series(loc0:loc1);
    derSection = dydx(loc0:loc1); % get section of derivative
    ddy =  gradient(derSection(:)) ./ gradient(timeSection(:)); % calculate inflection point in case there is no zero crossing
%     ddy = d2y(loc0:loc1);

    ft = fittype('gauss8');  % Gaussian curve fitting
    [coeffs, gof,output] = fit(timeSection', double(respSection)', ft);
%     derSection = differentiate(coeffs,timeSection);
    obj = @(x) filloutliers(coeffs(x),'spline');
    derSection =  gradient(obj(timeSection)) ./ gradient(timeSection(:)); % calculate inflection point in case there is no zero crossing
    [npeaks, peakLocs] = findpeaks(coeffs(timeSection),'MinPeakDistance',RT/400); % number of peaks in signal seperated by at least 1/2 breath
    
    




    [~,~,zeroCross]=zerocrossrate(derSection);
    extrazeroCross = abs(derSection)<0.1;
    zeroCross = zeroCross | extrazeroCross';
    if zeroCross(1)==1
        zeroCross(1)=[];
    end

    [~,~,ddyZero]=zerocrossrate(ddy);
    if ddyZero(1)==1
        ddyZero(1)=[];
    end
    
    ddyZero = zeros(size(respSection));
    [~,loc]=findpeaks(ddy,'NPeaks',20,'MinPeakHeight',0.5);
    ddyZero(loc)=1;
    ddyZero=logical(ddyZero);


    
    [~,peakLoc]=findpeaks(zeroCross*1); % finding clusters of zero crossings

    % Assigning zero crossing to clusters for derivative
    zero_crossing_idx = find(zeroCross);
    cluster_idx = ones(size(zero_crossing_idx));
    cluster_id = 1;
    for i = 2:length(zero_crossing_idx)
        if zero_crossing_idx(i) - zero_crossing_idx(i-1) > 1
            cluster_id = cluster_id + 1;
        end
        cluster_idx(i) = cluster_id;
    end
    
    % Find the closest point to zero for each cluster
    unique_clusters = unique(cluster_idx);
    closest_zero_idx = zeros(length(unique_clusters), 1);
    for i = 1:length(unique_clusters)
        cluster_members = zero_crossing_idx(cluster_idx == unique_clusters(i));
        [~, min_idx] = min(abs(derSection(cluster_members)));
        closest_zero_idx(i) = cluster_members(min_idx);
    end
    
    % Check if there are 5 zeros, if not, use ddy to find them

    if length(closest_zero_idx)~=5

    end


%     respSection= respSection-mean(respSection);
%     toppeaks=maxk(respSection,2);
%     lowpeaks= maxk(-respSection,2);'
%     high_data=respSection>mean(toppeaks)*0.1;
%     low_data=respSection<mean(-lowpeaks)*.9;


%     close all
%     figure
    clf


   
%     hold on
%     plot(high_data,'linewidth',2,'color','b')
%     plot(low_data,'r--','linewidth',2)
%     hold on
%     yyaxis right
%     plot(respSection,'k')
%     [onsets,offsets] = breathTimes(double(respSection),10)
%     findchangepts(respSection,'MaxNumChanges',5,'MinDistance',5,'Statistic','mean','min')
%     pulsewidth(double(respSection),timeSection)
% 
    hold on
    % fitted data plots
    plot(timeSection,respSection,'LineWidth',3)
%     plot(timeSection, obj(timeSection)) % grubbs cleaned sig

    plot(timeSection,coeffs(timeSection)) % raw fit


%     plot(timeSection,movmean(respSection,5),'LineWidth',3)
%     plot(timeSection,derSection,'LineWidth',1.5)
%     plot(timeSection,ddy,'LineWidth',1.5,'color','r')
% 
    plot(t1,streams{idxRespiration}.time_series(cb),'b*','markersize',10,'LineWidth',3)
    plot(t0,streams{idxRespiration}.time_series(pb),'b*','markersize',10,'LineWidth',3)
%     plot(timeSection,zeros(size(timeSection)),'k--')
% %     yyaxis right
%     hold on
% %     plot(timeSection(zeroCross),derSection(zeroCross),'color','#D95319','marker','+','linestyle','none')
% %     yyaxis right 
% %     hold on
%     plot(timeSection,zeros(size(timeSection)),'k--')
%     plot(timeSection(ddyZero),respSection(ddyZero),'color','k','marker','+','linestyle','none','markersize',20,'LineWidth',3)
    plot(timeSection(closest_zero_idx),respSection(closest_zero_idx),'color','#D95319','marker','+','linestyle','none','markersize',20,'LineWidth',3)
    plot(timeSection(peakLocs), respSection(peakLocs),'r*','markersize',20)
%     legend('Respiration','dy','d2y','','','','','inflection point','peaks')

    waitforbuttonpress
%     RT(t) = t1-t0;
%     
end


%% Entire signal fit

obj = @(x) filloutliers(coeffs(x),'spline');
plot(timeSection, obj(timeSection))
%% Phase anlaysis
figure
[~,loc]=min(abs(timeSeries'-pressedTime));

slowAngles=respirationPhase(loc(slowPresses));
fastAngles=respirationPhase(loc(attendedPresses));
% histogram(slowAngles,15)
% hold on
% histogram(fastAngles,15)
% legend('slow','fast')
% xlabel('phase angle')
% ylabel('Count')

%
close all
% Define your vector of angles (in radians)
% Set the radius of the circle
radius = 1; % You can adjust this value as needed
% Create a polar plot
figure;
polarplot(slowAngles, radius*ones(size(slowAngles)), 'bo'); % slow
hold on;
polarplot(fastAngles, radius*ones(size(fastAngles)), 'ro'); % fast
title('Angles of Slow and Attended Trials');
grid on;
set(gca,'rticklabels','')
rlim([0 1.01])
legend('slow','fast')


%% Pre processing EEG
[EEGset, ProcessVars] = preProcessing_EEGO_twotap(rawFile, behFile);
%%% Need to redo this with new behTbl
% slowTrialsNeural  = ProcessVars.behtabledata.Accuracy==1  & ismember(1:size(ProcessVars.behtabledata,1),ProcessVars.trialset(2).ind)';
% fastTrialsNeural  = ProcessVars.behtabledata.Accuracy==1  & ismember(1:size(ProcessVars.behtabledata,1),ProcessVars.trialset(1).ind)';