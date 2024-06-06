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

for id = 1
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
    madRT = mad(behTbl.ResponseTime);
    
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
    
    %%% Phase Stuff 
    respirationPhase = angle(hilbert(respiration)); % calculate phase 
    [~,loc]=min(abs(timeSeries'-pressedTime));
    slowAngles=respirationPhase(loc(slowPresses));
    fastAngles=respirationPhase(loc(attendedPresses));
    avgFast(count) = circ_mean(fastAngles'); % save phase data
    avgSlow(count) = circ_mean(slowAngles'); % save phase data
    %%%

    % Plotting respiration data
    pressesLoc = find(markers==11); % all button presses
    pressedTime=markersTime(pressesLoc);
    slowTimes = pressedTime(slowPresses);
    attendedTimes = pressedTime(attendedPresses);
    [sm,sl]=min(abs(timeSeries-slowTimes')');
    [am,al]=min(abs(timeSeries-attendedTimes')');

    % Single trial
    attendedTrials = find(attendedPresses);
    distractedTrials = find(slowPresses);

    discardedTrials=0;
    for t=1:length(attendedTrials)
        toplot=1;
        try
            RT = behTbl.ResponseTime(attendedTrials(t)); % get response time from behTbl
            t1=pressedTime(attendedTrials(t)); % time stamp of button press
            t0 = t1-RT/1000;% get previous button press
            [~,pb]=min(abs(timeSeries'-t0)); % index of previous button press relative to all time
            [~,cb]=min(abs(timeSeries'-t1)); % index of current button press relative to all time
        
            Buffer = round(RT/80); % go back half a breath
            loc0=pb-Buffer;
            loc1=cb+Buffer;
            timeSection = timeSeries(loc0:loc1);% get section of time stamps
            respSection = movmean(streams{idxRespiration}.time_series(loc0:loc1),10);
            [npeaks, peakLocs] = findpeaks(respSection,'MinPeakDistance',RT/200); % initial peak
        
            respSection = interp1(timeSection,respSection, timeSection(1):0.02:timeSection(end));
            timeSection = timeSection(1):0.02:timeSection(end);
        
            ft = fittype('gauss8');  % Gaussian curve fitting
            [coeffs, gof,output] = fit(timeSection', double(respSection)', ft);
            obj = @(x) filloutliers(coeffs(x),'spline');
            derSection = differentiate(coeffs,timeSection);
            [npeaks, peakLocs] = findpeaks(coeffs(timeSection),'MinPeakDistance',RT/200,'MinPeakHeight',mean(respSection)); % number of peaks in signal seperated by at least 1/2 breath
            
            [~,~,zeroCross]=zerocrossrate(derSection);
            % Get rid of first and last points 
            if zeroCross(1)==1
                zeroCross(1)=0;
            end
            if zeroCross(end)==1
                zeroCross(end)=0;
            end
        
            % Assigning zero crossing to clusters within 1 second of each other
            zero_crossing_idx = find(zeroCross);
            cluster_idx = ones(size(zero_crossing_idx));
            cluster_id = 1;
            for i = 2:length(zero_crossing_idx)
                if zero_crossing_idx(i) - zero_crossing_idx(i-1) > 20
                    cluster_id = cluster_id + 1;
                end
                cluster_idx(i) = cluster_id;
            end
            % Find clusters of points and takes average loc
            unique_clusters = unique(cluster_idx);
            closest_zero_idx = zeros(length(unique_clusters), 1);
            for i = 1:length(unique_clusters)
                cluster_members = zero_crossing_idx(cluster_idx == unique_clusters(i));
                closest_zero_idx(i) = round(mean(cluster_members));
            end
            [~,is_peak]=findpeaks([zeros(1,10) respSection(closest_zero_idx) zeros(1,10)]);
            is_peak=is_peak-10;
            isPeak = false(size(closest_zero_idx));
            isPeak(is_peak) = true;
            % t0:t1 range of current trial in timestamps
        
            %%% Find the p2 p3 that are within t0 and t1
            peakStamps = timeSection(peakLocs);
            trialPeaks = find(peakStamps < t1 & peakStamps > t0 & npeaks'>0.8*mean(npeaks));
            peak2 = peakStamps(trialPeaks(1));
            peak3 = peakStamps(trialPeaks(2));
            
            peakTimes = peakStamps(trialPeaks);
            if length(trialPeaks)~=2
%                 error('multiple peaks')
                toplot=1;
            end
            [~,closestLoc] = min(abs(timeSection(closest_zero_idx)'-peakTimes)); % time stamps of peaks and troughs from dy
            inhalPkidx = closest_zero_idx(closestLoc); % Peak of inhalation during trial
        
            %%% Find flanker peaks p1 and p4
        
        %     trialPeaks = find((peakStamps > t1 | peakStamps < t0) & npeaks'>0.8*mean(npeaks));
            trialPeaks = [trialPeaks(1)-1 trialPeaks(2)+1];
        
            peak1 = peakStamps(trialPeaks(1));
            peak4 = peakStamps(trialPeaks(2));
            flankerPeaks = peakStamps(trialPeaks);
            if length(flankerPeaks)~=2
                toplot =1;
            end
            [~,closestLoc] = min(abs(timeSection(closest_zero_idx)'-flankerPeaks)); % time stamps of peaks and troughs from dy
            flankerPkidx = closest_zero_idx(closestLoc); % Peak of inhalation during trial
        
            %%% Finding troughs between flanker peaks
            pks4 = sort([inhalPkidx ;flankerPkidx]);
        
            troughIdx = closest_zero_idx(closest_zero_idx > flankerPkidx(1) & ...
                closest_zero_idx < flankerPkidx(2) & ~ismember(closest_zero_idx,inhalPkidx));
            troughIdx = average_troughs(troughIdx, pks4);
        
            %%% Finding troughs outisde of flanker breaths closest one to flanker
            d1 = troughIdx(1)-flankerPkidx(1); % distance from flanker to trough
            d2 = flankerPkidx(2)-troughIdx(2); % other distance
        
            searchLim1 = round(flankerPkidx(1) - 1.5*d1); % make symmetrical and go back 20% further
            searchLim2 = round(flankerPkidx(2) + 1.5*d2); % make other side symmetrical 
        
            lastPeak = find(isPeak,1,'last');
            while lastPeak - closestLoc(2) > 2
                lastPeak =lastPeak-2;
            end
            if closest_zero_idx(lastPeak) == flankerPkidx(2)
                lastPeak = inf;
            end
            firstPeak = find(isPeak,1,'first');
            while closestLoc(1) - firstPeak > 2
                firstPeak =firstPeak+2;
            end
        
            if closest_zero_idx(firstPeak) == flankerPkidx(1)
                firstPeak = -inf;
            end
        
            firstTrough = find(closest_zero_idx > searchLim1 & closest_zero_idx< flankerPkidx(1) & ~isPeak);
            firstTrough = closest_zero_idx(firstTrough(firstTrough>firstPeak));
            lastTrough = find(closest_zero_idx < searchLim2 & closest_zero_idx> flankerPkidx(2)& ~isPeak);
            lastTrough =closest_zero_idx(lastTrough(lastTrough<lastPeak));
            if length(firstTrough)~=1 | length(lastTrough)~=1
                errCd = 'outside bounds not found';
                toplot=1;
            end
            if toplot
                clf
                hold on
                % fitted data plots
                plot(timeSection,respSection,'LineWidth',3)
                plot(timeSection,coeffs(timeSection)) % raw fit
            
                plot(t1,streams{idxRespiration}.time_series(cb),'b*','markersize',10,'LineWidth',3) % end of trial
                plot(t0,streams{idxRespiration}.time_series(pb),'b*','markersize',10,'LineWidth',3)  % start of trial
            
                plot(timeSection(inhalPkidx),respSection(inhalPkidx),'color','#D95319','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
                plot(timeSection(flankerPkidx),respSection(flankerPkidx),'color','g','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
                plot(timeSection(troughIdx),respSection(troughIdx),'color','m','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
                plot(timeSection(lastTrough),respSection(lastTrough),'color','r','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
                plot(timeSection(firstTrough),respSection(firstTrough),'color','r','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
                waitforbuttonpress
            end
        catch
            discardedTrials=discardedTrials+1;
        end
    end
    discardedTrials
end

%% Phase anlaysis
figure
[~,loc]=min(abs(timeSeries'-pressedTime));

slowAngles=respirationPhase(loc(slowPresses));
fastAngles=respirationPhase(loc(attendedPresses));

[pval table] = circ_wwtest(fastAngles,slowAngles); % ttest between slow and fast angles
avgFast = circ_mean(fastAngles');
avgSlow = circ_mean(slowAngles');

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