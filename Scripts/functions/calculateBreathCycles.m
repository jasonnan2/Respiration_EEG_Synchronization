function [discardedTrials, avgI, avgE, avgBreathCycle,totalRate]=calculateBreathCycles(pressedTime,timeSeries,streams, idxRespiration, behTbl, trialNums, plotErr)

avgE=[];avgI=[];avgBreathCycle=[];
discardedTrials=0;
for t=1:length(trialNums)
    toplot=0;
    errcd=0;
    try
        RT = behTbl.ResponseTime(trialNums(t)); % get response time from behTbl
        t1=pressedTime(trialNums(t)); % time stamp of button press
        t0 = t1-RT/1000;% get previous button press
        [~,pb]=min(abs(timeSeries'-t0)); % index of previous button press relative to all time
        [~,cb]=min(abs(timeSeries'-t1)); % index of current button press relative to all time
    
        Buffer = round(RT/80); % go back half a breath
        loc0=pb-Buffer;
        loc1=cb+Buffer;
        timeSection = timeSeries(loc0:loc1);% get section of time stamps
        respSection = movmean(streams{idxRespiration}.time_series(loc0:loc1),5);
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
            errcd=1;
%             error('Too many peaks')
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
            errcd=1;
%             error('Flaker Peaks Wrong')
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
            errcd=1;
%             error('No outer bounds')
            toplot=1;
        end

        if errcd 
            error('Some Error')
        end
        

        %%% Caulating inhalation and expiration data
        taskTroughTimes = timeSection(troughIdx);
        taskPeakTimes = timeSection(inhalPkidx);
        flankerPeakTimes = timeSection(flankerPkidx);
        flankerTroughTimes = timeSection([firstTrough lastTrough]);

        I1 = taskPeakTimes(1)-taskTroughTimes(1);
        I2 = taskPeakTimes(2)-taskTroughTimes(2);
        E1 = taskTroughTimes(2)-taskPeakTimes(1);
        E2 = taskTroughTimes(3)-taskPeakTimes(2);
        
        avgE(end+1) = mean([E1,E2]);
        avgI(end+1) = mean([I1 I2]);
        %%%
        preBreathRange = [flankerTroughTimes(1) taskTroughTimes(1)];
        breath1Range = taskTroughTimes(1:2);
        breath2Range = taskTroughTimes(2:3);
        postBreathRange = [taskTroughTimes(3) flankerTroughTimes(2)];
        
        % time Stamps of trial
        press0 = t0;
        press1=t1;

        if press0>breath1Range(1) & press0<breath1Range(2)
            % within first breath
            breath1Cycle = (breath1Range(2)-press0)/diff(breath1Range);
        else
            breath1Cycle = 1+(preBreathRange(2)-press0)/diff(preBreathRange);
        end

        if press1>breath2Range(1) & press1<breath2Range(2)
            % within second breath
            breath2Cycle = (press1-breath2Range(1))/diff(breath2Range);
        else
            breath2Cycle = 1+(press1-postBreathRange(1))/diff(postBreathRange);
        end
        
        avgBreathCycle(end+1)=nansum([breath1Cycle breath2Cycle]);
    catch
        if plotErr
            clf
            hold on
            % fitted data plots
            plot(timeSection,respSection,'LineWidth',3)
            plot(timeSection,coeffs(timeSection)) % raw fit
        
            plot(t1,streams{idxRespiration}.time_series(cb),'b*','markersize',10,'LineWidth',3) % end of trial
            plot(t0,streams{idxRespiration}.time_series(pb),'b*','markersize',10,'LineWidth',3)  % start of trial
        
            plot(timeSection(inhalPkidx),respSection(inhalPkidx),'color','#D95319','marker','+','linestyle','none','markersize',20,'LineWidth',3) % peaks during trial
            plot(timeSection(flankerPkidx),respSection(flankerPkidx),'color','g','marker','+','linestyle','none','markersize',20,'LineWidth',3) % flanker peaks
            plot(timeSection(troughIdx),respSection(troughIdx),'color','m','marker','+','linestyle','none','markersize',20,'LineWidth',3) % troughs during trial
            plot(timeSection(lastTrough),respSection(lastTrough),'color','r','marker','+','linestyle','none','markersize',20,'LineWidth',3) % flanker trough
            plot(timeSection(firstTrough),respSection(firstTrough),'color','r','marker','+','linestyle','none','markersize',20,'LineWidth',3) % flanker trough
            waitforbuttonpress
        end
        discardedTrials=discardedTrials+1;
    end
end

totalRate = sum(avgBreathCycle)/(sum(avgI)+sum(avgE));

avgI=nanmean(avgI);
avgE= nanmean(avgE);
avgBreathCycle=nanmean(avgBreathCycle);



end
