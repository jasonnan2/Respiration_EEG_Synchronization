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
% id: 18 - error isn't normal behavior 
subIds = [1,3,5,6,7,9,10,11,13,14,15,16,17,18,19,20];
count=0;

id=3
rawFile = "A:\TwoTap\twotap"+string(id)+"\twotap"+string(id)+".xdf";
behFile = pickfiles({"A:\TwoTap\twotap"+string(id)}, {"trial_summary.csv"});
codeFile = pickfiles({"A:\TwoTap\twotap"+string(id)}, {".csv"},{".csv"},{'trial_summary'});

behTbl = readtable(behFile);
codeTbl = readtable(codeFile);
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
respiration=streams{idxRespiration}.time_series; % respiration is sampled at 10Hz
timeSeries=streams{idxRespiration}.time_stamps; 

markers=streams{idxMarker}.time_series;
markersTime=streams{idxMarker}.time_stamps;
% if size(behTbl,1)~= sum(markers==11)
%     count=count+1;
% end



%% Debugging behTable file
clc
clear responseTbl
if size(behTbl,1)~= sum(markers==11) % Sometimes behTbl will have an extra marker at start
    nDiff=size(behTbl,1)-sum(markers==11);
    responseTbl = codeTbl(codeTbl.code==11,:); % response generated from codes table
    responseTbl{end+1:end+nDiff,:}=nan(nDiff,4);
    responseTbl.behTblDuration=behTbl.ResponseTime;
    responseTbl.behTbllocalTime=behTbl.Timestamp;
    responseTbl.behTblAccuracy=behTbl.Accuracy;
    responseTbl = responseTbl(:,[1,7,5,2,4,6]);
end


%% Fixing behTbl

responseTbl = codeTbl(codeTbl.code==11,:); % response generated from codes table

for i=1:size(responseTbl,1) % go through all codes in responseTbl
    while abs(responseTbl.duration(i) - behTbl.ResponseTime(i)) >5 % buffer of 5mS
        % as long as the trial in behTbl is not the same as the response,
        % remove that row. "i" stays the same so it will check the next row once the
        % previous one is removed
        behTbl(i,:)=[];
    end
end

%% Pre processing EEG
[EEGset, ProcessVars] = preProcessing_EEGO_twotap(rawFile, behFile);
% for t=1:length(pressedTime)
%    [~,loc]= min(abs(pressedTime(t)-timeSeries)); % 10Hz sample rate
%     epochedResp(:,t)=respiration(loc-10:loc+5);
% end
slowTrialsNeural  = ProcessVars.behtabledata.Accuracy==1  & ismember(1:size(ProcessVars.behtabledata,1),ProcessVars.trialset(2).ind)';
fastTrialsNeural  = ProcessVars.behtabledata.Accuracy==1  & ismember(1:size(ProcessVars.behtabledata,1),ProcessVars.trialset(1).ind)';

%%
respiration = highpass(respiration,0.01,10); % cutoff at 0.01Hz = lets in anything under 100 seconds/ breath cycle 
respiration = movmean(respiration,5); % Very simple low pass filter

dydx = gradient(respiration(:)) ./ gradient(timeSeries(:)); % calculate derivitive 
respirationPhase = angle(hilbert(respiration)); % calculate phase 

% getting markers from xdf


%%
pressesLoc = find(markers==11); % all button presses
pressedTime=markersTime(pressesLoc);

close all
subplot(211)
hold on
plot(timeSeries,(respiration))
[m,l]=min(abs(timeSeries-pressedTime')');
plot(pressedTime,respiration(l),'r*')
ylabel('Respiration Signal')

subplot(212)
hold on
plot(timeSeries,(respirationPhase))
[m,l]=min(abs(timeSeries-pressedTime')');
plot(pressedTime,respirationPhase(l),'r*')
ylabel('Respiration Phase')

plot(timeSeries,zeros(size(timeSeries)),'k--')

%%



%%
figure
[~,loc]=min(abs(timeSeries'-pressedTime));

slowAngles=respirationPhase(loc(slowTrials));
fastAngles=respirationPhase(loc(fastTrials));
histogram(slowAngles,15)
hold on
histogram(fastAngles,15)
legend('slow','fast')
xlabel('phase angle')
ylabel('Count')

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
title('Angles of Slow and Fast Trials');
grid on;
set(gca,'rticklabels','')
rlim([0 1.01])
legend('slow','fast')
