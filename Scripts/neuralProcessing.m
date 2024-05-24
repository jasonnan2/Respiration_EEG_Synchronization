%% Twotap revision s with 64 chan
% Jason Nan

templatefile='headModel_templateFile_63chan.mat';
windowSize = (40/1000)*500; % 40 ms
overlaping = 25;
solverType = 'bsbl';
saveFull = false;
account4artifacts = true;
postprocCallback = [];
src2roiReductionType = 'hist';
conductivity = [0.33, 0.022, 0.33];
orientation = false;
%%
task='twotap';
results='A:\TwoTap_results\';
subIds = [1,3,5,6,7,9,10,11,13,14,15,16,17,18,19,20];

for s=subIds
    resultsFolder = [results  task num2str(s)];

    if ~ exist(resultsFolder,"dir")
        mkdir(resultsFolder)
    end

    rawFile = char("A:\TwoTap\twotap"+string(s)+"\twotap"+string(s)+".xdf");
    behFile = pickfiles({"A:\TwoTap\twotap"+string(s)}, {"trial_summary.csv"});

    [EEGset, ProcessVars] = preProcessing_EEGO_twotap(rawFile,behFile)
    k=1;
    EEGset{k}.etc.pipelineSettingsFile = fullfile('',['processvars_twotap_' task '_' rawFile '_' EEGset{k}.setname,'.mat']);
    pop_saveset(EEGset{k},'filepath',resultsFolder,'filename',['epoch_twotap_' task num2str(s) '_' EEGset{k}.setname],'savemode','onefile');
    save([resultsFolder '\ProcessVars.mat'],'ProcessVars')
    
    epochFiles = pickfiles(resultsFolder,{'epoch_twotap_','.set'},{'epoch_twotap_','.set'},{'peb'});
    nepochFiles = size(epochFiles,1);
    [~, fileName] = fileparts(deblank(epochFiles(k,:)));
    EEG = pop_loadset(deblank(epochFiles(k,:)));
    windowSize = (40/1000)*EEG.srate; % 40 ms
    % Run PEB and save cleaned data
    EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
    EEG = pop_inverseSolution(EEG,windowSize, overlaping, solverType, saveFull, account4artifacts, src2roiReductionType);
    pop_saveset(EEG,'filepath',resultsFolder,'filename',['Newpeb_cleaned_' fileName],'savemode','onefile');
end
