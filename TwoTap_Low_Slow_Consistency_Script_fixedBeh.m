%% This step extract individul subjects from .set format and save them into .mat files format
Put_EEGLAB_inPath = addpath ('/eeglab');
eeglab

clc
clear

clustypeNames = {'thetaClus','alphaClus','betaClus'}; % To be used in Frequency band specific Source Mapping
labelnames = {'theta','alpha','beta'}; % To be used in Frequency band specific Source Mapping

%%% src localization things for neatlabs1
%src localization
solverType = 'bsbl';
saveFull = true;
account4artifacts = true;

src2roiReductionType = 'mpower';
templatefile = '\eeglab\plugins\BrainEAnalysis\headModel_templateFile_32chan.mat'; %headModel.getDefaultTemplateFilename();
%         templatefile = headModel.getDefaultTemplateFilename();
conductivity = [0.33, 0.022, 0.33];
orientation = true;
windowSize = (40/1000)*250; %EEG.srate;

hm = headModel.loadFromFile('\eeglab\plugins\BrainEAnalysis\headModel_templateFile_32chan.mat');
T = hm.indices4Structure(hm.atlas.label);
T = double(T)';
P = sparse(bsxfun(@rdivide,T, sum(T,2)))';

EEGdata = [];
EEGtimes = [];
EEGbaselinedata = [];
EEGbasetimes = [];

%% Initialization of data files
files_all = [];
files = cell(1);

Finalresults = '\Completed Projects\TwoTap_Project\TwoTap_Results\Source_TT_Results_May2024';
IndividualResults = '\Completed Projects\TwoTap_Project\TwoTap_Results\Source_TT_Results_May2024';

taskNames = {'go_green','middle_fish','lost_star','lucky_door','face_off','two_tap'};
Ntasks = length(taskNames);

cndnNames = {'stim'};
cndn = 1;


%% Envelope extraction

for taski = [6] %  % Can run for all the task in one go

    %pick the appropriate peb+ cleaned files
    rawData = 'D:\Completed Projects\TwoTap_Project\TwoTap_Raw_Data\Input_Files_for_Preprocessing';
    results = 'D:\Completed Projects\TwoTap_Project\TwoTap_Results\All_Preprocessedfiles_TwoTap';
    files2 = pickfiles(rawData, {'.xdf' 'pre'});
    files_all = char(squeeze(files2));

    for subii = 1:size(files_all,1)
        %extract the subject id here / pilot id or the pilot name here such as
        tabledata2 = [];

        rawFile1 = deblank(files_all);
        rawFile = deblank(rawFile1(subii,:));
        [pathName, rawFileName] = fileparts(rawFile);
        Save_Particp_ID = {rawFileName}; % Save this variable in trial info excelsheet
        Particp_ID (subii,1) = Save_Particp_ID;

        disp(['Processing subject ' rawFileName ' ' num2str(subii) '/' num2str(size(files_all,1))]);

        try

            Sessionname = 'pre'; %% will have to change depending on study type (i.e. CoCo will change to post/pre or ucmed
            subjectindex = Particp_ID(subii,1);

        end
        meantscoreAct = nan(68,7);
        taskset_invalid = [];

        for taski= 6 %1:Ntasks

            tabledata_hi = [];
            theta = nan(68,1);
            alpha = nan(68,1);
            beta = nan(68,1);

            try
                tscoreAct = nan(3,68);

                task = taski;
                block = 1;

                if task == 6
                    pcfile.tfinfo{1} = [-4000 4000 3 7];
                    pcfile.tfinfo{2} = [-4000 4000 8 12];
                    pcfile.tfinfo{3} = [-4000 4000 13 30];

                end

                disp(['Processing task' num2str(task)]);

                subjectTaskFolder = fullfile(results,rawFileName,taskNames{taski});
                behFile = pickfiles(pathName, {rawFileName Sessionname 'twotap' 'trial_summary.csv'});

                EEGmain = struct;
                EEGmain_baseline = struct;

                specName = 'tmapStr_newpeak_flank';
                if taski == 1
                    specName1 = 'nogo1_';
                elseif taski == 5
                    specName1 = '.';
                elseif taski == 4
                    specName1 = '_ch_';
                else
                    specName1 = '.';
                end

                %pipeline settings, block info is present
                processvarFiles = pickfiles(subjectTaskFolder, {'processvars', specName1, taskNames{task}, rawFileName, '_stim.mat'});
                load(processvarFiles,'pipelineSettings');

                %%% Fixing the twotap behavior trial summary
               codeFile = replace(behFile,'_trial_summary','');
                codeTbl = readtable(codeFile); % origianl lsl csv
                behTbl = readtable(behFile); % original trial summary

                responseTbl = codeTbl(codeTbl.code==11,:); % response generated from codes table
                diffRows = abs(size(responseTbl,1)-size(behTbl,1));
                count=0;
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
                end
                if count ~= diffRows
                    error('Check behTbl')
                end

                % based on all trials
                rtmedian = median(behTbl.ResponseTime);
                rtstd = std(behTbl.ResponseTime);
                madRT = mad(behTbl.ResponseTime);

                % these are ok, they are just 1:#trials
                trialsetind = pipelineSettings.trialsetind{block};
                blockind = (pipelineSettings.blockind{block});
                Orig_Trial_Numb_Aft_Preprocs = size(trialsetind,1); % Save this for each subject

                % On focusing on good neural trials now
                behTbl(pipelineSettings.rejTrials,:)=[]; % Remove bad neural trials
                RTall = behTbl.ResponseTime;
                accall = behTbl.Accuracy;
                adtrial = abs(RTall- rtmedian);


                %%% Find all low trials > 1 mad of all trial RT
                trialsind{2} = (trialsetind(blockind(find(adtrial(blockind) > madRT)))); % For low consistency, trials outside of 1 median deviation
                % But still avoiding extreme outliers

                low = trialsind{2};

                %%% Find all low Slow trials using following code
                trialsind{1} = low(find(RTall(low) > (rtmedian - madRT))) % New Version Code for Low_Slow trials

                faci_tilda = 1;
                trialno = size(trialsind{1},1);
                Excluded_Trial_nos = size(trialsind{2},1) - size(trialsind{1},1)
                Percntg_Excluded_trials = (Excluded_Trial_nos./Orig_Trial_Numb_Aft_Preprocs)*100

                LowSlowConsis_trialnos =  trialno;
                Percntg_LowSlowCons_trials = (LowSlowConsis_trialnos./Orig_Trial_Numb_Aft_Preprocs)*100

                pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_', specName1 '.set',cndnNames{1}});
                EEGmain = pop_loadset(deblank(pebFiles));

                if EEGmain.trials ~=size(behTbl,1)
                    error('Checl behTbl')
                end

                filename = ['brainmap_'];

                for clustype = 1:3 %pos, neg


                    EEG = EEGmain;

                    tfinfoarr = pcfile.tfinfo{clustype};

                    if isempty(tfinfoarr) == 0
                        timesind1 = find(EEG.times > tfinfoarr(1));
                        timesind1 = timesind1(1);
                        timesind2 = find(EEG.times < tfinfoarr(2));
                        timesind2 = timesind2(end);

                        %adding the following to save the pop
                        %inverse solution
                        if (timesind2 - timesind1) < 80
                            timesind2 = timesind1 + 80;
                        end

                        freqind1 = tfinfoarr(3); freqind2 = tfinfoarr(4);

                        fileEEGsrc =  pickfiles(results, {filename});
                        %                 if isempty(fileEEGsrc)

                        EEG = pop_eegfiltnew(EEGmain, freqind1, freqind2, 826,0,[],0);

                        EEG.data = permute([permute(EEG.data(:,timesind1:timesind2,trialsind{1}(trialsind{1} < size(EEG.data,3))),[3 1 2]).*1],[2 3 1]);
                        EEG.data = nanmean(EEG.data,3);

                        EEG.times = EEGmain.times(timesind1:timesind2);
                        EEG.trials = size(EEG.data,3); %length(trialsind{facii});
                        EEG.pnts = length(EEG.times);

                        EEG.epoch = [];
                        EEG.etc.subjectNames = subjectindex;

                        EEG = pop_forwardModel(EEG, templatefile, conductivity, orientation);
                        EEG = pop_pebp(EEG,saveFull, account4artifacts, src2roiReductionType);


                        % Need flanks to avoid edge artifacts while taking envelope (each side flank is of 500 msec in length(125 samples) in 256 sampling rate
                        % Final Section genrates envelopes of each frequency bands
                        %need flanks to avoid edge
                        %artifacts while taking
                        %envelope (each side flank
                        %is of 500 msec in
                        %length(125 samples) in 256
                        %sampling rate
                        flanklen = 125;
                        %                                                     flankact = reshape(randsample(EEG.etc.src.act(:),size(EEG.etc.src.act,1)*(flanklen)),[size(EEG.etc.src.act,1), flanklen]);
                        %                                                     flankactFull = reshape(randsample(EEG.etc.src.actFull(:),size(EEG.etc.src.actFull,1)*(flanklen)),[size(EEG.etc.src.actFull,1), flanklen]);
                        %                                                 flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);

                        EEG1 = EEG; EEG2 = [];
                        flanki_s = randi(size(EEG.etc.src.act,2)*size(EEG.etc.src.act,3)-flanklen-1,[68,1]); flankact = EEG.etc.src.act(:,flanki_s:flanki_s+flanklen-1);

                        %                         if subi == 1
                        EEG1.etc.src.act(:,1:end+flanklen*2,1) = envelope(cat(2,flankact,EEG.etc.src.act(:,:,1),flankact)',1,'peak')';
                        %                         else
                        %                             EEG1.etc.src.act(:,:,subi) = envelope(cat(2,flankact,EEG.etc.src.act(:,1:end,subi),flankact)',1,'peak')';
                        %                         end
                        %restoring the time
                        %length
                        tempvar = EEG1.etc.src.act(:,flanklen+1:end-flanklen,1);
                        EEG2.etc.src.act(:,:,1) = tempvar;

                        EEG.etc.src.act = EEG2.etc.src.act; clear EEG1 EEG2;

                        if size(tempvar,2) == 1999
                            tabledata(clustype,:,:) = tempvar;
                            test(subii) = 0;
                        else
                            tabledata(clustype,:,:) = nan(68,1999);
                            test(subii) = 1;

                        end

                        %record all ROI for each clustype

                        if clustype == 1
                            theta = tempvar;
                        elseif clustype == 2
                            alpha = tempvar;
                        else
                            beta = tempvar;
                        end

                    end
                end
            catch
                warning(deblank(pebFiles)+" failed")


            end
            % Per subject source localized results

            save(fullfile(IndividualResults,[rawFileName, '_TwoTap_src_LowSlowConsis.mat']),'tabledata');


            tabledata2(:,:,:,taski) = tabledata;

        end

        % Collation of all subjects in the group
        if ~isempty(tabledata2)

            Pre_TT_srcemap(:,:,:,:,subii) = tabledata2;
            Pre_TT_srcmap_LowSlowConsis = squeeze(Pre_TT_srcemap(:,:,:,6,:)); % Get rid of the task dimension

        end

        filename_Trials = [Finalresults '/Trial_Numbers_LowSlowConsis_Info_CoCo.xlsx'];
        export_trials = cat(2,Orig_Trial_Numb_Aft_Preprocs, LowSlowConsis_trialnos, Percntg_LowSlowCons_trials);
        temp = array2table(export_trials, 'VariableNames',{'Orig_Trial_Numb_Aft_Preprocs', 'LowSlowConsis_trialnos','Percntg_LowSlowConsis_trials'});
        temp2(subii,:)= temp;

    end

    temp2= addvars (temp2, Particp_ID, 'After', 'Percntg_LowSlowConsis_trials')


    writetable(temp2,filename_Trials,'Sheet','Pre_Trials');
    save(fullfile(Finalresults,['Pre_TwoTap_src_LowSlowConsis.mat']),'Pre_TT_srcmap_LowSlowConsis');

end


% size(Post_TT_srcmap_Low_SlowConsis)
% size(Pre_TT_srcmap_Low_SlowConsis)
