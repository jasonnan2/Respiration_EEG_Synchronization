%% Copy all Resting state data over to data3

taskFolderName = 'rest';
destDir = '/media/owner/data3/RestingState';

% Get a list of all subject folders in the project directory
subjectFolders = dir(fullfile(results, '*'));
subjectFolders = subjectFolders([subjectFolders.isdir]);  % keep only directories
subjectFolders = subjectFolders(~ismember({subjectFolders.name}, {'.', '..'}));  % remove '.' and '..'

% Loop over all subject folders
for i = 1:length(subjectFolders)
    % Construct the full path to the task folder
    taskFolder = fullfile(results, subjectFolders(i).name, taskFolderName);
    
    % Check if the task folder exists
    if isfolder(taskFolder)
        % Construct the full path to the destination folder
        destFolder = fullfile(destDir, subjectFolders(i).name, taskFolderName);
        
        % Create the destination folder if it does not exist
        if ~isfolder(destFolder)
            mkdir(destFolder);
        end
        
        % Copy the task folder to the destination folder
        copyfile(taskFolder, destFolder);
    end
end
%%
results = '/media/owner/data3/RestingState';
saveDir = "/media/owner/data3/Jason/Active/Resting/data/";

allFolders = dir(results);
allFolders = allFolders([allFolders.isdir]);  % keep only directories
allFolders = allFolders(~ismember({allFolders.name}, {'.', '..'}));  % remove '.' and '..'
newtbl=table();
allSubs={}; count=0;
for i=1:length(allFolders)
    folderName = allFolders(i).name;
    subName = strsplit(folderName,'_');
    subName=subName{1};
    subjectFolder=fullfile(results, folderName);
    pebFiles = pickfiles(subjectFolder,{'peb_cleaned_','rest', 'na.set'});
    if ~isempty(pebFiles)
        count=count+1;
        allSubs{count}=subName;
        EEG=pop_loadset(pebFiles);
        EEG.etc.src=[];
        save(saveDir+subName+"_rest.mat",'EEG');
    end
    % Removing missing data
%     if isempty(pebFiles)
%         rmdir(subjectFolder,'s');
%     end
    %subIdx = find(strcmp(lower(tbl.SubjectID),lower(subName)));
    %newtbl=[newtbl ; tbl(subIdx,:)];
end



%%

results = '/media/owner/data3/RestingState';
tbl = readtable('Total_Data_324_for_Correlation.xlsx');

i=1;
%for i=1:size(tbl,1)

subjectname=tbl.SubjectID{i};
subjectTaskFolder = fullfile(results,subjectname,'rest');

pebFiles = pickfiles(subjectTaskFolder,{'peb_cleaned_','rest', subjectname, 'na.set'});
%pop_loadset(pebFiles)
%%
function allFolders = getAllFolders(dirName)
    % Get the data for the current directory
    dirData = dir(dirName);
    
    % Get the indices for directories
    dirIndex = [dirData.isdir];
    
    % Get a list of subdirectories
    subDirs = {dirData(dirIndex).name};
    validIndex = ~ismember(subDirs,{'.','..'});
    
    % Prepend path to subdirectories
    allFolders = cellfun(@(x) fullfile(dirName,x),...
        subDirs(validIndex),'UniformOutput',false);
    
    % Loop over valid subdirectories
    for iDir = 1:numel(allFolders)
        nextDir = allFolders{iDir};
        allFolders = [allFolders; getAllFolders(nextDir)];
    end
end
