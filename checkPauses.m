rawData = 'D:\Completed Projects\TwoTap_Project\TwoTap_Raw_Data';
files2 = pickfiles(rawData, {'trial_summary.csv' 'twotap','pre'});
files_all = char(squeeze(files2));
pauses=[];
for subii=1:size(files_all)
    behFile=deblank(files_all(subii,:));
    codeFile = replace(behFile,'trial_summary','');
    codeTbl=readtable(codeFile);
    if sum(codeTbl.code == 502 | codeTbl.code ==500502) & sum(codeTbl.code == 503 | codeTbl.code ==500503)
        pauses(end+1)=subii;
    end
%     
%     if sum(codeTbl.code == 503 | codeTbl.code ==500503)
%         behTbl
%         codeTbl
%         input(rawFile)
%     end
end
