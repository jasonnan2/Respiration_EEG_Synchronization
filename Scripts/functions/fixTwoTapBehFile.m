function [behTbl, responseTbl]=fixTwoTapBehFile(behTbl,codeTbl)

    responseTbl = codeTbl(codeTbl.code==11,:); % response generated from codes table
    diffRows = abs(size(responseTbl,1)-size(behTbl,1));
    if size(responseTbl,1)~=size(behTbl,1)
        count=0;
        for i=1:size(responseTbl,1)-1% go through all codes in responseTbl
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
        % Fixing Nth trial
        if size(responseTbl,1)~=size(behTbl,1)
            % Finding the index where the trials match up
            lastDuration = responseTbl.duration(i+1); % last duration
            trailingTrials = behTbl.ResponseTime(i+1:end);
            [~,behAlignedTrial]=min(abs(trailingTrials-lastDuration));
            behTbl.ResponseTime(i+1) = sum(trailingTrials(1:behAlignedTrial));
            behTbl(i+2:end,:)=[]; % removing all the trailing trials
        end
    end
    % Checking if the number of RT differences is the same as number of
    % different rows initially
    if size(responseTbl,1)~=size(behTbl,1)
        error('Check behTbl')
    end
end

