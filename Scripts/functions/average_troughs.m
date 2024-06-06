function troughIdx_avg = average_troughs(troughIdx, pks4)
    for i=1:length(pks4)-1

        between = find(troughIdx > pks4(i) & troughIdx< pks4(i+1));
        if length(between)==1
            troughIdx_avg(i)=troughIdx(between);
        else
            troughIdx_avg(i)=round(mean(troughIdx(between)));
        end
    end
end