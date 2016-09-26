function result = combineConds(averages,indices,positions,flag)
maxlen=max(indices(:,end));
nconds=length(averages);
npositions=length(positions);
result = zeros(maxlen,nconds,npositions);
for i = 1:nconds
    if flag == 1
        result(:,i,:)=cat(1,zeros(maxlen-indices(i,end),npositions),averages{i}(:,positions));
    else
        result(:,i,:)=cat(1,averages{i}(:,positions),zeros(maxlen-indices(i,end),npositions));
    end       
end