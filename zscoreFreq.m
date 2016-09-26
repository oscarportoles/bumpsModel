function data =zscoreFreq(data,subjects)
nsubjs=max(subjects);
for i = 1:nsubjs;
    a=find(subjects==i);
    data(a,:)=zscore(data(a,:));
end

