function [meantimes durations times]=estimateTimes(eventprobs,subjects,conds,x,y)
lens=y-x+1;
nsubjs=max(subjects);
subjects=subjects(x);
ncases=length(x);
nconds=max(conds);
nstates=size(eventprobs,3);
times=zeros(ncases,nstates);
for j =1:ncases
        times(j,:)=[1:lens(j)]*reshape(eventprobs(1:lens(j),j,:),lens(j),nstates);
end
times=cat(2,times,lens);
durations=zeros(ncases,nstates+1);
durations(:,1)=times(:,1);
durations(:,2:end)=times(:,2:end)-times(:,1:end-1);
durations(:,nstates+1)=durations(:,nstates+1)+2;
durations(:,1)=durations(:,1)-2;
meantimes=zeros(nconds,nstates+1,nsubjs);
for j = 1:nsubjs
       for k = 1:nconds
           meantimes(k,:,j)=mean(durations((subjects==j).*(conds==k)==1,:));
       end
end

    
        
    
    
