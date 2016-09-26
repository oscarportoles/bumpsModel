function [params events temp] = gammaParamsCond(eventprobs,lens,conds,map)
offset1=3;
spacing=5;
nconds=max(conds);
map=cat(2,map,ones(nconds,1));
nstates=size(eventprobs,3);
temp=zeros(300,nstates,nconds);
params=zeros(nstates+1,2,nconds);
for i = 1:nconds
    temp(:,:,i)=reshape(mean(eventprobs(:,conds==i,:),2),300,nstates);
end
maxobs=301-(2*(offset1-1)+(nstates-1)*spacing);
events=zeros(nstates+2,nconds);
for j = 1:nconds
    offset=-2;
    for i = 1:nstates
        if map(j,i)~=0
            offset=offset+spacing;
            temp(1:maxobs,i,j)=temp(offset:offset+maxobs-1,i,j);
            events(i+1,j)=[1:maxobs]*temp(1:maxobs,i,j);
        else
            events(i+1,j)=events(i,j);
        end
    end
    events(i+2,j)=lens(j)-offset-1;
end
temp(maxobs+1:end,:,:)=[];
for j = 1:nstates+1
    uniques=unique(map(:,j));
    uniques=uniques(uniques>0);
    for i = 1:length(uniques)
        subs=find(map(:,j)==i);
        params(j,:,subs)=repmat([2,mean(events(j+1,subs)-events(j,subs))/2],[1,1,length(subs)]);
    end
end
params(1,2,:)=params(1,2,:)-.5/2;
params(2:nstates,2,:)=params(2:nstates,2,:)+.5/2;