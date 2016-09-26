function [average result maxes]=averagefunction10(data,eventprobs,x,y)
ncases=length(x);
lens=y-x+1;
nstates=size(eventprobs,3);
avlen=round(mean(lens));
maxes=zeros(ncases,nstates);
ndims=size(data,2);
result=zeros(avlen,ncases,ndims);
for i = 1:nstates
    [~,maxes(:,i)]=max(eventprobs(:,:,i));
end
maxes=cat(2,-1*ones(ncases,1),maxes,lens);
indices=round(mean(maxes));
for i = 1:ncases
    if min(maxes(i,2:end-1)-maxes(i,1:end-2))>1
        temp=data(x(i):y(i),:);
        for j = 1:nstates
            if maxes(i,j+1)-maxes(i,j)>4
                result(indices(j)+2:indices(j+1)-2,i,:)=warp(temp(maxes(i,j)+2:maxes(i,j+1)-2,:),indices(j+1)-indices(j)-3);
                result(indices(j)+2:indices(j+1)-2,i,:)=warp(temp(maxes(i,j)+2:maxes(i,j+1)-2,:),indices(j+1)-indices(j)-3);
            end
        end
        for j = 1:nstates
            result(indices(j+1)-2:indices(j+1)+2,i,:)=temp(maxes(i,j+1)-2:maxes(i,j+1)+2,:);
        end
        if (maxes(i,end)-maxes(i,end-1))>2  && (indices(end)-indices(end-1))>2
            result(indices(end-1)+3:end,i,:)=warp(temp(maxes(i,j+1)+3:end,:),indices(end)-indices(end-1)-2);
        end
    end
end
average=reshape(mean(result(:,sum(result(:,:,1)==0)==0,:),2),indices(end),ndims);
        