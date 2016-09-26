function [average result means]=warpMeans(data,eventprobs,x,y)
ncases=length(x);
lens=y-x+1;
nstates=size(eventprobs,3);
avlen=round(mean(lens));
means=zeros(ncases,nstates);
ndims=size(data,2);
result=zeros(avlen,ncases,ndims);
for i = 1:nstates
    means(:,i)=[1:300]*eventprobs(:,:,i);
end
means=cat(2,-1*ones(ncases,1),round(means),lens);
indices=round(mean(means));
for i = 1:ncases
    if min(means(i,2:end-1)-means(i,1:end-2))>1
        temp=data(x(i):y(i),:);
        for j = 1:nstates
            if means(i,j+1)-means(i,j)>4
                result(indices(j)+2:indices(j+1)-2,i,:)=warp(temp(means(i,j)+2:means(i,j+1)-2,:),indices(j+1)-indices(j)-3);
                result(indices(j)+2:indices(j+1)-2,i,:)=warp(temp(means(i,j)+2:means(i,j+1)-2,:),indices(j+1)-indices(j)-3);
            end
        end
        for j = 1:nstates
            result(indices(j+1)-2:indices(j+1)+2,i,:)=temp(means(i,j+1)-2:means(i,j+1)+2,:);
        end
        if (means(i,end)-means(i,end-1))>2  && (indices(end)-indices(end-1))>2
            result(indices(end-1)+3:end,i,:)=warp(temp(means(i,j+1)+3:end,:),indices(end)-indices(end-1)-2);
        end
    end
end
average=reshape(mean(result,2),indices(end),ndims);
        