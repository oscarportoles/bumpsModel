function vals = warp(data,len)
ndims=size(data,2);
vals=zeros(len,ndims);
n = size(data,1);
data=reshape(repmat(reshape(data,1,n,ndims),[len 1 1]),n,len,ndims);
for i = 1:len
    vals(i,:)=mean(data(:,i,:),1);
end