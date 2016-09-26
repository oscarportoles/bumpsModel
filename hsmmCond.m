function [lkh1 mags1 params1 eventprobs1]=hsmmCond(zscores,conds,mags,params,thresh,x,y,map)
bumps=calcBumps(zscores);
nconds=max(conds);
if length(size(params))==2;
    params=repmat(params,[1,1,nconds]);
end
lkh1=-Inf;
nstates=size(mags,2);
ncases=length(x);
ndims=size(bumps,2);
lens=y-x+1;
means=zeros(300,ncases,ndims);
lengths=zeros(nconds,1);
for i = 1:nconds
    lengths(i)=mean(lens(conds==i));
end
for i = 1:ncases
    means(1:lens(i),i,:)=bumps(x(i):y(i),:);
end
[lkh eventprobs]=calcEEG50Condsh(bumps,conds,mags,params,x,y,map);
while (lkh-lkh1)>thresh
    lkh1=lkh;
    mags1 = mags;
    params1 =params;
    eventprobs1=eventprobs;
    for i = 1:nstates
        a=find(map(conds,i)>0);
        for j = 1:ndims
            mags(j,i)=mean(sum(eventprobs(:,a,i).*means(:,a,j)));
        end
    end
    params=gammaParamsCond(eventprobs,lengths,conds,map);
    for i = 1:nstates+1
        for j = 1:nconds
            if prod(params(i,:,j))< 5
                params(i,:,j)=params1(i,:,j);
            end
        end
    end
    [lkh eventprobs]=calcEEG50Condsh(bumps,conds,mags,params,x,y,map);
    %reshape(params(:,2,:),nstates+1,nconds)'
    %lkh-lkh1
    %lkh
end
    
    
