function [lkh1 mags1 eventprobs1]=hsmmEEGfixParams(zscores,mags,params,thresh,x,y)
bumps=calcBumps(zscores);
lkh1=-Inf;
nstates=size(mags,2);
ncases=length(x);
ndims=size(bumps,2);
lens=y-x+1;
means=zeros(300,ncases,ndims);
for i = 1:ncases
    means(1:lens(i),i,:)=bumps(x(i):y(i),:);
end
[lkh eventprobs]=calcEEG50h(bumps,mags,params,x,y);
while (lkh-lkh1)>thresh
    lkh1=lkh;
    mags1 = mags;
    params1 =params;
    eventprobs1=eventprobs;
    for i = 1:nstates
        for j = 1:ndims
            mags(j,i)=mean(sum(eventprobs(:,:,i).*means(:,:,j)));
        end
    end
    [lkh eventprobs]=calcEEG50h(bumps,mags,params,x,y);
    %params(:,2)'
    %lkh-lkh1
    %lkh
end
    
    
