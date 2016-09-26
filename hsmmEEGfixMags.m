function [lkh1 params1 eventprobs1]=hsmmEEGfixMags(zscores,mags,params,thresh,x,y)
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
    params1 =params;
    eventprobs1=eventprobs;
    params=gammaParams(eventprobs,lens,2);
    for i = 1:nstates+1
        if prod(params(i,:)) < 5
            params(i,:)=params1(i,:);
        end
    end
    [lkh eventprobs]=calcEEG50h(bumps,mags,params,x,y);
    %params(:,2)'
    %lkh-lkh1
    %lkh
end
    
    
