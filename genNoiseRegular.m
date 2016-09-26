function [signal prefix pos pure LP] = genNoiseRegular(ntrials,mags,params, SNR )
%take mags and params from a known fit
    template=[0.3090    0.8090    1.0000    0.8090    0.3090]';
    ndims=size(mags,1);
    d = size(mags,2);
    nstates=size(mags,2);
    bumps=zeros(5,ndims,nstates);
    for i = 1:nstates
        bumps(:,:,i)=template*mags(:,i)';
    end
    nobs = ntrials*300;
    signal = zeros(nobs,ndims);
    prefix=zeros(nobs,1);
    pos = zeros(ntrials,nstates);
    LP=zeros(300,nstates+1);
    for j = 1:nstates+1
            if j==nstates+1 || params(j,2)~=0 
                LP(:,j)=gammaEEG(params(j,1),params(j,2),300);
            end
    end
    LP=cumsum(LP);
    index=1;
    j = 1;
    while j <= ntrials 
        previous=index;
        for t = 1:nstates
            if params(t,2) ~= 0
                next = find(LP(:,t)>rand(),1);
                pos(j,t) = index+next-previous+2;
                signal(index+next-1:index+next+3,:) = bumps(:,:,t);
                index=index+next+4;
            end
        end
        index=index+round(gamrnd(params(nstates+1,1),params(nstates+1,2)));
        if index-previous>300 || index-previous < 40
            signal(previous:index,:)=0;
            index=previous;
        else
            prefix(previous:index-1) = [1:index-previous]'; 
            j = j+1;
        end
    end  
    x=find(prefix(:,1)==1);
    y=cat(1,x(2:end,1),index)-1;
    lens=y-x+1;
    prefix = repmat(prefix(1:index-1),1,4);
    signal = signal(1:index-1,:);
    pure = signal;
    load meanpower
    for k = 1:ntrials
        Noi = noise3(lens(k), 100,10,ndims,meanpower);
        Sig = signal(x(k):y(k),:);
        signal(x(k):y(k),:) = zscore(Sig + Noi * (bandpower(reshape(Sig,1,size(Sig,1)*size(Sig,2))) / bandpower(reshape(Noi,1,size(Noi,1)*size(Noi,2))))/SNR);
    end
end