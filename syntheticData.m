function [signal x y pos pure LP] = syntheticData(conds,mags,params, SNR ,meanpower)
%take mags and params from a known fit
    nconds=max(conds);
    ntrials=length(conds);
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
    LP=zeros(300,nstates+1,nconds);
    for j = 1:nstates+1
        for i = 1:nconds
            if j==nstates+1 || params(j,2,i)~=0 
                LP(:,j,i)=gammaEEG(params(j,1,i),params(j,2,i),300);
            end
        end
    end
    LP=cumsum(LP);
    index=1;
    j = 1;
    while j <= ntrials 
        cond = conds(j);
        previous=index;
        for t = 1:nstates
            if params(t,2,cond) ~= 0
                next = find(LP(:,t,cond)>rand(),1);
                pos(j,t) = index+next-previous+2;
                signal(index+next-1:index+next+3,:) = bumps(:,:,t);
                index=index+next+4;
            end
        end
        index=index+find(LP(:,nstates+1,cond)>rand(),1);
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
    signal = signal(1:index-1,:);
    pure = signal;
    for k = 1:ntrials
        Noi = noise3(lens(k), 100,10,ndims,meanpower);
        Sig = signal(x(k):y(k),:);
        signal(x(k):y(k),:) = zscore(Sig + Noi * (bandpower(reshape(Sig,1,size(Sig,1)*size(Sig,2))) / bandpower(reshape(Noi,1,size(Noi,1)*size(Noi,2))))/SNR);
    end
end


function [signal] = noise3(frames, srate,n,ndims,meanpower)
    signal = zeros (frames,ndims);
    for j = 1:ndims
       freq=0;
       for i = 1:n
          freq = min(125,freq + (4*rand(1)));
          freqamp = meanpower(ceil(freq)) / meanpower(1);
          phase = rand(1)*2*pi;
          signal(:,j)  = signal(:,j) + sin([1:frames]'/srate*2*pi*freq + phase) * freqamp;
       end
    end
end