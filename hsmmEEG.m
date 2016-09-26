function [lkh1 mags1 params1 eventprobs1]=hsmmEEG(zscores,mags,params,thresh,x,y)
% @ zsocres: zscored time series (10 PCA)
% @ mags: initial conditions for bumps magnitudes [nPCAs*nBumps]
% @ params: initial conditions for Gamma distribution scale parameter
% @ thresh: Thereshold or flag for the alfgorithm (1: do HSMM)
% @ x: begining of each trial
% @ y: end of each trial
% ---> output


bumps=calcBumps(zscores); 
% @ bumps: [Nsamples X nPCs]
% correlation at each sample of the sample and the
% previous 5 samples on time domain with a bump morphology.
lkh1=-Inf;
nstates=size(mags,2);   % number of bumps
ncases=length(x);       % number of trials
ndims=size(bumps,2);    % number of PCs from PCA
lens=y-x+1;             % length of each trial + 1

[lkh eventprobs]=calcEEG50h(bumps,mags,params,x,y);
% @ eventprobs: [samples(375)*ncases*nstates] = [375*trials*nBumps]
% % lkh: float, likelihood

if thresh==0
    lkh1=lkh;
    mags1 = mags;
    params1 =params;
    eventprobs1=eventprobs;
else
    means=zeros(375,ncases,ndims);
    for i = 1:ncases
        means(1:lens(i),i,:)=bumps(x(i):y(i),:); % arrange bumps dimensions by trials [375*trial*PCs]
    end
    while (lkh-lkh1)>thresh
        lkh1=lkh;
        mags1 = mags;
        params1 =params;
        eventprobs1=eventprobs;
        for i = 1:nstates
            for j = 1:ndims
                mags(j,i)=mean(sum(eventprobs(:,:,i).*means(:,:,j)));
                % mean across trials of the sum of all samples in a trial lenght as: 
                % (estimated bump prob * correlation of bump morphology and EEG)
                % repeated for each PC
            end
        end
        params=gammaParams(eventprobs,lens,2); % update gamma-2 parameter
        
        for i = 1:nstates+1
            if prod(params(i,:)) < 6 % originaly < 5
            % multiply scale and shape parameters to get the mean distance
            % of the gamma-2 pdf. I guess it constains that bumps are separated
            % at least a bump length
                params(i,:)=params1(i,:); 
            end
        end
        [lkh eventprobs]=calcEEG50h(bumps,mags,params,x,y); 
        %params(:,2)'
        %lkh-lkh1
        %lkh
    end
end
    
    
