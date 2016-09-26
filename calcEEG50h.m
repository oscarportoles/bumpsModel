function [likelihood eventprobs forward backward]=calcEEG50h(bumps,mags,params,x, y)
% @ bumps: [nsamples*nPCA] crosscorrelation at each point and the last 5 points
% with the Bump morphology on the time domain
% @ mags: estimated or initial Bump magnitudes [nPCA*nBumps]
% @ params: estimated or initial 2-gamma paramter

spacing=6;              % distance between peaks, originaly, spacing = 5
offset = 3;             % how soon the first peak can be or how late the last,originaly offset = 2
nstates=size(mags,2);   % number of bumps
lens=y-x+1;             % duration of each trial
ndims=size(bumps,2);    % number of PCs (PCA)
nsamples=size(bumps,1); % total data lenght N
ntrials=length(x);      % number of trials
gains=zeros(nsamples,nstates); %

for i = 1:ndims
    gains=gains+bumps(:,i)*mags(i,:)-repmat(mags(i,:).^2,nsamples,1)/2;
    % repmat->append vertically the (estimated bump-magnitudes)^2 of one PC 
    % for all samples. Then all matrix values are divided by 2.
    % bump*mags-> gives [nsamples*nBumps] It multiplies each sample of each
    % PCA from bumps' correlations with y each sample estimated of each PCA
    % of each bump. n -> Total N of samples
    % gain(n,sum(pca)) = gain(n,pca) + corrBump(n,pca) * estBumpsMorph(pca,bumps) - (estBumpMorph(pca,bumps)^2)/2
    % sum for all PCs of the 'normalized' correlation of P(having a sin) and bump morphology
end
gains=exp(gains); % gains exponential
probs=zeros(375,ntrials,nstates);
probsB=zeros(375,ntrials,nstates);
for i = 1:ntrials
    probs(offset+1:y(i)-x(i)+1-offset,i,:)=gains(x(i)+offset:y(i)-offset,:); % assign gains per trial to 'probs'
    for j = 1:nstates
        probsB(offset+1:y(i)-x(i)+1-offset,i,j)=wrev(gains(x(i)+offset:y(i)-offset,nstates+1-j)); 
        % assign the reverse of gains per trial to 'probsB'
    end
end
LP=zeros(375,nstates+1); % Gamma pdf with each stage parameters
for j = 1:nstates+1
    LP(:,j)=gammaEEG(params(j,1),params(j,2),375); % Compute Gamma pdf from 1 to 375 with parameters 'params'
end
%LP(:)=1/375;
BLP(:,:)=LP(:,nstates+1:-1:1); % reverse Gamma pdf with each stage parameter
forward = zeros(375,ntrials,nstates);
forwardB = zeros(375,ntrials,nstates);
backward = zeros(375,ntrials,nstates);
% eq1 in Appendix, first definition of likelyhood
forward(offset+1:375,:,1)=repmat(LP(1:375-offset,1),1,ntrials).*probs(offset+1:375,:,1); % gains * Gamma pdf
forwardB(offset+1:375,:,1)=repmat(BLP(1:375-offset,1),1,ntrials); % reversed Gamma pdf

% start Viturbi algorithm??

for i = 2:nstates
    next=cat(1,zeros(spacing,1),LP(1:375-spacing,i));
    nextB=cat(1,zeros(spacing,1),BLP(1:375-spacing,i));
    addB=forwardB(:,:,i-1).*probsB(:,:,i-1); % reversed gains * reversed Gamma pdf
    for j=1:ntrials
        temp=conv(forward(:,j,i-1),next);
        forward(:,j,i)=temp(1:375);
        temp=conv(addB(:,j),nextB);
        forwardB(:,j,i)=temp(1:375);
    end
    forward(:,:,i)=forward(:,:,i).*probs(:,:,i); % gains
end
forwardB=forwardB(:,:,nstates:-1:1);
for j=1:ntrials
    for i = 1:nstates
        backward(1:lens(j),j,i)=wrev(forwardB(1:lens(j),j,i));
    end
end
backward(1:offset,:,:)=0;
temp=forward.*backward; % [375,ntrials,nstates] .* [375,ntrials,nstates];
likelihood=sum(log(sum(temp(:,:,1)))); %sum(log(sum of 'temp' by columns, samples in a trial))
eventprobs=temp./repmat(sum(temp),[375,1,1]); % normalization [-1, 1] divide each trial and state by the sum of the n points in a trial




