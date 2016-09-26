function [likelihood eventprobs forward backward]=calcEEG50Condsh(bumps,conds,mags,params,x, y,map)
% @ bumps: bump values 2D [Nxch] calculated in 'calcBumps(zscores)'
% @ conds: experimental condition of each trial
% @ mags: zeros(10, 5);
% @ params: flattend version or same as 'actparams' gamma-2 scale and factor parameters if flat duration is
% constrained to mean duration equal to ACT-R model
% @ x: identifies samples where a trial starts
% @ y: identifies samples where a trial ends
% @ map: 'mapping' indicates what bumps to calculate for the ACT-R model
% per condition
% ------------------> OUT
% @ likelihood: log-likelyhood of he resulting fit
% @ eventprobs: probavility that  the five bumps are on particular points
% @ 

spacing=5;
offset = 2;
%offset is how soon the first peak can be or how late the last
%space is the distance between peaks
nconds=max(conds);
nstates=size(mags,2);
lens=y-x+1;
ndims=size(bumps,2);
nsamples=size(bumps,1);
ntrials=length(x);
gains=zeros(nsamples,nstates);

for i = 1:ndims
    gains=gains+bumps(:,i)*mags(i,:)-repmat(mags(i,:).^2,nsamples,1)/2;
end
% gains = zeros; with the given mags parameter

gains=exp(gains); % gains = 1
probs=zeros(300,ntrials,nstates);
probsB=zeros(300,ntrials,nstates);
for i = 1:ntrials
    probs(offset+1:y(i)-x(i)+1-offset,i,:)=gains(x(i)+offset:y(i)-offset,:);
    % probs = 1
    for j = 1:nstates
        probsB(offset+1:y(i)-x(i)+1-offset,i,j)=wrev(gains(x(i)+offset:y(i)-offset,nstates+1-j));
        % probsB = 1
    end
end
LP=zeros(300,nstates+1,nconds);
for j = 1:nstates+1
    for i = 1:nconds
        if j==nstates+1 || map(i,j)~=0 
            LP(:,j,i)=gammaEEG(params(j,1,i),params(j,2,i),300);
            % returns gamma probability density function normalized to one
        end
    end
end 
BLP=LP(:,nstates+1:-1:1,:);
forward = zeros(300,ntrials,nstates);
backward = zeros(300,ntrials,nstates);
for j = 1:nconds
    nj =sum(conds==j);
    lensT=lens(conds==j);
    forwardT=zeros(300,nj,nstates);
    forwardBT=zeros(300,nj,nstates);
    backwardT=zeros(300,nj,nstates);
    
    forwardBT(offset+1:300,:,1)=repmat(BLP(1:300-offset,1,j),1,nj);
    
    forwardT(offset+1:300,:,1)=repmat(LP(1:300-offset,1,j),1,nj).*probs(offset+1:300,conds==j,1);
    probsT=probs(:,conds==j,:);
    probsBT=probsB(:,conds==j,:);
    for i = 2:nstates
        if map(j,i)>0
            next=cat(1,zeros(spacing,1),LP(1:300-spacing,i,j));
            for k=1:nj
                temp=conv(forwardT(:,k,i-1),next);
                forwardT(:,k,i)=temp(1:300);
            end
            forwardT(:,:,i)=forwardT(:,:,i).*probsT(:,:,i);  
        else
            forwardT(:,:,i)=forwardT(:,:,i-1);
        end
        if map(j,nstates+2-i)>0
            nextB=cat(1,zeros(spacing,1),BLP(1:300-spacing,i,j));
            addB=forwardBT(:,:,i-1).*probsBT(:,:,i-1);
            for k=1:nj
                temp=conv(addB(:,k),nextB);
                forwardBT(:,k,i)=temp(1:300);
            end
        else
            forwardBT(:,:,i)=forwardBT(:,:,i-1);
        end
    end
    forwardBT=forwardBT(:,:,nstates:-1:1);
    for k=1:nj
        for i = 1:nstates
            backwardT(1:lensT(k),k,i)=wrev(forwardBT(1:lensT(k),k,i));
        end
    end
    backwardT(1:offset,:,:)=0;
    backward(:,conds==j,:)=backwardT;
    forward(:,conds==j,:)=forwardT;
end
temp=forward.*backward;
likelihood=sum(log(sum(temp(:,:,1))));
eventprobs=temp./repmat(sum(temp),[300,1,1]);




