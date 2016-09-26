function [lkh1 mags1 eventprobs1]=hsmmEEGModel(zscores,conds,mags,modelparams,thresh,x,y,map)
% @ zscores: normalized PCA of all trials concatenated
% @ conds: experimental condition of each trial
% @ mags: zeros(10, 5); initialization
% @ modelparams: 'actparams' gamma-2 scale parameter if flat duration is
% constrained mean duration equal to ACT-R model
% @ thershold: 1, 
% @ x: identifies samples where a trial starts
% @ y: identifies samples where a trial ends
% @ map: 'mapping' indicates what bumps to calculate for the ACT-R model
% per condition
% -----------> OUT
% @ lkh1: log-likelyhood of the resulting fit
% @ mags5: it gives the means of the bumps [5x (ch or Cs)]
% @ eventprobs1: It gives for each trial(dim-2) the probabilities that the
% 5 bumps (dim-3) are centered on particular locations (dim-1)

bumps=calcBumps(zscores); % compute bumps template by time seris (components) There is kind of regularization

nconds=max(conds);
% resiz modelparam to a flat 3D matrix
if length(size(modelparams))==2;
    modelparams=repmat(modelparams,[1,1,nconds]);
end
lkh1=-Inf;
nstates=size(mags,2);
ncases=length(x);
ndims=size(bumps,2);
lens=y-x+1;
means=zeros(300,ncases,ndims);
lengths=zeros(nconds,1);

for i = 1:nconds
    lengths(i)=mean(lens(conds==i)); % mean lenght of one condition
end

for i = 1:ncases
    means(1:lens(i),i,:)=bumps(x(i):y(i),:); % segment bumps
end

[lkh eventprobs]=calcEEG50Condsh(bumps,conds,mags,modelparams,x,y,map); % initializa first iteration

if thresh==0
    lkh1 = lkh;
    eventprobs1=eventprobs;
    mags1=mags;
else
    while (lkh-lkh1)>thresh
        lkh1=lkh;
        mags1 = mags;
        eventprobs1=eventprobs;
        for i = 1:nstates
            a=find(map(conds,i)>0);
            for j = 1:ndims % Cs
                mags(j,i)=mean(sum(eventprobs(:,a,i).*means(:,a,j)));
            end
        end
        [lkh eventprobs]=calcEEG50Condsh(bumps,conds,mags,modelparams,x,y,map);
    end
end
    
    
