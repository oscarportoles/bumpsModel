function [likes mags params]=testSubjects(zscores,subjects,magstart,x,y)
% @ zscores: z-scored time series (10 PCA)
% @ subjects: [Nsamples] vector with subject ID
% @ magstart: [nPC, nBumps] Initial bump magnitudes, best from all subject study?
% @ x: begining of each trial
% @ y: end of each trial

nstates=size(magstart,2);           % number of bups == stages
nsubjs=max(subjects);               % number fo subjects      
ndims=size(zscores,2);              % number of PCs
likes=zeros(nsubjs,1);              
mags=zeros(ndims,nstates,nsubjs);
params=zeros(nstates+1,2,nsubjs);

parpool(nsubjs)
parfor i = 1:nsubjs
    [likes(i) mags(:,:,i) params(:,:,i)]=testSubject(zscores,subjects,nstates,magstart,1,x,y,i);
end
end

    function [like1 mags1 params1]=testSubject(zscores,subjects,nstates,mags,thresh,x,y,i)
    indicesx=zeros(length(subjects),1);
    indicesy=zeros(length(subjects),1);
    indicesx(x)=1;  % puts a 1 at the begining of a trial
    indicesy(y)=1;  % puts a 1 at the end of a trial
    %%% gamma parameter = 1  ??????!!!!!
    params=repmat([2 188/nstates],nstates+1,1); % originaly 188 -> 150, half a trial length
    a=find(subjects~=i);        % Idx of all but one subject (out)
    ax=find(indicesx(a)==1);    % Idx begining of trials of subjects in
    ay=find(indicesy(a)==1);    % Idx end of trials of subjetcs in
    
    % Do HSMM of subjects In
    [~,mags1,params1]=hsmmEEG(zscores(a,:),mags,params,thresh,ax,ay);
    
    % takes subject Out
    a=find(subjects==i);
    ax=find(indicesx(a)==1);
    ay=find(indicesy(a)==1);
    % Do likelyhood of fit of subject In to the subject Out
    bumps=calcBumps(zscores(a,:));
    like1=calcEEG50h(bumps,mags1,params1,ax, ay); % fit of subject out
    end
    
    