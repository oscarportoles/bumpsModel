function [likes mags params]=testSubjectsCond(zscores,conds,subjects,magstart,x,y,map)
nstates=size(magstart,2);
nsubjs=max(subjects);
nconds=max(conds);
ndims=size(zscores,2);
likes=zeros(nsubjs,1);
mags=zeros(ndims,nstates,nsubjs);
params=zeros(nstates+1,2,nconds,nsubjs);
parfor i = 1:nsubjs
    [likes(i) mags(:,:,i) params(:,:,:,i)]=testSubject(zscores,conds,subjects,nstates,magstart,1,x,y,i,map);
end
end
    function [like1 mags1 params1]=testSubject(zscores,conds,subjects,nstates,mags,thresh,x,y,i,map)
    indicesx=zeros(length(subjects),1);
    indicesy=zeros(length(subjects),1);
    indicesx(x)=1;
    indicesy(y)=1;
    params=repmat([1 150/nstates],nstates+1,1);
    a=find(subjects~=i);
    b=find(subjects(x)~=i);
    ax=find(indicesx(a)==1);
    ay=find(indicesy(a)==1);
    [~,mags1,params1]=hsmmCond(zscores(a,:),conds(b),mags,params,thresh,ax,ay,map);
    a=find(subjects==i);
    b=find(subjects(x)==i);
    ax=find(indicesx(a)==1);
    ay=find(indicesy(a)==1);
    bumps=calcBumps(zscores(a,:));
    like1=calcEEG50Condsh(bumps,conds(b),mags1,params1,ax, ay,map);    
    end
  