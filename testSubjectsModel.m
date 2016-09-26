function [likes mags]=testSubjectsModel(zscores,conds,subjects,magstart,modelparams,x,y,map)
nstates=size(magstart,2);
nsubjs=max(subjects);
ndims=size(zscores,2);
likes=zeros(nsubjs,1);
mags=zeros(ndims,nstates,nsubjs);
parfor i = 1:nsubjs
    [likes(i) mags(:,:,i)]=testSubjectModel(zscores,conds,subjects,magstart,modelparams,x,y,i,map);
end
end
    function [like1 mags1]=testSubjectModel(zscores,conds,subjects,mags,modelparams,x,y,i,map)
    indicesx=zeros(length(subjects),1);
    indicesy=zeros(length(subjects),1);
    indicesx(x)=1;
    indicesy(y)=1;
    a=find(subjects~=i);
    ax=find(indicesx(a)==1);
    ay=find(indicesy(a)==1);
    b=find(subjects(x)~=i);
    [~,mags1]=hsmmEEGModel(zscores(a,:),conds(b),mags,modelparams,1,ax,ay,map);
    a=find(subjects==i);
    ax=find(indicesx(a)==1);
    ay=find(indicesy(a)==1);
    bumps=calcBumps(zscores(a,:));
    b=find(subjects(x)==i);
    like1=calcEEG50Condsh(bumps,conds(b),mags1,modelparams,ax,ay,map);
    end
    
    