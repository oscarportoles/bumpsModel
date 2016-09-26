function bumps = calcBumps(data)
% it puts on each sample the correlation of that sample and the previous
% five samples with a Bump morphology on time domain.
% @ data: 'zscores' [nxch] a normalizad PCA from all trials concatenated
% This function has been modified to accpet 125 Hz sampling frequency

% Bump normalization sum(P) = 1.294.  
%template=[0.3090    0.8090    1.0000    0.8090    0.3090]';
template = [0.309016994374947;0.728968627421412;0.968583161128631;0.968583161128631;0.728968627421411;0.309016994374948];
template=template/sum(template.^2);

width = size(data,2); % number of EEG channels or PCA components
nsamples = size(data,1); % total samples
bumps=zeros(size(data));

for j = 1:width % PCA components
    temp=zeros(nsamples,6);
    temp(:,1)=data(:,j);
    for i = 2:6
        temp(:,i)=cat(1,temp(2:end,i-1),0); % concatante along 1D
        % puts the a component in a [nsamples X length(bump)] matrix shifted.
        % each column is a copy of the first one but shifted one sample
        % upwards
    end
    bumps(:,j)=temp*template;
end

%bumps(3:end,:)=bumps(1:end-2,:); old
%bumps([1,2,nsamples-1,nsamples],:)=0; old
bumps(4:end,:)=bumps(1:end-3,:); % the edges are centered (changed)
bumps([1,2,3,nsamples-2,nsamples-1,nsamples],:)=0; % the leakage from the edges gets zeros (changed)