
function [lkhst2, magst2, paramst2 ] = DiscoverDiffStates2(SNR,mags,params,N,meanpower)
    num = length(mags);
    ncases=num^2;
    indices=reshape(cat(2,repmat(1:num,num,1),repmat([1:num]',1,num)),ncases,2);
    lkhst2 = zeros(ncases,1);
     magst2=cell(ncases,1);
    paramst2=cell(ncases,1);
    for i = 1:num
        [signals{i} x{i} y{i}] = syntheticData(ones(N,1),mags{i},params{i}, SNR ,meanpower);
    end   
    parfor i = 1:ncases
        [lkhst2(i), magst2{i}, paramst2{i}]=fitsignal(signals,x,y,mags,params,indices(i,:));
    end
    lkhst2=reshape(lkhst2,num,num)';
end

function [lkhst, magst, paramst] = fitsignal(signals,xs,ys,mags,params,indices)
    i = indices(1);
    h =indices(2);
    signal = signals{i};
    x = xs{i};
    y = ys{i};
    chunk = floor(length(x)/20);
    [~, magst, paramst]=hsmmEEG(signal(x(chunk+1):end,:),mags{h},params{h},.01,x(chunk+1:end)-y(chunk),y(chunk+1:end)-y(chunk));
    lkhst=hsmmEEG(signal(1:y(chunk),:),magst,paramst,0,x(1:chunk),y(1:chunk));
end
