function d = gammaEEG(a,b,maxlength)
% returns gamma probability density function normalized to one
d=zeros(maxlength,1);
for i = 1:maxlength
    d(i)=gampdf(i-.5,a,b);
end
d = d/sum(d);

