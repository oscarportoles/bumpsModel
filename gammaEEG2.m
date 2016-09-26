function d = gammaEEG2(a,b,maxlength)
d=zeros(maxlength,1);
for i = 1:maxlength
    d(i)=gampdf(i,a,b);
end
d = d/sum(d);

