function [condsB, x5, y5] = getVals(conds,x,y)
condsB=zeros(y(end),1);
ncases=length(x);
for i = 1:ncases
    condsB(x(i):y(i))=conds(i);
end
x5=cell(5,1);
y5=cell(5,1);
for i = 1:5
    lens=y(conds==i)-x(conds==i)+1;
    x5{i}=zeros(size(lens));
    y5{i}=zeros(size(lens));
    x5{i}(1)=1;
    y5{i}(1)=lens(1);
    for j =2:length(lens)
        x5{i}(j)=y5{i}(j-1)+1;
        y5{i}(j)=y5{i}(j-1)+lens(j);
    end
end