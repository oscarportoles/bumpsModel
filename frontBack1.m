function [fronts backs]= frontBack1(data,x,y,conds)
nchan=size(data,2);
nconds=max(conds);
fronts=zeros(100,nchan,nconds);
backs=zeros(100,nchan,nconds);
lens=y-x+1;
for i = 1:nconds
    xa=x(conds==i);
    ya=y(conds==i);
    lensa=lens(conds==i);
    ncases=length(xa);
    temp=zeros(100,nchan,3);
    for j = 1:ncases
        for k = 1:100
            if lensa(j)>=k
                temp(k,:,1)=temp(k,:,1)+data(xa(j)+k-1,:);
                temp(k,:,2)=temp(k,:,2)+data(ya(j)-k+1,:);
                temp(k,:,3)=temp(k,:,3)+1;
            end
        end
    end
    fronts(:,:,i)=temp(:,:,1)./temp(:,:,3);
    backs(:,:,i)=temp(:,:,2)./temp(:,:,3);
    for j = 1:nchan
        backs(:,j,i)=wrev(backs(:,j,i));
    end
end
    
    