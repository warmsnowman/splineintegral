function [xy,ii]=SlefxSelf(a)
%a is 2*n matrix
%suppose the n th edge never intersect with the n+1 th edge
xx=[];
yy=[];
ii=[];

i=1;
for j=i+2:size(a,2)-2
    [tx,ty]=polyxpoly(a(1,i:i+1), a(2,i:i+1),a(1,j:j+1), a(2,j:j+1));
    if ~isempty(tx)
        ti=[i;j];
        xx=[xx,tx];
        yy=[yy,ty];
        ii=[ii,ti];
    end
end

for i=2:size(a,2)-3
    for j=i+2:size(a,2)-1
        [tx,ty]=polyxpoly(a(1,i:i+1), a(2,i:i+1),a(1,j:j+1), a(2,j:j+1));
        if ~isempty(tx)
            ti=[i;j];
            xx=[xx,tx];
            yy=[yy,ty];
            ii=[ii,ti];
        end
    end
end
xy=[xx;yy];
end