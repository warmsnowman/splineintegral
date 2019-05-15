function [pts] = DeleteSamePts(pts,minDist)
% delete the adjacent same points in the polygonal representation of curve.

%input preconditioning
if exist('minDist','var')==0
    minDist=eps;%let minDist be the machine number if input is missing.
end
%%
i=1;
while i<size(pts,1) %execute the loop if i less than the number of points.
    if norm(pts(i,:)-pts(i+1,:))<minDist %the adjacent points are same.
        %ie. the distance between them is less the minDist.
        %if sum((pts(:,i)-pts(:,i+1)).^2)<eps
            pts(i+1,:)=[];%delete the latter one.
      
    else
        i=i+1;
    end
end

end
