function xi = findxi(control_points, control_points_down,point)
%since it is down chain, we have p_i>p_{i+1}
xiset=[];
for i=1:length(control_points_down)
if  control_points(control_points_down(i),2)>y ...
        && y>control_points(control_points_down(i)+1,2)   
    
               x1=control_points(control_points_down(i),1);
            x2=control_points(control_points_down(i)+1,1);
            
            y1=control_points(control_points_down(i),2);
            y2=control_points(control_points_down(i)+1,2); 
    
    
    temp=x1+(y1-point(2))/(y1-y2)*(x2-x1);
    
    xiset=[xiset; temp];
    
    
    
end

%choice the most close one.
[~, index]=min(abs(xiset-point(1)));
xi=xxiset(index);

end

