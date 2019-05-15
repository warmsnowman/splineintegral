function [nodes_x,nodes_y,weights]=gauss2019(poly,cubature_degree)

%%
%preconddition
%delete same point
poly = DeleteSamePts(poly);
if isempty(SlefxSelf(poly'))
else
    warning('polygon is self interset, the cubature rule may not right.') 
end
%%
%We assume the orientation is anticlockwise, if not, convert it to anticlockwise

[maxy,index]=max(poly(:,1));
if poly(index,2)<poly(index+1,2)%y coordinate
  orientation=-1;%anticlockwise
else
    orientation=1;%clockwise
end

% if orientation is anticlockwise, we don't need to do anything 
%otherwise, filp the vector.
if orientation==1
   poly= udflip(poly);
end

%%
%find the ascending chain, descending chain, horizontal chain
polyy=poly(:,2);
polyx=poly(:,1);
signdify=sign(diff(polyy));
udchain=[];
for i=1:length(signdify)-1
    if signdify(i)~=signdify(i+1)
            udchain=[udchain; i+1, signdify(i);];  
    end       
end
udchain=[udchain;length(signdify)+1, signdify(end);];


udchain=[zeros(size(udchain,1),1),udchain];
udchain(1,1)=1;
for i=2:size(udchain,1)
    udchain(i,1)=udchain(i-1,2);
end
% udchain(i,3) indicates the type of chain.
% 1 ascending chain; -1 descending chain; 0 horizontal chain
% The chain consists of the points from udchain(i,1) to udchain(i,2)
% is udchain(i,3) chain for i=1:length(udchain)
achain=udchain(udchain(:,3)==1,:);
dchain=udchain(udchain(:,3)==-1,:);


%%
%assume there is only one ascending chain and descending chain. 
%use y coordinate as the dependent variable. 
tempy=unique([poly(achain(1,1):achain(1,2),2);poly(dchain(1,1):dchain(1,2),2)]);
%tempy is a ascending sequence. 
achainpp=[interp1(poly(achain(1,1):achain(1,2),2),poly(achain(1,1):achain(1,2),1),tempy),tempy];
%%
cubature_type=4; %GAUSS-LEGENDRE (TENSORIAL)

% cubature_manager(2,1,4)
%The second argument is the degree of spline used to approxiamting the boundary.
%If the boundary is polygon, the degree is 1. 
nodes_x=[];
nodes_y=[];
weights=[];
spline_order_vett=[2,size(achainpp,1)];
L=size(spline_order_vett,1);
%It will be useful when we study the polynomial approximated boundary.
control_points=achainpp;
for block_index=1:L
    spline_block_order=spline_order_vett(block_index,1);
    if spline_block_order==1
        % 2 is straight line.
        % 1 is nothing. 
        continue
    end
    spline_block_degree=spline_block_order-1;
    
    % Initial and final indices of "control points" in the block.
    if (block_index ==1)
        initial_point_index=1;
    else
        initial_point_index=spline_order_vett(block_index-1,2);
    end
    final_point_index=spline_order_vett(block_index,2);
    
    % Spline order in the block.
    
    % Control points (x_loc,y_loc) in the block.
    pts_loc=control_points(initial_point_index:final_point_index,:);
    
    % Parametrical description of the block.
    
    s_loc=[0;cumsum(edgeLength(pts_loc))];
    
    
    
    
    
    % Computing the spline parametrical description of the block.
    % "ppx", "ppy" describe S_i1, S_i2 in the paper, while "ppy1"
    % describe S'_i2.
    
    switch spline_block_order
        case 4            
            % CUBIC SPLINES BY CSAPE. AS DEFAULT WE USE PERIODIC CUBIC SPLINES.
            % Derivatives parameters are computed as well.
            ppx=csape(s_loc,pts_loc(:,1),SPLtypestring);
            ppy=csape(s_loc,pts_loc(:,2),SPLtypestring);
            [breaks_y,coeffs_y]=unmkpp(ppy);
            N_y=size(coeffs_y,1);
            dcoeffs_y=[zeros(N_y,1) 3*coeffs_y(:,1) 2*coeffs_y(:,2) ...
                coeffs_y(:,3)];
            ppy1=mkpp(breaks_y,dcoeffs_y);
            
        otherwise
            
            ppx=spapi(spline_block_order,s_loc,pts_loc(:,1));
            ppy=spapi(spline_block_order,s_loc,pts_loc(:,2));
            ppy1=fnder(ppy,1);
            
    end
    
    
    
    % Every block is subdivided in "number_of_subblocks" curves determined
    % by successive control points.
    number_of_subblocks=final_point_index-initial_point_index;
    
    % Cubature rule on the square [-1,1] x [-1,1].
    % Padua Points: 0, Gauss-Legendre: 4.
    [x_pts, y_pts, wpd]=cubature_manager(cubature_degree,spline_block_degree,...
        cubature_type);
 
    % Computing quadrature points from a general sub-block. The cases in
    % which the order is 2 is a little different from other spline orders.
    % Consequently, we distinguish between them.
    for index_control_point=1:number_of_subblocks
        
        if (spline_block_order == 2)
            
            x1=pts_loc(index_control_point,1);
            x2=pts_loc(index_control_point+1,1);
            
            y1=pts_loc(index_control_point,2);
            y2=pts_loc(index_control_point+1,2);
            
               %yi is monotonic sequence.                     
               % Computing nodes.
                    
                    xeta=(x2-x1)/2.*y_pts+(x2+x1)/2;
                    y=(y2-y1)/2.*y_pts+(y2+y1)/2;
                    g=lowerlimit(y,poly,dchain );
                    x=(xeta-g)/2.*x_pts+(xeta+g)/2;
                    
                    nodes_x=[nodes_x; x];
                    nodes_y=[nodes_y; y];
                    
                    diff_y=(y2-y1)/2;
                    diff_lambda=(xeta-g)/2;
                    local_weights=diff_y.*...
                        diff_lambda.*wpd;
                    
                    weights=[weights; local_weights];
                    

  
            
        else
            % spline_block_order ~= 2.
            
       error('We have not finished this subroutine for the higher degree polynomial.')
        end
    end
    
end

%delete zero weights
zeroweights=(weights==0);
nodes_x(zeroweights)=[];
nodes_y(zeroweights)=[];
weights(zeroweights)=[];


end