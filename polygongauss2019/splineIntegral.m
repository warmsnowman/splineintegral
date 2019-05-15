curve=[0,0;
     0.25,0;
     0.5,0;
     0.75,0;
     1,0;
     1,0.25;
     1,0.5;
     1,0.75;
     1,1;
     0.75,1;
     0.5,1;
     0.25,1;
     0,1;
     0,0.75;
     0,0.5;
     0,0.25;
     0,0;];
t=[0:0.1:2*pi, 2*pi]';
curve=[cos(t),sin(t)];

% minlen=1/nSeg;
order=6;
passivefunction=@(x,y,t) ones(size(x));
% Compute product integral on the donating region from SPLINEGAUSS.
% see testDonatingRegion for its usage.

% constants for SPLINEGAUSS_2009
splType = 'not-a-knot';
cubatureDegree = order;
cubature_type=4;% guass legendre.
% -----------------------------------

spline_order_vett=[3,size(curve,1)];

[xNodes, yNodes, weights] = splinegauss(cubatureDegree, curve,...
  spline_order_vett,  splType,cubature_type);
fNodes = passivefunction(xNodes, yNodes,0);
productInt = weights'*fNodes;

if isempty(productInt)
    productInt=0;
end

 productInt-pi



% 
% 
%                 S_i1_q_ij_u=fnval(ppx,t);
%                 S_i2_q_ij_u=fnval(ppy,t);