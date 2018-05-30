curve=[0,0.25,0.5,0.75,1;0,0,0,0,0];
curve=fliplr(curve)';
% curve=[0,1;0,0];
nSeg=60;

minlen=1/nSeg;
order=4;
passivefunction=@(x,y,t) ones(size(x));
% Compute product integral on the donating region from SPLINEGAUSS.
% see testDonatingRegion for its usage.

% constants for SPLINEGAUSS_2009
splType = 'not-a-knot';
cubatureDegree = order+18;
cubature_type=4;% guass legendre.
% -----------------------------------

spline_order_vett=[3,5];

[xNodes, yNodes, weights] = splinegauss(cubatureDegree, curve,...
  spline_order_vett,  splType,cubature_type);
fNodes = passivefunction(xNodes, yNodes,0);
productInt = weights'*fNodes

if isempty(productInt)
    productInt=0
end

 