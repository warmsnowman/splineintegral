addpath '../polygongauss2019';
addpath '../SPLINEGAUSS_2009';
poly=[1,0;
    2,0;
    3,4;
    1,5;
    -1,2;];

poly=[poly;
    poly(1,:)];
cubature_degree=2;
[nodes_x,nodes_y,weights]=gauss2019(poly,cubature_degree);
figure (1)
plot(poly(:,1),poly(:,2))
hold on
plot(nodes_x,nodes_y,'.')
hold off
pfun=@(x,y)x.*y;
pfun(nodes_x,nodes_y)'*weights
%%
%splinegauss_2009b
P = [0 0]; Q = [0 1];
rotation = 0;
cumulative = 1;
splType = 'not-a-knot';
cubature_type = 4;
spline_order_vett=[2,size(poly,1)];

[nodes_x1,nodes_y1,weights1]=splinegauss_2009b(cubature_degree,poly,...
    rotation,P,Q,spline_order_vett,cumulative,splType,cubature_type);
figure (2)
plot(poly(:,1),poly(:,2))
hold on
plot(nodes_x1,nodes_y1,'.')
hold off
pfun(nodes_x1,nodes_y1)'*weights1
