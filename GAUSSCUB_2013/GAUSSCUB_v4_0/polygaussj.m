
function [cubature_val,nodes_x,nodes_y,weights]=polygaussj(intfcn,N,polygon_sides,rotation,P,Q,quadf,cubature_type,varargin)

%----------------------------------------------------------------------------------------------------
% REFERENCE PAPER [1] A. SOMMARIVA and M. VIANELLO "Gauss-like and triangulation-free cubature over 
%                                                   polygons".
%
% INPUT:
%
% [intfcn] : IT IS AN "INLINE" FUNCTION OR A FUNCTION DEFINED BY THE USER. FOR EXAMPLE, IF THE m-FILE 
%            OF SUCH A FUNCTION IS NAMED "fct2d.m" ONE HAS TO PUT "intfct" EQUAL TO "@fct2d".
%
% [N]      : DEGREE OF THE 1 DIMENSIONAL QUADRATURE RULE.
%
% [polygon_sides]: IF THE POLYGON HAS "L" SIDES, "boundary.pts" IS A VARIABLE CONTAINING ITS VERTICES,
%            ORDERED COUNTERCLOCKWISE. AS LAST ROW MUST HAVE THE COMPONENTS OF THE FIRST VERTEX. 
%            IN OTHER WORDS, THE FIRST ROW AND LAST ROW ARE EQUAL. 
%            "polygon_sides" IS A "L+1 x 2" MATRIX.
%
% [rotation]: SUPPOSE x_min, x_max, y_min, y_max ARE THE MINIMUM AND MAXIMUM VALUES IN x
%            AND y REACHED BY THE POLYGON, I.E. "R=[x_min,x_max] x [y_min,y_max]" IS THE  
%            SMALLEST RECTANGLE WITH SIDES PARALLEL TO THE AXIS x AND y, CONTAINING THE POLYGON.
%
%            [rotation=0]: GENERAL CASE, BUT OBSERVE THAT THE INTEGRAND MUST BE DEFINED
%                          IN THE RECTANGLE "R" DESCRIBED ABOVE.
%            [rotation=1]: GOOD FOR CONVEX POLYGONS. FOR CONVEX POLYGONS IT SUFFICES 
%                          THAT THE INTEGRAND BE DEFINED IN THE POLYGON.
%            [rotation=2]: THE USER CHOOSES A SPECIAL REFERENCE SEGMENT "PQ". CHECK [1]  
%                          FOR FURTHER DETAILS. SEE THE VARIABLES "P", "Q" BELOW.
%
% [P, Q]   : IF [rotation=2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT "PQ" 
%            HAVING "P", "Q" AS EXTREMA. 
%            "P" AND "Q" ARE "1 x 2" ROW VECTORS.
%
% [quadf]  : IF ONE USES "inline" FUNCTIONS PUT "[]".
%
% [cubature_type]: IT CHOOSES THE 1D QUADRATURE RULE .
%
%            [cubature_type=1]: FEJER 1.
%            [cubature_type=2]: FEJER 2.
%            [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%            [cubature_type=4]: GAUSS-LEGENDRE.
%
% OUTPUT:
%
% [cubature_val]: IS THE VALUE OF THE CUBATURE OF THE FUNCTION SPECIFIED IN "intfcn" ON THE BOUNDARY 
%            DESCRIBED IN "polygon_sides" BY A QUADRATURE RULE OF DEGREE "N".
%
% [nodes_x, nodes_y]: THE CUBATURE RULE PRODUCES THE NODES "(nodes_x,nodes_y)".
%            IF THE CUBATURE RULE GENERATES "K" NODES, "nodes_x" AND "nodes_y" ARE "K x 1" COLUMN
%            VECTORS.
%
% [weights]: THE  CUBATURE RULE PRODUCES THE WEIGHTS "weights" (RELATIVE TO THE NODES 
%            "(nodes_x,nodes_y)".
% 
%
%
% EXAMPLE 1: USING "inline" FUNCTIONS.
% [cubature_result,x,y,weights]=polygaussj(inline('x+3*y'),N,polygon_sides,rotation,P,Q);
%
% EXAMPLE 2: USING AN "m-FUNCTION" "fct2D.m" THAT REQUIRES AS INPUT THE VARIABLE "function_type".
% [cubature_result,x,y,weights]=polygaussj(@fct2D,N,polygon_sides,rotation,P,Q,[],cubature_type,function_type);
%
%----------------------------------------------------------------------------------------------------
%% Copyright (C) 2007 Marco Vianello and Alvise Sommariva
%%
%% This program is free software; you can redistribute it and/or modify
%% it under the terms of the GNU General Public License as published by
%% the Free Software Foundation; either version 2 of the License, or
%% (at your option) any later version.
%%
%% This program is distributed in the hope that it will be useful,
%% but WITHOUT ANY WARRANTY; without even the implied warranty of
%% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%% GNU General Public License for more details.
%%
%% You should have received a copy of the GNU General Public License
%% along with this program; if not, write to the Free Software
%% Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA  02111-1307  USA

%% Authors:  
%% Marco Vianello    <marcov@euler.math.unipd.it>
%% Alvise Sommariva  <alvise@euler.math.unipd.it>
%% Date: March 7, 2007
%----------------------------------------------------------------------------------------------------


%----------------------------------------------------------------------
% BOUNDARY PTS.
%----------------------------------------------------------------------
x_bd=polygon_sides(:,1);
y_bd=polygon_sides(:,2);

%----------------------------------------------------------------------
% "MINIMUM" RECTANGLE CONTAINING POLYGON.
%----------------------------------------------------------------------
x_min=min(x_bd); x_max=max(x_bd);
y_min=min(y_bd); y_max=max(y_bd);

%----------------------------------------------------------------------
% SOME AUTOMATIC SETTINGS.
%----------------------------------------------------------------------
if nargin == 3 rotation=0; P=[x_min y_max]; Q=[x_min y_min]; end
if nargin == 6 quadf = @fct2D; end
if nargin == 8 cubature_type=4; end
if isempty(quadf), quadf = @fct2D; end

%--------------------------------------------------------------------------
% POLYGON ROTATION (IF NECESSARY).
%--------------------------------------------------------------------------
fprintf('\n \t ----------------------------------------------------------------------------');
switch rotation
case 0
    fprintf('\n \t [ROTATION]: NO.');
    rot_matrix=eye(2);
    axis_abscissa=[x_min y_max]-[x_min y_min];
case 1
    fprintf('\n \t [ROTATION]: AUTOMATIC');
    [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=auto_rotation(polygon_sides,[],[]);
    fprintf(' [ANGLE CLOCKWISE (RESPECT Y, IN DEGREES)]: %5.5f',rot_angle*180/pi);
case 2
    fprintf('\n \t [ROTATION]: PREFERRED DIRECTION');
    nrm_vect=norm(Q-P);
    if (nrm_vect > 0)
        direction_axis=(Q-P)/nrm_vect;
        [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=auto_rotation(polygon_sides,P,Q);
    else
        fprintf('\n \t [WARNING]: THE DIRECTION VECTOR IS NULL, I USE AUTOMATIC ROTATION.');
        [polygon_sides,rot_matrix,rot_angle,axis_abscissa,P,Q]=auto_rotation(polygon_sides,P,Q);
    end
    fprintf(' [ANGLE CLOCKWISE (RESPECT Y)]: %5.5f',rot_angle*180/pi);
end

%--------------------------------------------------------------------------
% COMPUTE NODES AND WEIGHTS OF 1D QUADRATURE RULE.
%--------------------------------------------------------------------------

% DEGREE "N".
[s_N,w_N]=cubature_rules_1D((N-1),cubature_type);
N_length=length(s_N);
M=N+1;

%----------------------------------------------------------------------
% L: NUMBER OF SIDES OF THE POLYGON.
% N: DEGREE 1D QUADRATURE PRIMITIVE.
%----------------------------------------------------------------------
L=length(polygon_sides(:,1))-1;

%----------------------------------------------------------------------
% COMPUTE 2D NODES (nodes_x,nodes_y) AND WEIGHTS "weights".
%----------------------------------------------------------------------

nodes_x=[];
nodes_y=[];
weights=[];

for index_side=1:L
    x1=polygon_sides(index_side,1); x2=polygon_sides(index_side+1,1);
    y1=polygon_sides(index_side,2); y2=polygon_sides(index_side+1,2);
    
    a=axis_abscissa(1);
    
    if ~(x1 == a & x2 == a)
        if (y2-y1) ~=0  
            
            if (x2-x1) ~=0 
                % THE SIDE OF THE POLYGON IN x VARIABLE IS
                %         A_coeff*t+B_coeff
                % FOR t IN [-1,1].
                A_coeff=(x2-x1)/2; B_coeff=(x1+x2)/2; 
                
                % THE SIDE OF THE POLYGON IN y VARIABLE IS
                %    C_coeff*t+D_coeff
                % FOR t IN [-1,1].
                C_coeff=(y2-y1)/2; D_coeff=(y1+y2)/2; 
                
                
                c=(a-B_coeff)/A_coeff;                % MODIFICATION COEFF. OF
                % GAUSS-JACOBI.
                
                [x_modj,w_modj]=mod_jacobi(M,0,0,c);  % MODIFIED GAUSS-JACOBI: 
                % a=0, b=0 (IN [-1,1]).
                s_M_loc=x_modj;                       % G.-J. NODES.
                w_M_loc=w_modj;                       % G.-J. WEIGHTS.
                M_length=length(s_M_loc);
                
                %--------------------------
                % COMPUTING CUBATURE NODES.
                %--------------------------        
                x_quad_side=B_coeff+A_coeff*s_M_loc;  % SIZE: (M_length,1)
                y_quad_side=D_coeff+C_coeff*s_M_loc;  % SIZE: (M_length,1)    
                
                scaling_fact_plus=(x_quad_side+a)/2;          % SIZE: (M_length,1)
                scaling_fact_minus=(x_quad_side-a)/2;         % SIZE: (M_length,1)
                
                term_1=repmat(scaling_fact_plus,1,N_length);   % SIZE: (M_length,N_length)
                term_2=repmat(scaling_fact_minus,1,N_length);  % SIZE: (M_length,N_length)
                
                rep_s_N=repmat(s_N',M_length,1);
                x=term_1+term_2.*rep_s_N;
                y=repmat(y_quad_side,1,N_length);
                
                %----------------------------
                % COMPUTING CUBATURE WEIGHTS.
                %----------------------------                  
                leading_term=A_coeff*C_coeff/2;       % I=leading_term*\int_-1^1 \int_-1^1 ...
                half_length_y=C_coeff;
                local_weights=leading_term*w_M_loc;
                
            else
                s_M_loc=s_N;
                w_M_loc=w_N;
                M_length=length(s_M_loc);
                
                half_pt_x=(x1+x2)/2; half_pt_y=(y1+y2)/2;
                half_length_x=(x2-x1)/2; half_length_y=(y2-y1)/2;
                
                % QUADRATURE POINTS ON THE SIDE.
                x_quad_side=half_pt_x+half_length_x*s_M_loc;  % SIZE: (M_length,1)
                y_quad_side=half_pt_y+half_length_y*s_M_loc;  % SIZE: (M_length,1)
                
                scaling_fact_plus=(x_quad_side+a)/2;          % SIZE: (M_length,1)
                scaling_fact_minus=(x_quad_side-a)/2;         % SIZE: (M_length,1)
                
                local_weights=(half_length_y*scaling_fact_minus).*w_M_loc; % SIZE: (M_length,1)
                
                term_1=repmat(scaling_fact_plus,1,N_length);   % SIZE: (M_length,N_length)
                term_2=repmat(scaling_fact_minus,1,N_length);  % SIZE: (M_length,N_length)
                
                rep_s_N=repmat(s_N',M_length,1);
                
                % x, y ARE STORED IN MATRICES. A COUPLE WITH THE SAME INDEX IS A POINT,
                % i.e. "P_i=(x(k),y(k))" FOR SOME "k".
                x=term_1+term_2.*rep_s_N;
                y=repmat(y_quad_side,1,N_length);
            end          
            
            number_rows=size(x,1);
            number_cols=size(x,2);
            
            x=x(:); x=x';
            y=y(:); y=y';
            
            rot_gauss_pts=rot_matrix'*[x;y]; % THE INVERSE OF A ROTATION MATRIX IS ITS TRANSPOSE.
            
            x_rot=rot_gauss_pts(1,:);        % GAUSS POINTS IN THE ORIGINAL SYSTEM.
            y_rot=rot_gauss_pts(2,:);
            
            x_rot=reshape(x_rot',number_rows,number_cols);
            y_rot=reshape(y_rot',number_rows,number_cols);
            
            
            nodes_x=[nodes_x; x_rot];
            nodes_y=[nodes_y; y_rot];
            weights=[weights; local_weights];
            
        end
    end
end

f = fcnchk(intfcn);

f_xy = feval(f,nodes_x,nodes_y, varargin{:});

cubature_val=(weights'*f_xy)*w_N; % COMPUTING CUBATURE.







%----------------------------------------------------------------------
% FUNCTIONS USED IN THE ALGORITHM.
%----------------------------------------------------------------------


%----------------------------------------------------------------------
% 1. "auto_rotation"
%----------------------------------------------------------------------
function [polygon_bd_rot,rot_matrix,rot_angle,axis_abscissa,vertex_1,vertex_2]=auto_rotation(polygon_bd,vertex_1,vertex_2)


% AUTOMATIC ROTATION OF A CONVEX POLYGON SO THAT "GAUSSIAN POINTS", AS IN THE PAPER ...
% ARE ALL CONTAINED IN THE CONVEX POLYGON. SEE THE PAPER FOR DETAILS.


% FIND DIRECTION AND ROTATION ANGLE.
if length(vertex_1) == 0
    % COMPUTING ALL THE DISTANCES BETWEEN POINTS.A LITTLE TIME CONSUMING AS PROCEDURE.
    distances = points2distances(polygon_bd);  
    [max_distances,max_col_comp]=max(distances,[],2);
    [max_distance,max_row_comp]=max(max_distances,[],1);
    vertex_1=polygon_bd(max_col_comp(max_row_comp),:);
    vertex_2=polygon_bd(max_row_comp,:);
    direction_axis=(vertex_2-vertex_1)/max_distance;
else
    direction_axis=(vertex_2-vertex_1)/norm(vertex_2-vertex_1);
end

rot_angle_x=acos(direction_axis(1));
rot_angle_y=acos(direction_axis(2));

if rot_angle_y <= pi/2
    if rot_angle_x <= pi/2
        rot_angle=-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
else
    if rot_angle_x <= pi/2
        rot_angle=pi-rot_angle_y;
    else
        rot_angle=rot_angle_y;
    end
end


% CLOCKWISE ROTATION.
rot_matrix=[cos(rot_angle) sin(rot_angle); 
    -sin(rot_angle) cos(rot_angle)];

number_sides=size(polygon_bd,1)-1;

polygon_bd_rot=(rot_matrix*polygon_bd')';

axis_abscissa=rot_matrix*vertex_1';



%----------------------------------------------------------------------
% 3. "cubature_rules_1D"
%----------------------------------------------------------------------

function [nodes,weights]=cubature_rules_1D(n,cubature_type)

% SEE WALDVOGEL PAPER. ADDED NODES

% Weights of the Fejer2, Clenshaw-Curtis and Fejer1 quadrature by DFTs
% n>1. Nodes: x_k = cos(k*pi/n)

N=[1:2:n-1]'; l=length(N); m=n-l; K=[0:m-1]';

switch cubature_type
    
case 1 % FEJER 1.
    v0=[2*exp(i*pi*K/n)./(1-4*K.^2); zeros(l+1,1)];
    v1=v0(1:end-1)+conj(v0(end:-1:2)); 
    weights=ifft(v1);
    k=(1/2):(n-(1/2)); nodes=(cos(k*pi/n))'; 
    
case 2 % FEJER 2.
    v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
    v2=-v0(1:end-1)-v0(end:-1:2); 
    wf2=ifft(v2); weights=[wf2;0];
    k=0:n; nodes=(cos(k*pi/n))'; 
    
case 3 % CLENSHAW CURTIS.
    g0=-ones(n,1); g0(1+l)=g0(1+l)+n; g0(1+m)=g0(1+m)+n;
    g=g0/(n^2-1+mod(n,2));
    v0=[2./N./(N-2); 1/N(end); zeros(m,1)];
    v2=-v0(1:end-1)-v0(end:-1:2);
    wcc=ifft(v2+g); weights=[wcc;wcc(1,1)];
    k=0:n; nodes=(cos(k*pi/n))'; 
    
case 4 % GAUSS LEGENDRE
    beta=0.5./sqrt(1-(2*(1:n)).^(-2));
    T=diag(beta,1)+diag(beta,-1);
    [V,D]=eig(T);
    x=diag(D); [x,i]=sort(x); x=x';
    w=2*V(1,i).^2; 
    nodes=x';
    weights=w';
    
end








%----------------------------------------------------------------------
% 3. "points2distances"
%----------------------------------------------------------------------

function distances = points2distances(points)  

% Create euclidean distance matrix from point matrix.

% Get dimensions.
[numpoints,dim]=size(points);                    

% All inner products between points.
distances=points*points';                     

% Vector of squares of norms of points.
lsq=diag(distances);                            

% Distance matrix.
distances=sqrt(repmat(lsq,1,numpoints)+repmat(lsq,1,numpoints)'-2*distances); 





%----------------------------------------------------------------------
% 4. "mod_jacobi"
%----------------------------------------------------------------------


function [x_modj,w_modj]=mod_jacobi(N,a,b,c)

% IT COMPUTES THE NODES AND WEIGHTS OF A MODIFIED GAUSS-JACOBI RULE
% w(x)=(x-c) * (1-x)^a * (1+x)^b, x \in [-1,1].

% INPUT:
%
% N: DEGREE.
% a: JACOBI EXPONENT.
% b: JACOBI EXPONENT.
%
% OUTPUT:
% x_modj: MODIFIED JACOBI NODES.
% w_modj: MODIFIED JACOBI WEIGHTS.

ab_legendre=r_jacobi(N,a,b);         % GAUSS-LEGENDRE RECURSION TERMS.
ab_modjacobi=chri1(N-1,ab_legendre,c); % GAUSS-JACOBI MOD. RECURSION TERMS.
xw_modjacobi=gauss(N-1,ab_modjacobi);  % GAUSS-JACOBI MOD. WEIGHTS AND NODES.

x_modj=xw_modjacobi(:,1);            % GAUSS-JACOBI MOD. NODES.
w_modj=xw_modjacobi(:,2); 




%----------------------------------------------------------------------
% 5. "chri1"
%----------------------------------------------------------------------

% CHRI1 Modification by a linear factor.
%
%    Given a weight function w(t) through the first n+1 recurrence
%    coefficients ab0 of its orthogonal polynomials, ab=CHRI1(n,
%    ab0,z) generates the first n recurrence coefficients for
%    the orthogonal polynomials relative to the modified weight
%    function (t-z)w(t). The alpha- and beta-coefficients of the
%    original weight function are to be provided in the first and
%    second column of the (n+1)x2 array ab0; those of the modified
%    weight function are returned in the first and second column 
%    of the nx2 array ab.
%
function ab=chri1(N,ab0,z)
N0=size(ab0,1); if N0<N+1, error('input array ab0 too short'), end
r=zeros(N+1,1);
r(1)=z-ab0(1,1); r(2)=z-ab0(2,1)-ab0(2,2)/r(1);
ab(1,1)=ab0(2,1)+r(2)-r(1); ab(1,2)=-r(1)*ab0(1,2);
if N==1, return, end
for k=2:N
    r(k+1)=z-ab0(k+1,1)-ab0(k+1,2)/r(k);
    ab(k,1)=ab0(k+1,1)+r(k+1)-r(k);
    ab(k,2)=ab0(k,2)*r(k)/r(k-1);
end


%----------------------------------------------------------------------
% 6. "gauss"
%----------------------------------------------------------------------

% GAUSS Gauss quadrature rule.
%
%    Given a weight function w encoded by the nx2 array ab of the 
%    first n recurrence coefficients for the associated orthogonal
%    polynomials, the first column of ab containing the n alpha-
%    coefficients and the second column the n beta-coefficients, 
%    the call xw=GAUSS(n,ab) generates the nodes and weights xw of
%    the n-point Gauss quadrature rule for the weight function w.
%    The nodes, in increasing order, are stored in the first 
%    column, the n corresponding weights in the second column, of
%    the nx2 array xw.
%
function xw=gauss(N,ab)
N0=size(ab,1); if N0<N, error('input array ab too short'), end
J=zeros(N);
for n=1:N, J(n,n)=ab(n,1); end
for n=2:N
    J(n,n-1)=sqrt(ab(n,2));
    J(n-1,n)=J(n,n-1);
end
[V,D]=eig(J);
[D,I]=sort(diag(D));
V=V(:,I);
xw=[D ab(1,2)*V(1,:)'.^2];


%----------------------------------------------------------------------
% 7. "r_jacobi"
%----------------------------------------------------------------------
% R_JACOBI Recurrence coefficients for monic Jacobi polynomials.
%
%    ab=R_JACOBI(n,a,b) generates the first n recurrence 
%    coefficients for monic Jacobi polynomials with parameters 
%    a and b. These are orthogonal on [-1,1] relative to the
%    weight function w(t)=(1-t)^a(1+t)^b. The n alpha-coefficients
%    are stored in the first column, the n beta-coefficients in
%    the second column, of the nx2 array ab. The call ab=
%    R_JACOBI(n,a) is the same as ab=R_JACOBI(n,a,a) and
%    ab=R_JACOBI(n) the same as ab=R_JACOBI(n,0,0).
%
%    Supplied by Dirk Laurie, 6-22-1998; edited by Walter
%    Gautschi, 4-4-2002.
%
function ab=r_jacobi(N,a,b)
if nargin<2, a=0; end;  if nargin<3, b=a; end
if((N<=0)|(a<=-1)|(b<=-1)) error('parameter(s) out of range'), end
nu=(b-a)/(a+b+2);
mu=2^(a+b+1)*gamma(a+1)*gamma(b+1)/gamma(a+b+2);
if N==1, ab=[nu mu]; return, end 
N=N-1; n=1:N; nab=2*n+a+b;
A=[nu (b^2-a^2)*ones(1,N)./(nab.*(nab+2))];
n=2:N; nab=nab(n);
B1=4*(a+1)*(b+1)/((a+b+2)^2*(a+b+3));
B=4*(n+a).*(n+b).*n.*(n+a+b)./((nab.^2).*(nab+1).*(nab-1));
ab=[A' [mu; B1; B']];


