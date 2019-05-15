
clear all; clear plot;

fprintf('\n \t [GAUSS LIKE FORMULA FOR POLYGONS][VERS. 4.0] [MARCH 7, 2006] \n');

%--------------------------------------------------------------------------
% REFERENCE PAPER 
% [1] A. SOMMARIVA and M. VIANELLO 
%     "Gauss-like and triangulation-free cubature over polygons".
%
% THIS DEMO SHOWS HOW TO USE THE PROCEDURE "polygauss" FOR COMPUTING CUBATU-
% RE OVER POLYGONS OF CONTINUOUS FUNCTIONS.
%
%--------------------------------------------------------------------------
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
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% SEVERAL POLYGONS ARE DEFINED IN THE m-FILE "define_polygon.m". 
%
% [polygon_type=1]: UNIT SQUARE [0,1]^2.
% [polygon_type=2]: A 6 SIDES CONVEX POLYGON USED IN PAPER [1].
% [polygon_type=3]: A 9 SIDES NONCONVEX POLYGON USED IN THE PAPER [1].
%
% TO DEFINE YOUR OWN POLYGON, PLEASE OPEN THE m-FILE "define_polygon.m" 
% AND ADD IT TO THE LIST. OBSERVE THAT THE FIRST VERTEX IS NOT REPEATED. 
% FOR EXAMPLE, THE UNIT SQUARE IS DEFINED AS "[0 0; 1 0; 1 1; 0 1]" 
% AND NOT "[0 0; 1 0; 1 1; 0 1; 0 0]".
%--------------------------------------------------------------------------
polygon_type=2;

%--------------------------------------------------------------------------
% SEVERAL FUNCTIONS ARE DEFINED IN THE m-FILE "fct2D.m". MORE PRECISELY
%
% [function_type=1]: z=franke(x,y);
% [function_type=2]: z=( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
% [function_type=3]: z=(x+y).^19;
% [function_type=4]: z=exp(loc_arg); loc_arg=(x-0.5).^2+(y-0.5).^2; 
% [function_type=5]: z=exp(log_arg); loc_arg=-100*((x-0.5).^2+(y-0.5).^2); 
% [function_type=6]: z=cos(30*(x+y));
% [function_type=7]: z=ones(size(x));
%
% TO DEFINE YOUR OWN FUNCTION, PLEASE OPEN THE m-FILE "fct2D.m" AND ADD IT 
% TO THE LIST. 
%--------------------------------------------------------------------------
function_type=5;

%--------------------------------------------------------------------------
% WRITE THE DEGREE OF THE 1-DIMENSIONAL GAUSS-LEGENDRE RULE.
%--------------------------------------------------------------------------
N=10; 

%--------------------------------------------------------------------------
% FOR THE DESCRIPTION OF HOW A ROTATION WORKS SEE PAPER [1].
%
% SUPPOSE x_min, x_max, y_min, y_max ARE THE MINIMUM AND MAXIMUM VALUES IN x
% AND y REACHED BY THE POLYGON, I.E. "R=[x_min,x_max] x [y_min,y_max]" IS THE  
% SMALLEST RECTANGLE WITH SIDES PARALLEL TO THE AXIS x AND y, CONTAINING THE
% POLYGON.
%
% [rotation=0]: GENERAL CASE, BUT OBSERVE THAT THE INTEGRAND MUST BE DEFINED
%               IN THE RECTANGLE "R" DESCRIBED ABOVE.
% [rotation=1]: GOOD FOR CONVEX POLYGONS. FOR CONVEX POLYGONS IT SUFFICES 
%               THAT THE INTEGRAND BE DEFINED IN THE POLYGON.
% [rotation=2]: THE USER CHOOSES A SPECIAL REFERENCE SEGMENT "PQ". CHECK [1]  
%               FOR FURTHER DETAILS. SEE THE VARIABLES "P", "Q" BELOW.
%--------------------------------------------------------------------------
rotation=1; 

%--------------------------------------------------------------------------
% IF [rotation = 2] THEN THE ALGORITHM CHOOSES A PREFERRED SEGMENT "PQ" 
% HAVING "P", "Q" AS EXTREMA.
%--------------------------------------------------------------------------
P=[0 0.75]; Q=[1 0.5];

%--------------------------------------------------------------------------
% [plot_polygon = 0] NO PLOTS.
% [plot_polygon = 1] IT PLOTS THE DOMAIN AND THE CUBATURE NODES.                   
%--------------------------------------------------------------------------
plot_polygon = 1;

%--------------------------------------------------------------------------
% [cubature_type]: IT CHOOSES THE 1D QUADRATURE RULE .
%
%           [cubature_type=1]: FEJER 1.
%           [cubature_type=2]: FEJER 2.
%           [cubature_type=3]: CLENSHAW CURTIS (VIA WEIGHTS).
%           [cubature_type=4]: GAUSS-LEGENDRE.
%--------------------------------------------------------------------------
cubature_type=4;


%--------------------------------------------------------------------------
%                    THE MAIN PROGRAM STARTS HERE.
%--------------------------------------------------------------------------

%--------------------------------------------------------------------------
% DEFINE POLYGON.
%--------------------------------------------------------------------------
% "define_polygon.m" DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT
% THE END THE FIRST POINT). THE VARIABLE "polygon_type" IS DEFINED ABOVE. 
% "polygon_sides" IS A VECTOR HAVING AS COMPONENTS THE VERTICES OF THE 
% COMPONENT, DESCRIBED COUNTERCLOCKWISE. 
%--------------------------------------------------------------------------
polygon_sides=define_polygon(polygon_type);

%--------------------------------------------------------------------------
% COMPUTING INTEGRAL 
%--------------------------------------------------------------------------
% "cubature_result", NODES "(x_nodes,y_nodes)" AND WEIGHTS "weights".
% "fct2D.m" IS AN m-file DEFINED BY THE USER, HAVING "function_type" AS
% VARIABLE THAT SWITCHES BETWEEN FUNCTIONS. SEE THE FILE "fct2D.m" FOR 
% DETAILS.
%--------------------------------------------------------------------------
cputime1=cputime;
[cubature_result,x_nodes,y_nodes,weights]=polygauss(@fct2D,N,polygon_sides,rotation,P,Q,[],cubature_type,function_type);
cputime2=cputime;

%--------------------------------------------------------------------------
% SOME STATS.
%--------------------------------------------------------------------------
fprintf('\n \t ----------------------------------------------------------------------------');
fprintf('\n \t [SIDES   ]: %2.0f [CUBATURE TYPE]: %2.0f [DEGREE ]: %5.0f [PTS]: %5.0f',...
    length(polygon_sides(:,1))-1,cubature_type, N,size(x_nodes,1)*size(x_nodes,2));
fprintf('\n \t [DOMAIN  ]: %2.0f [FUNCTION     ]: %2.0f [CPUTIME]: %2.2e',polygon_type,function_type,cputime2-cputime1);

fprintf('\n \t [PGAUSS RESULT ]: %5.15f',cubature_result); 
exact_result=exact_integrals(polygon_type,function_type);
fprintf('\n \t [EXACT RESULT  ]: %5.15f',exact_result); 
fprintf('\n \t [ABSOLUTE ERROR]: %5.15e',abs(cubature_result-exact_result)); 
if abs(exact_result) > 0
    relerr=abs(cubature_result-exact_result)/abs(exact_result);
    fprintf('\n \t [RELATIVE ERROR]: %5.15e',relerr); 
end
fprintf('\n \t ----------------------------------------------------------------------------');

%--------------------------------------------------------------------------
% PLOT POLYGONS.
%--------------------------------------------------------------------------
if plot_polygon == 1        
    fill(polygon_sides(:,1),polygon_sides(:,2),'w'); hold on; 
    plot(x_nodes,y_nodes,'k.'); 
end
 


