
function polygon_sides=define_polygon(polygon_type)

%-------------------------------------------------------------------------------
% POLYGON DEFINITION. 
%-------------------------------------------------------------------------------
%
% INPUT:
%
% [polygon_type]: PARAMETER THAT CHOOSES THE POLYGON. 
%
% OUTPUT:
% [polygon_sides]: DEFINES THE COORDINATES OF THE VERTICES (REPEATING AT THE END 
%                  THE FIRST POINT). THE VARIABLE "polygon_type" IS DEFINED ABOVE;
%                  "boundary_pts" IS A VECTOR HAVING AS COMPONENTS THE
%                  VERTICES OF THE COMPONENT, DESCRIBED COUNTERCLOCKWISE. 
%                  
%-------------------------------------------------------------------------------
% OBSERVE THAT THE FIRST VERTEX IS NOT REPEATED. 
% FOR EXAMPLE, THE UNIT SQUARE IS DEFINED AS "[0 0; 1 0; 1 1; 0 1]" 
% AND NOT "[0 0; 1 0; 1 1; 0 1; 0 0]".
%-------------------------------------------------------------------------------

switch polygon_type
    
case 1
    fprintf('\n \t [POLYGON]: UNIT SQUARE [0,1]^2'); 
    polygon_sides=[0 0; 1 0; 1 1; 0 1]; 
case 2
    fprintf('\n \t [POLYGON]: CONVEX POLYGON'); 
    polygon_sides=[0.1 0; 0.7 0.2; 1 0.5; 0.75 0.85; 0.5 1; 0 0.25];
case 3    
    fprintf('\n \t [POLYGON]: NON CONVEX POLYGON'); 
    polygon_sides=(1/4)*[1 0; 3 2; 3 0; 4 2; 3 3; 3 0.85*4; 2 4; 0 3; 1 2];
end

%-------------------------------------------------------------------------------
% "CLOSE" THE POLYGON (IT IS NEEDED BY INPOLYGON AND BY DBLQUAD (FOR TESTS)).
%-------------------------------------------------------------------------------
polygon_sides=[polygon_sides; polygon_sides(1,:)];