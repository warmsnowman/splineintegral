function z=fct2D(x,y,function_type)

%-------------------------------------------------------------------------------
% FUNCTION DEFINITION. 
%-------------------------------------------------------------------------------
%
% INPUT:
%
% [x,y]: EVALUATE THE FUNCTION IN THE POINTS (x,y). "x", "y" ARE COLUMN VECTORS.
%
% [function_type]: PARAMETER THAT CHOOSES THE FUNCTION.
%
% OUTPUT:
%
% [z]: VALUE OF THE FUNCTION IN THE POINTS (x,y). IT IS A COLUMN VECTOR.
%
%-------------------------------------------------------------------------------

switch function_type
case 1
    z=franke(x,y);
case 2
    z=( (x-0.5).^2 +(y-0.5).^2 ).^(1/2);
case 3
    z=(x+y).^19;
case 4
    s=(x-0.5).^2+(y-0.5).^2;
    z=exp(s);
case 5
    s=(x-0.5).^2+(y-0.5).^2; s=-100*s;
    z=exp(s);
case 6
    z=cos(30*(x+y));
case 7
    z=ones(size(x));
end

