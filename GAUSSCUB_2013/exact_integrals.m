
function exact_result=exact_integrals(polygon_type,function_type)

% ALL THESE INTEGRALS HAS BEEN PRECOMPUTED BY HIGH DEGREE GAUSS LEGENDRE CUBATURE
% OVER POLYGONS.

switch polygon_type
    
case 1
    switch function_type
    case 1
        exact_result=0.406969589491556; 
    case 2
        exact_result=0.382597858232106; 
    case 3
        exact_result=4993.214285714256200; 
    case 4
        exact_result=1.188043774905800; 
    case 5
        exact_result=0.031415926535801; 
    case 6
        exact_result=0.000289906533545;
    case 7
        exact_result=1.000000000000000;
    end
    
case 2
    switch function_type
    case 1
        exact_result=0.260067709012133;
    case 2
        exact_result=0.156825125547501; % 10 DIGITS PRECISION.
    case 3
        exact_result=169.704343403127520; % 12 DIGITS PRECISION.
    case 4
        exact_result=0.593459365720563;
    case 5
        exact_result=0.031414528632393;
    case 6
        exact_result=0.008421180941490;
    case 7
        exact_result=0.535000000000000;
    end
    
case 3
    switch function_type
    case 1
        exact_result=0.175565704080668;
    case 2
        exact_result=0.139381456868220;   % 12 DIGITS ARE CORRECT.
    case 3
        exact_result=130.841234986796820; % 13 DIGITS ARE CORRECT.
    case 4
        exact_result=0.531841554503018;
    case 5
        exact_result=0.031220838971539;
    case 6
        exact_result=0.014222050981512;
    case 7
        exact_result=0.481250000000000;
    end
    
end