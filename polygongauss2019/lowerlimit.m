function x = lowerlimit(y,poly,dchain )
x=interp1(poly(dchain(1,1):dchain(1,2),2),poly(dchain(1,1):dchain(1,2),1),y);
end

