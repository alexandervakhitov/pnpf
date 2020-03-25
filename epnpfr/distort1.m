function proj = distort1(projn, k1, f) 
    r = norm(projn);
    proj = f*projn*(1+k1*r^2);    
end