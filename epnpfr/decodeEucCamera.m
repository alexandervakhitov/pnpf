function [R,t]=decodeEucCamera(camparam)
    rotpars = zeros(3,1);
    rotpars(:) = camparam(1:3);
    R=rodrigues(rotpars);
    t=zeros(3,1);
    t(1:3)=camparam(4:6);    
end