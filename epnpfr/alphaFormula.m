function [vMy, alphaMy] = alphaFormula(x1, x2)
    cosG = x1'*x2 / norm(x1) / norm(x2);
    
    
    %tg(2*alpha) (r1^2 + r2^2 cos(2G)) = r2^2 sin(2G)
    if (1-cosG ^2<0)
        sinG = 0;
    else
        sinG = sqrt(1-cosG^2);    
    end
    ang = 0.5*atan2(norm(x2)^2*2*sinG*cosG, norm(x1)^2+norm(x2)^2*(2*cosG^2-1));  
    
    x2o = (x2-x2'*x1/(norm(x1)^2)*x1);
    if (norm(x2o) > 0)
        x2o = x2o / norm(x2o);
    end
    vNormed = x1/norm(x1)*cos(ang) + x2o*sin(ang);
    vMy = (x1'*vNormed)*vNormed;
    alphaMy = x2'*vMy/norm(vMy)^2;
    if (alphaMy<0)
        alphaMy;
    end
end