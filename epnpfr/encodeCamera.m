function camparam=encodeCamera(R,t)  %#codegen
    camparam=zeros(6,1);
    camparam(1:3)=rotMatToVect(R);    
    camparam(4:6)=t(1:3);
end