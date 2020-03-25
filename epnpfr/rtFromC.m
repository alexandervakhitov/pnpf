function [errFlag, R, t] = rtFromC(C, BPans, tarPtNum)
    Ch = [C; ones(1, tarPtNum)];
    errFlag = 0;
    BPansh = [BPans; ones(1, tarPtNum)];
    Rth = BPansh / Ch;
    if (det(Rth(1:3, 1:3)) < 0)
        Rth = -Rth;
    end
    if (det(Rth(1:3, 1:3)) == 0)
        errFlag = 1;
    end
    Rth = Rth / (det(Rth(1:3, 1:3)))^(1/3);
    Rt = Rth(1:3, 1:4);
    R = Rt(:,1:3);
    t = Rt(:,4);

end