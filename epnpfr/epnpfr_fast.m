function  [f1,R1,t1] = epnpfr_fast(pts, Uc)
    isFast = 1;
    pnpOpts.fMax = 1e4;
    pnpOpts.fMin = 30;
    pnpOpts.errThr = 7.5;
    pnpOpts.isFastBA = 1;
    tarPtNum = 4;
    [f1,R1,t1] = pnpfmy(pts, Uc, tarPtNum, isFast, pnpOpts);
end