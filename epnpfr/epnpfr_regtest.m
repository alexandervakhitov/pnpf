function  [resObjNoReg resObjReg] = epnpfr_regtest(pts, Uc, R, t, f, N)
    isFast = 1;    
    tarPtNum = 4;
    [resObjNoReg resObjReg] = pnpfmy_regtest(pts, Uc, R, t, f, tarPtNum, isFast, N);
end