function res = LKerResMultiC(x, kerVect, lambda, initSol, varNum, L)
    nKV = size(kerVect, 2);
    sol = initSol + kerVect*x(1:nKV);
    [vMy, alphaMy] = alphaFormula(sol(1:varNum), sol(varNum+1:2*varNum));
    if (alphaMy < 0)
        alphaMy = 0;
    end
    x1est = vMy-sol(1:varNum);
    x2est = alphaMy*vMy-sol(varNum+1:2*varNum);                 
    res = L*[x1est; x2est];    
end