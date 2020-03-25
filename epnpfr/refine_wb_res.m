function res = refine_wb_res(x, fGuess, L, betaNum)    
    varNum = betaNum + betaNum*(betaNum-1)/2;
    sol1 = zeros(varNum, 1);
    ind = betaNum+1;
    for i1 = 1:betaNum
        sol1(i1) = x(i1)^2;
        for i2 = i1+1:betaNum
            sol1(ind) = x(i1)*x(i2);
            ind = ind+1;
        end
    end
    sol = [sol1; fGuess^2*sol1];
    res = L*sol;
end