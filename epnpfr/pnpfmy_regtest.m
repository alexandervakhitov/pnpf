function [resObjNoReg resObjReg] = pnpfmy(pts, Uc, Rans, tans, fans, tarPtNum, isFast, N)        
%      pnpOpts.errThr = 2.5;
    
    doTimeMes = 0;
            
    stopN = 3;
    if (isFast)
        stopN = 2;
    end
    if (tarPtNum == 3)
        stopN = 2;
    end
    Rf = zeros(3, 3, stopN);
    tf = zeros(3, stopN);
    ff = zeros(stopN, 1);
    
    errf = zeros(stopN, 1);
    
    noDist = 0;   
    if (doTimeMes)
        s1 = tic
    end
%     [C V D] = preprocessPNP(pts, Uc, tarPtNum);            
    [C V D] = preprocessPNP_mex(pts, Uc, tarPtNum);   
    if (tarPtNum == 3)
        C = C(1:3, 1:3);
    end

    red_arr = [];
    red_ls_arr = [];
    resObjNoReg.redArr = [];
    resObjNoReg.redLsArr = [];
    resObjNoReg.rankArr = [];
    resObjNoReg.mcArr = [];
    
    resObjReg.redArr = [];
    resObjReg.redLsArr = [];
    resObjReg.rankArr = [];
    resObjReg.mcArr = [];
    
    
        
        varNum = N + N*(N-1)/2;    

        Ker = zeros(3*tarPtNum, N);

        for i = 1:1:N
            Ker(:,i) = V(:,3*tarPtNum-(i-1));
        end

        eigDiag = diag(D);
        if (length(eigDiag) < 3*tarPtNum)
            dimdiff = 12 - length(eigDiag);
            eigDiag = [eigDiag; zeros(dimdiff, 1)];
        end
        eigvals = zeros(tarPtNum, 1);

        for ansind = 1:3*tarPtNum
            eigvals(ansind) = eigDiag(3*tarPtNum-ansind+1);        
        end

        [BP] = formBasePoints(Ker, tarPtNum);        
        alpha = 1e-3;
        beta = alpha/50^2;            
        
        
        for alIt = 1:5
            alpha = 0.03*(alIt-1);
            beta = alpha/10^2;       
            kerLen = 1;
            [red red_ls rankBetaSqs mc] = testAlphaBeta(alpha, beta, BP, C, N, eigvals, pts, Rans, tans, fans, Ker, tarPtNum, Uc, kerLen);
            rls(alIt, 1) = red_ls;
            reds(alIt, 1) = red;
%             if (red_ls > 6)
%                 kerLen = 2;
%                 [red red_ls rankBetaSqs mc] = testAlphaBeta(alpha, beta, BP, C, N, eigvals, pts, Rans, tans, fans, Ker, tarPtNum, Uc, kerLen);
%                 rls(alIt, 2) = red_ls;
%                 reds(alIt, 2) = red;
%             end
%             if (red_ls > 6)
%                 kerLen = 3;
%                 [red red_ls rankBetaSqs mc] = testAlphaBeta(alpha, beta, BP, C, N, eigvals, pts, Rans, tans, fans, Ker, tarPtNum, Uc, kerLen);
%                 rls(alIt, 3) = red_ls;
%                 reds(alIt, 3) = red;
%             end
            resObjReg.redArr = [resObjReg.redArr; red];
            resObjReg.redLsArr = [resObjReg.redLsArr; red_ls];
            resObjReg.rankArr = [resObjReg.rankArr; rankBetaSqs];
            resObjReg.mcArr = [resObjReg.mcArr; mc];
        end
        

        if (min(resObjReg.redLsArr) > 7)
            resObjReg.redLsArr;
        end
        
%         alpha = 0;
%         beta = 0;            
%         [redNoReg redLsNoReg rankBetaSqsNoReg mc] = testAlphaBeta(alpha, beta, BP, C, N, eigvals, pts, Rans, tans, fans, Ker);
%         resObjNoReg.redArr = [resObjNoReg.redArr; redNoReg];
%         resObjNoReg.redLsArr = [resObjNoReg.redLsArr; redLsNoReg];
%         resObjNoReg.rankArr = [resObjNoReg.rankArr; rankBetaSqsNoReg];
%         resObjNoReg.mcArr = [resObjNoReg.mcArr; mc];
    
end

function [red red_ls rankBetaSqs minResRel] = testAlphaBeta(alpha, beta, BP, C, N, eigvals, pts, Rans, tans, fans, Ker, tarPtNum, Uc, kerLen)
    [L, R] = formLRRegularized(BP, C, N, alpha, beta, eigvals, size(pts, 2));

    rankBetaSqs = rank(L'*L);
    betaSqs = pinv(L'*L)*L'*R;    
    [UL,SL,VL] = svd(L);
    if (size(VL, 2) > 2)
        if (N==2)
            kerVect = VL(:, end-kerLen+1:end);    
        end
        if (N==3)
            kerVect = VL(:, end-10:end);    
        end
        if (N>3)
            kerVect = VL(:, end);    
        end
    else
        kerVect = VL(:, end);    
    end

    BPansh = [Rans tans; 0 0 0 1] * [C; ones(1, 4)];
    BPans = diag([1 1 1/fans^2])*BPansh(1:3, 1:4);
    BPansVec = BPans(:);
    coef0 = Ker'*BPansVec;
    red = norm(BPansVec - Ker*coef0)/norm(BPansVec);  
    coef0Sqs = coef0.^2;
    for i = 1:N
        for j = i+1:N
            coef0Sqs = [coef0Sqs; coef0(i)*coef0(j)];
        end
    end
    coef0Sqs = [coef0Sqs; fans^2*coef0Sqs];
    coef0SqsPCoef = kerVect'*(coef0Sqs - betaSqs);
    betaSqsPlusKer = kerVect*coef0SqsPCoef + betaSqs;
    red_ls = norm(coef0Sqs - betaSqsPlusKer)^2 / norm(coef0Sqs)^2;   
    red_ls1 = norm(coef0.^2 - betaSqs(1:N)) / norm(coef0.^2)
    
    [minResRel betaGuess emin] = checkEigVectors(betaSqs, kerVect, coef0Sqs, C, BP, tarPtNum, Uc, pts);
    
    for i = 1:size(betaGuess, 1)
        bgn(i) = norm(betaGuess(i, :)' - coef0)/norm(coef0);
    end
    
%     red_ls = min(bgn);
    
%     if (red_ls > 0.25 && alpha > 0 && N == 3)
%         mc;
%     end
%     
%     red_ls = emin;
%     
%     if (emin > 1 && alpha>0)
%         emin;
%     end

end
