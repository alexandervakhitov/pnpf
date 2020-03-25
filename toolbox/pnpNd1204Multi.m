function [ Rest, test, fest, retVal] = pnpNd1204( C, V, D, tarPtNum, testvec, fans, N, stopN, pts, isFast, haveMex)    
    doTimeMes = 0;

    candsBest = [];
    
    ptDim = 3;
    kerVectDim = ptDim * tarPtNum;
    varNum = N + N*(N-1)/2;
    
    Ker = zeros(kerVectDim, N);    
    for i = 1:N
        Ker(:,i) = V(:,kerVectDim-(i-1));
    end

    eigDiag = diag(D);
    if (length(eigDiag) < kerVectDim)
        dimdiff = kerVectDim - length(eigDiag);
        eigDiag = [eigDiag; zeros(dimdiff, 1)];
    end
    eigvals = zeros(tarPtNum, 1);

    for ansind = 1:kerVectDim
        eigvals(ansind) = eigDiag(kerVectDim-ansind+1);        
    end
    
    [BP] = formBasePoints(Ker, tarPtNum);        
    
    if (tarPtNum == 3)
        alpha = 0;%0.00001;
        beta = 0;%alpha/50^2;    
    else
        alpha = 0.00001;
        beta = alpha/50^2;            
    end
    [L, R] = formLRRegularized(BP, C, N, alpha, beta, eigvals, size(pts, 2));
    
    rankBetaSqs = rank(L'*L);
    %betaSqs = pinv(L'*L)*L'*R;    
    betaSqs = pinv(L'*L)*L'*R; 
    %betaSqs = inv(L'*L)*L'*R;
    rankBetaSqs = rank(L'*L);
    [UL,SL,VL] = svd(L);
    kerVect = VL(:, end-(size(VL, 1)-rankBetaSqs)+1:end);    

    cands = [];   
    alphaMy = -1;
        
%     if (rankBetaSqs == size(kerVect, 1))
%         [vMy, alphaMy] = alphaFormula(betaSqs(1:varNum), betaSqs(varNum + 1 : 2*varNum));
%         w = vMy;
%         fGuess = sqrt(alphaMy);
%         cands = addToCands(cands, alphaMy, vMy, N, L, R);
%     end    
%          
    errSign = 0;
    if (doTimeMes)
        sv=tic
    end
    while ((size(cands, 2) == 0) && rankBetaSqs >= varNum)            
        if (rankBetaSqs == size(kerVect, 1))
            rankBetaSqs = rankBetaSqs - 1;
        end
        kerVect = VL(:, end-(size(VL, 1)-rankBetaSqs)+1:end);
        %test
        [xkNs] = adjustInLKerMulti(betaSqs, kerVect, N, L, isFast, haveMex);
%         [xkNs] = adjustInLKerMulti2(betaSqs, kerVect, N, L, isFast);
        xknLen = size(xkNs, 1)/2;
        for i = 1:size(xkNs, 2)
            [vMy, alphaMy] = alphaFormula(xkNs(1:xknLen, i), xkNs(xknLen+1:end, i));
            if (alphaMy > 0)
                cands = addToCands(cands, alphaMy, vMy, N, L, R);
            end
        end
        %[vMy alphaMy errSign] = adjustInLKer(betaSqs, kerVect, N, L);
        %[vMy alphaMy errSign] = adjustInLKer2(betaSqs, kerVect, N, L);
        
        rankBetaSqs = rankBetaSqs - 1;
    end
    if (doTimeMes)
        s1nd = toc(sv)
    end
    if (length(cands) == 0)
        [ Rest, test, fest, betafound, bguess, retVal] = generateErrorReturn();
        return;        
    end        
    resCosts = zeros(size(cands, 1), 1);
    badInd = [];
    for candsInd = 1:size(cands, 1)
        
        fGuess = cands(candsInd, N+1);        
        x = cands(candsInd, 1:N);
        if (isFast == 0)            
%         if (1)
            options=[1E-01, 1E-15, 1E-15, 1E-20, 1e-15];
            itmax = 200;
            [ret, xpopt, info1, covar]=levmar('refine_wb_res', x, R, itmax, options, fGuess, L, N);   
        else
            info1 = [0 0 ];
            ret = -1;
        end
        if (ret >= 0)
            sol1 = generateBetaSqsFromBetas(N, xpopt);
            resVal = info1(2);
        else
            sol1 = generateBetaSqsFromBetas(N, x);
            xpopt = x;            
            resVal = info1(1);
        end 

        resCosts(candsInd) = resVal;        
        if (isFast == 0)
%         if (1)
            [ret, poptf, info1, covar]=levmar('refine_f_res', fGuess, R, itmax, options, L, sol1);        
        else
            ret = 1;
            poptf = fGuess;
            info1 = ones(1,2)*norm(refine_f_res(fGuess, L, sol1)-R);
        end
        if (ret == -1)            
            poptf = fGuess;
        end
        if (ret >= 0)
            if (poptf < 50)
                badInd = [badInd; candsInd];
            end
            resCosts(candsInd) = info1(2); 
        else 
            resCosts(candsInd) = info1(1);
        end
        cands(candsInd, 1:N+1) = [xpopt poptf];
    end
    
    if (doTimeMes)
        s2nd = toc(sv)
    end
    [R, t, f, candsBest] = vote(cands, BP, C, resCosts, tarPtNum, pts);
    
    if (doTimeMes)
        s3nd = toc(sv)
    end
    
    if (f < 0 && N<stopN)
        [ Rest, test, fest, betafound, bguess, retVal] = generateErrorReturn();
        return;
    else
        if (f < 0)
            f = 50;
        end
        Rest = R;
        test = t;
        fest = f;
        betafound = [];
        bguess = []; 
        retVal = 1;        
    end
end

