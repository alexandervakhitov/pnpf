function [ Rest, test, fest, retVal] = pnp1(C, V, D, tarPtNum, testvec, fans, pts, isFast)
    
    
    N = 1;
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
    if (tarPtNum == 3)
        alpha = 0;%0.00001;
        beta = 0;%alpha/50^2;    
    else
        alpha = 0.00001;
        beta = alpha/50^2;            
    end
    
    [L, R] = formLRRegularized(BP, C, N, alpha, beta, eigvals, size(pts, 2));
    
    rankBetaSqs = rank(L'*L);
    betaSqs = pinv(L'*L)*L'*R;    
    [UL,SL,VL] = svd(L);
    kerVect = VL(:, end-(size(VL, 1)-rankBetaSqs)+1:end);    

    cands = [];
    errSign = 0;
    alphaMy = -1;
    if (rankBetaSqs == length(betaSqs))
        [vMy, alphaMy] = alphaFormula(betaSqs(1:varNum), betaSqs(varNum + 1 : 2*varNum));
        cands = addToCands(cands, alphaMy, vMy, N, L, R);
    end    
         
    while ((alphaMy < 0 || sum(vMy(1:N) > 0) ==0 || errSign == 1) && rankBetaSqs >= varNum)            
        if (rankBetaSqs == size(kerVect, 1))
            rankBetaSqs = rankBetaSqs - 1;
        end
        kerVect = VL(:, end-(size(VL, 1)-rankBetaSqs)+1:end);
        [vMy alphaMy errSign] = adjustInLKer(betaSqs, kerVect, N, L, isFast);
        if (alphaMy > 0)
            cands = addToCands(cands, alphaMy, vMy, N, L, R);
        end
        rankBetaSqs = rankBetaSqs - 1;        
    end
    if (alphaMy < 0)
        [ Rest, test, fest, betafound, bguess, retVal] = generateErrorReturn();
        return;
    end        

% 
%     if (betaSqs(1) < 0 || betaSqs(2) < 0)
%         Rest = [];
%         test = [];
%         fest = [];
%         betafound = [];
%         bguess = []; 
%         retVal = -1;
%         return;
%     else
%         fGuess = sqrt(betaSqs(2) / betaSqs(1));
%         fGuess = max([fGuess 50]);
%         cands = [cands; sqrt(betaSqs(2)) fGuess; -sqrt(betaSqs(2)) fGuess];
%     end
    resCosts = zeros(size(cands, 1), 1);
    for candsInd = 1:size(cands, 1)
        x = cands(candsInd, 1);
        fGuess = cands(candsInd, 2);
        if (isFast == 0)
            options=[1E-01, 1E-15, 1E-15, 1E-20, 1e-15];
            itmax = 200;               
            [ret, xpopt, info1, covar]=levmar('refine_wb_res', x, R, itmax, options, fGuess, L, N);   
        else
            ret = -1;        
        end
        if (ret >= 0)
            sol1 = xpopt^2;
        else
            sol1 = x^2;
            xpopt = x;
        end
        if (isFast == 0)
%         if(1)
            options=[1E-01, 1E-15, 1E-15, 1E-20, 1e-15];
            itmax = 200;                           
            [ret, poptf, info1, covar]=levmar('refine_f_res', fGuess, R, itmax, options, L, sol1);
        else
            info1 = ones(1,2)*norm(refine_f_res(fGuess, L, sol1)-R);
            ret = 1;
            poptf = fGuess;
        end
        if (ret >= 0)
            cands(candsInd, :) = [xpopt' poptf];
            resCosts(candsInd) = info1(2); 
        else 
            resCosts(candsInd) = info1(1);
        end
    end
    
    [R, t, f, candsBest] = vote(cands, BP, C, resCosts, tarPtNum, pts);
    if (f < 0)
        Rest = [];
        test = [];
        fest = [];
        betafound = [];
        bguess = []; 
        retVal = -1;
    else
        Rest = R;
        test = t;
        fest = f;
        betafound = [];
        bguess = []; 
        retVal = 1;        
    end
end