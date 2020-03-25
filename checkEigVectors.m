function [minResRel betaGuess emin] = checkEigVectors(x0, kerVect, xTrue, C, BP, tarPtNum, Uc, pts)    
    v = size(kerVect, 1);    
    A = kerVect(1:v/2, :);
    a = x0(1:v/2);
    B = -kerVect(v/2+1:end, :);
    b = -x0(v/2+1:end);
    %l(Ax + a) = Bx + b
    normCoef = 1e0;
    Abig = normCoef*[A a; zeros(1, size(A, 2)) 0];
    Bbig = [B b; zeros(1, size(B, 2)) 0];
    Qbig2 = pinv(Abig)*Bbig;
    [V, D] = eig(Qbig2);
    
    negSqFvars = diag(D);
    
    fTrueSq = norm(xTrue(v/2+1:end))/norm(xTrue(1:v/2));
    
    negSqFvars1  = abs(negSqFvars*normCoef)  - fTrueSq ;
%     qqFvars  = -negSqFvars  + fTrueSq ;
    
    minResRel = min(abs(negSqFvars1  )) / fTrueSq;
    
    betaGuess = [];
    cands = [];
    Qbig2
    negSqFvars
    if (sum(abs(imag(negSqFvars))) > 0)
        negSqFvars = real(negSqFvars );
        V = real(V);
    end
    for i = 1:length(negSqFvars)
        evGuess = negSqFvars(i);
        eVect = V(:, i);
        eVect = eVect  / eVect (2);        
        if (evGuess > 0)
            evGuess = -evGuess;
            eVect = - eVect;            
        end
        xGuess = [A a; B b]*eVect;
        fGuess = sqrt(-evGuess);        
        bg1 = makeBetaGuess3Z(xGuess, 2, fGuess);
        betaGuess = [betaGuess; bg1];
        cands = [cands; bg1 fGuess*ones(size(bg1, 1), 1)];
    end
    betaNum = size(cands, 2)-1;
    errs = [];
    for i = 1:size(cands, 1)
%         resrate(i) = size(P, 2);
        fest = cands(i, betaNum+1);
        betaVect = zeros(betaNum, 1);
        for betaInd = 1:betaNum
%             betaVect(betaInd) = cands(i, betaInd)*fest;
            betaVect(betaInd) = cands(i, betaInd);
        end        
        BPans = formBPAns(BP, betaVect, fest, tarPtNum);
        if (tarPtNum == 4)
            [errFlag, R, t] = rtFromC(C, BPans, tarPtNum);
            
            noDist = 0;
            isFastBA = 1;
            isFast = 1;
            if (errFlag == 1)
                continue;
            end
            if (sum(isnan(R(:))) > 0)
                continue;
            end
            cands
            [Rest, test, fest, dest, avgerr] = refine_pos_dist(R, t, fest, 0, Uc, pts, 0, 0, ...
                noDist, isFast, 10, 1e5, isFastBA);
            errs = [errs; avgerr];
            if (errFlag == 1)
                continue;
            end
        else
            BPans = formBPAns(BP, cands(i, 1:betaNum), fest, tarPtNum);            
            F = []; 
            G = [];
            for ptInd = 1:3
                for ptInd2 = ptInd+1:3
                    fcol = C(:, ptInd) - C(:, ptInd2);
                    gcol = BPans(:, ptInd) - BPans(:, ptInd2);
                    F = [F fcol];
                    G = [G gcol];                                            
                end
            end  
            F(:, 3) = cross(F(:, 1), F(:, 2));
            G(:, 3) = cross(G(:, 1), G(:, 2));
            R = G/F;
            t = BPans(:, 1) - R*C(:, 1);
        end
    end
    
    if (size(errs, 1) == 0)
        errs;
    end
    emin = min(errs);
%     eigValGuess = zeros(size(negSqFvars));
%     for eigInd = 1:length(negSqFvars)
%         itmax = 200;
%         options=[1E-03, 1E-15, 1E-15, 1E-15, 1e-9];
%         x00 = negSqFvars(eigInd);
%         imvect = zeros(5, 1);
%         [ret, popt, info, covar]=levmar('squareCost', x00, imvect, itmax, options, A, a, B, b, Qbig2);    
%         eigValGuess(eigInd) = popt;
%     end
%     
%     negSqFvarsNL  = abs(eigValGuess)  - fTrueSq ;
% %     qqFvars  = -negSqFvars  + fTrueSq ;
%     
%     minResRelNL = min(abs(negSqFvarsNL  )) / fTrueSq;
    
    
    coefs = [];
    
    zeroThr = 1;
    zeroEV = [];
    xCoefs = [];
    for j = 1:size(kerVect, 2)
        xCoefs  = [xCoefs; kerVect(:, j)'*(xTrue - x0)];
    end
    for i = 1:size(V, 1)
        if (abs(D(i, i)) < zeroThr)
            zeroEV = [zeroEV V(:, i)];
        end
    end
    ds = size(kerVect, 1)/2;
    for i = 1:size(V, 1)
        if (abs(real(D(i,i))) > zeroThr)
            Vm = [zeroEV V(:, i)];
            xcEst = Vm*(pinv(Vm'*Vm)*Vm'*[xCoefs; 1]);                        
%             coefs = [coefs; norm(xTrue - fullAns)/norm(xTrue)];                                    
            xEst = [kerVect x0]*xcEst;
            coefs = [coefs; norm(xEst(ds+1:2*ds)-xTrue(ds+1:2*ds))/norm(xTrue(ds+1:2*ds))];                        
        end
    end
    mc = min(coefs);
    if (mc > 0.1)
        mc;
        %todo
        %add eigvects for zero eigvals to system
    end
end