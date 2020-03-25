function [f1,R1,t1] = pnpfmy(pts, Uc, tarPtNum, isFast, pnpOpts)        
%      pnpOpts.errThr = 2.5;
    
    doTimeMes = 0;
    
    Cf = zeros(tarPtNum*3, 1);
    fans = 1;
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
    if (pnpOpts.haveMex)
        [C V D] = preprocessPNP_mex(pts, Uc, tarPtNum);   
    else
        [C V D] = preprocessPNP(pts, Uc, tarPtNum);            
    end
    if (tarPtNum == 3)
        C = C(1:3, 1:3);
    end
    if (doTimeMes)
        s2 = toc(s1)
    end
    
    N = 1;
    
    currentSolution.avgErr = -1;
    currentSolution.f = [];
    currentSolution.R = [];
    currentSolution.t = []; 
%     while (N <= stopN && (retVal < 0 || avgerr > errThr || bestf > pnpOpts.fMax))   
    while (~checkPnPSolution(currentSolution, pnpOpts) && N<= stopN)
%         tic
        if (N == 1)
            [ Rest, test, fest, retVal] = pnp1(C, V, D, tarPtNum, Cf, fans, pts, isFast);        
        else
            [ Rest, test, fest, retVal] = pnpNd1204Multi( C, V, D, tarPtNum, Cf, fans, N, stopN, pts, isFast, pnpOpts.haveMex);        
        end
        if (doTimeMes)
            s5 = toc(s1)
        end        
        currentSolution = collectResult(retVal, Rest, test, fest, Uc, pts, currentSolution, isFast, pnpOpts);
        if (doTimeMes)
            s6 = toc(s1)
        end
        N = N+1;
        if (N == stopN+1)            
            if (currentSolution.avgErr < 0)
                stopN = 5;
            end
        end
    end
    
    if ((currentSolution.avgErr  < 0 || norm(currentSolution.t) >100) && isFast == 1 )
        [f1,R1,t1] = pnpfmy(pts, Uc, tarPtNum, 0, pnpOpts);
    else        
        f1 = currentSolution.f;
        R1 = currentSolution.R;
        t1 = currentSolution.t;
    end
end

function currentSolution = collectResult(retVal, Rest, test, fest, Uc, pts, currentSolution, isFast, pnpOpts)
    if (retVal > 0)        
        noDist = 0;
        [R1, t1, f1, k1est, avgerr] = refine_pos_dist(Rest, test, fest, 0, Uc, pts, 0, ...
            0, noDist, isFast, pnpOpts.fMin, pnpOpts.fMax, pnpOpts.isFastBA, pnpOpts.haveMex);
        solution.R = R1;
        solution.t = t1;
        solution.f = f1;
        solution.avgErr = avgerr;
        currentSolution = findBetterSolution(solution, currentSolution);        
    end    
end