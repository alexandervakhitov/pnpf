function [vMy alphaMy errSign popt] = adjustInLKer(x0, kerVect, betaNum, L, isFast)
    errSign = 0;
    
    nKV = size(kerVect, 2);
    
    if (nKV > 1)
       nKV; 
    end
    Q = zeros(nKV, nKV);
    b = zeros(nKV, 1);       
    
    varNum = betaNum + betaNum*(betaNum-1)/2;
    x1 = x0(1:varNum);
    x2 = x0(varNum+1:2*varNum);
    c = x1'*x2;
    for kerInd = 1:nKV
        b(kerInd) = x1'*kerVect(varNum+1:2*varNum, kerInd) + kerVect(1:varNum, kerInd)'*x2;
    end
    for kerInd1 = 1:nKV
        for kerInd2 = 1:nKV
            Q(kerInd1, kerInd2) = 0.5*(kerVect(1:varNum, kerInd1)'*kerVect(varNum+1:2*varNum, kerInd2)+kerVect(1:varNum, kerInd2)'*kerVect(varNum+1:2*varNum, kerInd1));
        end
    end
%     if (rank(Q) < nKV)
%         err;
%     end
    [Vq, Sq] = eig(Q);    
    if (nKV == 1)
        Vq = 1;
        Sq = Q;
    end
    if (nKV == 2)
        nKV;
    end
    %transform to diagonalize is x_new = Vq' ( x_old + w)
    %x_old = Vq * x_new - w
    %inequality transformed to x_new' Sq x_new + cn > 0
    
    %x0 + k'*xo>0
    %x0+k'*Vq*xn-k'*w>0 => k'*Vq*xn>k'*w-x0
    w = 0.5*pinv(Q)*b;
    cn = -0.25*b'*pinv(Q)*b+c;

    svect = diag(Sq);
    
    [Ab, bb, AbF, bbF] = formConstraints(betaNum,  Vq, kerVect, x0, w, varNum);
    
    if (sum(sign(diag(Sq))) == -nKV)        
        
        xkAns = -w;        
        %errSign = 1;
        %do whatever possible....
    else
        smallStep = 1e-1;
        
        posSVect = (svect > 0);
        xnAns = zeros(size(svect));
        cnPos = 0;
        xnAnsBounds = xnAns;
        if (cn > 0)            
            xnAnsBounds = zeros(size(xnAnsBounds));
            cnPos = 1;
        else
            for diagAnsComp = 1:length(svect)
                if (svect(diagAnsComp) > 0)
                    xnAnsBounds(diagAnsComp) = sqrt((-cn)/svect(diagAnsComp)) + smallStep;
                end
            end
        end             
        
        
        workInds = find(posSVect);
        ineqHoldInds = (Ab*xnAns > bb);
        if (ineqHoldInds < size(Ab, 1))
            noGoodAns = 1;
            for quadInd = 0 : 2^length(workInds)-1
               binVec = de2bi(quadInd, length(workInds));
               A2 = zeros(length(workInds), length(workInds));
               b2 = zeros(length(workInds), 1);
               f = zeros(length(workInds), 1);
               for varInd = 1:length(workInds)
                   if (binVec(varInd) == 0)
                       f(varInd) = -1;
                       A2(varInd, varInd) = -1;
                       if (~cnPos)
                           b2(varInd) = xnAnsBounds(workInds(varInd));
                       end
                   else
                       f(varInd) = 1;
                       A2(varInd, varInd) = 1;                       
                       if (~cnPos)
                           b2(varInd) = xnAnsBounds(workInds(varInd));
                       end
                   end
               end
               for conInd = 1:size(Ab, 1)
                   Ab2 = -[Ab(conInd, workInds); AbF(conInd, workInds); A2];
                   bb2 = -[bb(conInd); bbF(conInd); b2];
                   options = optimset('LargeScale', 'off');
                   [xl,fval,exitflag,output,lambda] = linprog(f, Ab2, bb2, [], [], [], [],[], options);
                   if (exitflag == 1 || exitflag == -3 || exitflag == -7)
                       constrIds = (Ab2*xl-bb2<=0);
                       if (sum(constrIds) == length(constrIds))
                           xnAns(workInds) = xl;                       
                           noGoodAns = 0;
                           break;
                       end
                   end
               end
               if (~noGoodAns)
                   break;
               end
            end
            if (noGoodAns)
                vMy = [];
                alphaMy = -1;
                popt = [];
                return;
            end
        end        
        
        xkAns = Vq*xnAns-w;
    end
    xkN = x0+kerVect*xkAns;
    
    [vNMy, alphaNMy] = alphaFormula(xkN(1:varNum), xkN(varNum+1:2*varNum));
    if (alphaNMy < 0)
        vMy = vNMy;
        alphaMy = alphaNMy;
        popt = [];
        return;
    end
    if (isFast == 1)
        vMy = vNMy;
        alphaMy = alphaNMy;
        popt = [];
        return;
    end
        
    kerCoef = xkAns;

    options=[1e-3, 1E-20, 1E-20, 1E-20, 1e-6];
    itmax = 200;        
    lagrVar = 0;
    addVarAlpha = sqrt(alphaNMy);
    addVars = zeros(betaNum, 1);
    
    
    for betaInd = 1:betaNum
        if (xkN(betaInd) > 0)
            addVars(betaInd) = sqrt(xkN(betaInd));
        else
            addVars(betaInd) = 0.001;            
        end
    end
    xStart = [kerCoef; addVarAlpha; addVars];
    lambda = 1e10;
    
    [vMyStart, alphaMyStart] = alphaFormula(xkN(1:varNum), xkN(varNum+1:2*varNum));
    
    
%     C = [kerVect(1:betaNum, :); kerVect(varNum+1:varNum+betaNum, :)];
%     d = -[x0(1:betaNum);x0(varNum+1:varNum+betaNum)];
%     imvect = zeros(size(L, 1), 1);
    xStart = [kerCoef];
    popt1 = xStart;
%     if (isFast == 0)
%         [ret, popt1, info1, covar] = levmar('LKerResMultiC_mex', xStart, imvect, itmax, options, 'blic', ...
%             -1e100*ones(size(xStart)), 1e100*ones(size(xStart)), C, d, kerVect, lambda, x0, varNum, L);
% 
%         if (ret < 0 || sum(C*popt1-d>-1e-3) < size(C, 1))
%             C = [kerVect(1:betaNum, :)];
%             d = -[x0(1:betaNum)];        
%             [ret, popt1, info1, covar] = levmar('LKerResMultiC_mex', xStart, imvect, itmax, options, 'blic', ...
%             -1e100*ones(size(xStart)), 1e100*ones(size(xStart)), C, d, kerVect, lambda, x0, varNum, L);
%         end
%     end
    popt = popt1;
        
    xAns = x0 + kerVect*popt(1:nKV);
    [vMy, alphaMy] = alphaFormula(xAns(1:varNum), xAns(varNum+1:2*varNum));
    if (alphaMy < 0)
        vMy = vMyStart;
        alphaMy = alphaMyStart;
    end
end
function [Ab, bb, AbF, bbF] = formConstraints(betaNum, Vq, kerVect, x0, w, varNum)
%get bounds from constraints on variables equal to squares
    Ab = zeros(betaNum, size(Vq, 2));
    bb = zeros(betaNum, 1);
    AbF = Ab;
    bbF = bb;
    for betaInd = 1:betaNum        
        Ab(betaInd, :) = kerVect(betaInd, :)*Vq;        
        bb(betaInd) = kerVect(betaInd, :)*w - x0(betaInd);
        AbF(betaInd, :) = kerVect(varNum + betaInd, :)*Vq;        
        bbF(betaInd) = kerVect(varNum + betaInd, :)*w - x0(varNum + betaInd);        
    end
end
