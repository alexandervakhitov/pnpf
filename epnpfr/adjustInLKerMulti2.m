function [xkNs] = adjustInLKerMulti2(x0, kerVect, betaNum, L, isFast)
    errSign = 0;
    
    xkNs = [];
    
    xks = [];
    
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
        xkN = x0+kerVect*xkAns;
        xkNs = [xkNs xkN];
        
        xks = [xks xkAns];
        
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
            for quadInd = 1 : 2*length(workInds)
               conSign = mod(quadInd, 2);
               conInd = floor((quadInd-1) / 2)+1;
               A2 = zeros(length(workInds), length(workInds));
               b2 = zeros(length(workInds), 1);
               f = zeros(length(workInds), 1);

               if (conSign == 0)
                   f = -1;
                   A2(1, conInd) = -1;
                   if (~cnPos)
                       b2 = xnAnsBounds(workInds(conInd));
                   end
               else
                   f = 1;
                   A2(1, conInd) = 1;                       
                   if (~cnPos)
                       b2 = xnAnsBounds(workInds(conInd));
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
                           xkAns = Vq*xnAns-w;
                           xkN = x0+kerVect*xkAns;
                           xkNs = [xkNs xkN];
                           xks = [xks xkAns];
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
    end

    if (isFast == 1)
%     if (1)
        return;
    end
    xkNs = [];
    for candInd = 1:size(xks, 2)
        C = [kerVect(1:betaNum, :); kerVect(varNum+1:varNum+betaNum, :)];
        d = -[x0(1:betaNum);x0(varNum+1:varNum+betaNum)];
        imvect = zeros(size(L, 1), 1);
        xStart = xks(:, candInd);
        popt1 = xStart;
        lambda = 0;
        options=[1e-3, 1E-20, 1E-20, 1E-20, 1e-6];
        itmax = 200;
        [ret, popt1, info1, covar] = levmar('LKerResMultiC_mex', xStart, imvect, itmax, options, 'blic', ...
            -1e100*ones(size(xStart)), 1e100*ones(size(xStart)), C, d, kerVect, lambda, x0, varNum, L);
        
        if (ret < 0 || sum(C*popt1-d>-1e-3) < size(C, 1))
            C = [kerVect(1:betaNum, :)];
            d = -[x0(1:betaNum)];        
            [ret, popt1, info1, covar] = levmar('LKerResMultiC_mex', xStart, imvect, itmax, options, 'blic', ...
            -1e100*ones(size(xStart)), 1e100*ones(size(xStart)), C, d, kerVect, lambda, x0, varNum, L);
        end
        xkN = x0+kerVect*popt1;
        xkNs = [xkNs xkN];
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
