function [fFin RFin tFin goodInds] = pnpfmyRANSAC(p3d, p2d, Ng, thr)
    n = size(p3d, 2);
    tarPtNum = 4;
    isFast = 1;
        
    
    pnpOpts.fMax = 2e4;
    pnpOpts.fMin = 30;
    pnpOpts.errThr = 7.5;
    pnpOpts.isFastBA = 0;
    pnpOpts.haveMex = 0;

    
    indx = zeros(6,1);  
    maxInliers = 0;
    for s = 1:Ng
        j = 1;
        taken = zeros(n,1);
        while j <= 6
            rnd = rand(1,1,'double');
            c = floor(rnd*n)+1;
            if (~taken(c))
                p3dc(1:3, j) = p3d(1:3, c);
                p2dc(1:2, j) = p2d(1:2, c);
                indx(j) = c;
                taken(c) = 1;                
                j = j+1;
            end
        end            
        
        [f1,R1,t1] = pnpfmy(p3dc, p2dc, tarPtNum, isFast, pnpOpts);
        if (length(f1) > 0)
%         [f1,R1,t1] = GPnP_f_GN(p3dc, p2dc);
            K = [f1 0 0; 0 f1 0; 0 0 1];
            P = K*[R1 t1];
            inliers = 0;
            for ptInd = 1:n            
                [dist1,proj1] = triang.findProjectionError(P, p3d(1:3, ptInd), p2d(1:2, ptInd));
                if (dist1 < thr)
                    inliers = inliers + 1;
                end
            end
        else
            f1;
        end
        if (inliers > maxInliers)
            maxInliers = inliers;
            maxR = R1;
            maxt = t1;
            maxf = f1;
        end
        s
    end
    goodInds = [];
    RFin = [];
    tFin = [];
    fFin = -1;
    if (maxInliers > 6)
        K = [maxf 0 0; 0 maxf 0; 0 0 1];
        P = K*[maxR maxt];
        inliers = 0;        
        for ptInd = 1:n            
            [dist1,proj1] = triang.findProjectionError(P, p3d(1:3, ptInd), p2d(1:2, ptInd));
            if (dist1 < thr)
                goodInds = [goodInds ptInd];
            end
        end
        isFast = 1;
        [fFin,RFin,tFin] = pnpfmy(p3d(:, goodInds), p2d(:, goodInds), tarPtNum, isFast, pnpOpts);
    end
end