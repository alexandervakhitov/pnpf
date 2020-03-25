function [R, t, f, candsBest] = vote(cands, BP, C, resCosts, tarPtNum, P)
    resrate = zeros(size(cands, 1), 1);
    minrate = size(P, 2);
    betaNum = size(cands, 2) - 1;
    
    candNum = size(cands, 1);
    Ta = zeros(3, 4, candNum);
    Tfa = zeros(candNum, 1);
    
    if (candNum == 0)
        f = -1;
        R = [];
        t = [];
        return;
    end
    
    validCandNum = 0;
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
        
        Ta(:, :, i) = [R t];
        Tfa(i) = fest;
        f = fest;
        Pm = [R t]*[P; ones(1, size(P, 2))];
        resrate(i) = 0;
        for ptInd = 1:size(Pm, 2)
            if (Pm(3, ptInd) < 0)
                resrate(i) = resrate(i)+1;
            end
        end
        if (resrate(i) <= minrate)
            minrate = resrate(i);
            Ra = R;
            ta = t;
            festa = fest;
        end
        validCandNum  = validCandNum +1;
    end
%     goodInds = [1:size(cands, 1)];%find(resrate == minrate);
    goodInds = find(resrate == minrate);
    res4good = resCosts(goodInds);
    [minval, minind] = min(res4good);
  
    Tta = Ta(:, :, goodInds(minind));
    candsBest = cands(goodInds(minind), 1:betaNum);
    R = Tta(1:3, 1:3);
    t = Tta(1:3, 4);
    f = Tfa(goodInds(minind));    
    if (validCandNum  == 0)
        R = [];
        t = [];
        f = -1;
    end
end