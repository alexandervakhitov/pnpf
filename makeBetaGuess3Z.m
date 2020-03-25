function betaGuess = makeBetaGuess3Z(wguess, betaNum, fGuess)    
    betaGuess = [];
    posInds = find(wguess(1:betaNum) > 0);
    posIndNum = length(posInds);
    if (posIndNum == 0)
        wguess = -wguess;
        posInds = find(wguess(1:betaNum) > 0);
        posIndNum = length(posInds);
%         
%         return;
    end    
    for i = 1:2^posIndNum
        binVec = de2bi(i-1, posIndNum);
        row = zeros(1, betaNum);
        for varInd = 1:posIndNum
            if (wguess(posInds(varInd)) < 0)
                row(posInds(varInd)) = 0;
            else
                row(posInds(varInd)) = (-1)^binVec(varInd)*sqrt(wguess(posInds(varInd)));
            end
        end
        betaGuess = [betaGuess; row];
    end
    
    if (betaNum > 1)
        dists = zeros(size(betaGuess, 1), 1);
        for i = 1:size(betaGuess, 1)
            sol1 = generateBetaSqsFromBetas(betaNum, betaGuess(i, 1:betaNum));
            dists(i) = norm([sol1; fGuess^2*sol1] - wguess);
        end
        [minVal] = min(dists);
        distInds = find(abs(dists-minVal) < 1e-10);
        betaGuess = betaGuess(distInds, :);
    end
end