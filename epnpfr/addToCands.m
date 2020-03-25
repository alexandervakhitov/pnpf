function cands = addToCands(cands, alphaMy, vMy, N, L, R)
    if (alphaMy > 0)
        fGuess = sqrt(alphaMy);
        fGuess = max([fGuess 50]);
        wguess = [vMy; alphaMy*vMy];
        %betaGuess = makeBetaGuess(wguess, N);   
        betaGuess = makeBetaGuess3(wguess, N, fGuess);   
        cands = [cands; betaGuess fGuess*ones(size(betaGuess, 1), 1)];        
%         [betaGuess ] = makeBetaGuessNL(wguess, N, L, R, fGuess);        
%         cands = [cands; betaGuess];
    end
end