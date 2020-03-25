function res = squareCost(fSqGuess, A, a, B, b, Q)
    eigVecGuess = findEigComb(Q, fSqGuess);
    xVec = [A a; -B -b]*eigVecGuess;
    res = [xVec(1:3)-xVec(4:6)/fSqGuess; (xVec(4)*xVec(5) - xVec(6)^2); (xVec(1)*xVec(2) - xVec(3)^2)*fSqGuess];%
end