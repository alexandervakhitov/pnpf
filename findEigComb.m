function eigVecGuess = findEigComb(A, eigValGuess)
    %Ax = lx, l known
    %(A-lE)x = 0
    B = A - eigValGuess*eye(size(A, 1));
    [U, S, V] = svd(B);
    eigVecGuess = V(:, end);    
end