function [L, R] = formLRRegularized(BP, C, N, alpha, beta, eigvals, ptNum)
    [L, R] = formLR(BP, C, N);
    varNum = N + N*(N-1)/2;
%     varNum = N;
    L2 = zeros(2*varNum, size(L, 2));
    R2 = zeros(2*varNum, 1);    
    rowInd = 1;
    for varInd = 1:N
        L2(rowInd, varInd) = alpha*eigvals(varInd)^2;
%         L2(rowInd+varNum, varNum + varInd) = beta*eigvals(varInd)^2;
        rowInd = rowInd+1;
    end
    colInd = N+1;
    for i1 = 1:N
        for i2 = i1+1:N
%             L2(rowInd, colInd) = alpha*eigvals(i1)^2*eigvals(i2)^2;
%             L2(rowInd+varNum, colInd+varNum) = beta*eigvals(i1)^2*eigvals(i2)^2;
            L2(rowInd, colInd) = alpha*abs(eigvals(i1)*eigvals(i2));
%             L2(rowInd+varNum, colInd+varNum) = beta*abs(eigvals(i1)*eigvals(i2));

            rowInd = rowInd+1;
            colInd = colInd+1;
        end
    end
        
        
%     L = [1e3*L; L2];
%     R = [1e3*R; R2];
%     L = [ptNum*L; L2];
%     R = [ptNum*R; R2];
    if (alpha > 1e-10)
        L = [L; L2];
        R = [R; R2];
    end

end
