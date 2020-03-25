function [L, R] = formLR(BP, C, N)
    maxptnum = size(BP, 2);
    cnt = 1;
    eqnum = floor(maxptnum*(maxptnum-1)/2);
    varNum = 2*(N+N*(N-1)/2);
    L = zeros(eqnum, varNum);
    R = zeros(eqnum, 1);
    cnt = 1;

    for i = 1:maxptnum-1
        for j = i+1:maxptnum  
            colCnt = 1;
            for bpInd = 1:N                                
                L(cnt, colCnt) = norm(BP(1:2, i, bpInd) - BP(1:2, j, bpInd))^2;
                colCnt = colCnt+1;
            end
            for bpInd1 = 1:N
                for bpInd2 = bpInd1+1:N
                    L(cnt, colCnt) = 2*(BP(1:2, i, bpInd1) - BP(1:2, j, bpInd1))'*(BP(1:2, i, bpInd2) - BP(1:2, j, bpInd2));            
                    colCnt = colCnt+1;
                end
            end
            for bpInd = 1:N
                L(cnt, colCnt) = (BP(3, i, bpInd) - BP(3, j, bpInd))^2;
                colCnt = colCnt+1;
            end
            for bpInd1 = 1:N
                for bpInd2 = bpInd1+1:N
                    L(cnt, colCnt) = 2*(BP(3, i, bpInd1) - BP(3, j, bpInd1))*(BP(3, i, bpInd2) - BP(3, j, bpInd2));            
                    colCnt = colCnt+1;
                end
            end
            R(cnt) = norm(C(:, i) - C(:, j))^2;
            cnt = cnt+1;
        end
    end
end
