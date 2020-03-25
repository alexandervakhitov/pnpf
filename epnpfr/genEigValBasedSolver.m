function xkNs = genEigValBasedSolver(x0, kerVect, betaNum)    
    v = size(kerVect, 1);    
    A = kerVect(1:v/2, :);
    a = x0(1:v/2);
    B = kerVect(v/2+1:end, :);
    b = x0(v/2+1:end);
    At = pinv(A);
    at = pinv(A)*a;
    Q = B*At;
    q = -B*at+b;
    Qbig = [Q q; zeros(1, size(Q, 2)) 0];
    Ebig = [eye(size(Q)) zeros(size(Q, 1), 1);
            zeros(1, size(Q, 2)) 0];
    W = pinv(Qbig)*Ebig;
    [V, D] = eig(Qbig, Ebig);
    Abig = [A a; zeros(1, size(A, 2)) 0];
    Bbig = [B b; zeros(1, size(B, 2)) 0];
    Qbig2 = pinv(Bbig)*Abig;
    [V, D] = eig(Qbig2);
    d = diag(D);
    xkNs = [];
    for dInd = 1:length(d)
        ktr = real(d(dInd));
        M = A-ktr*B;
        m = ktr*b-a;
        for conInd = 1:betaNum
            PC = -A(conInd, :);
            pcvec = a(conInd);
            PC = [PC; -B(conInd, :)];
            pcvec = [pcvec; b(conInd)];
           [X,RESNORM,RESIDUAL,EXITFLAG] = lsqlin(M,m,PC, pcvec);
           if (EXITFLAG > 0)
              xkNs = [xkNs kerVect*X+x0]; 
           end
        end
    end    
    
end