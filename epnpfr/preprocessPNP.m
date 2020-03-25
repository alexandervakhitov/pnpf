function [C V D] = preprocessPNP(P, U0, tarPtNum)
    %P = [x1...xn; y1...yn; z1...zn]; - points
    %U = [u1...un; v1...vn]; - projections
    % (u,v) - principle point
    % N: X = B1*V1+...+BN*VN;
    %type=1:normal, type=2:planar, type=3:distortion
    
    doTimeMes = 1;
    
    U = U0;
    
    % n - number of points
    [m n] = size(P);        
    % create C - new coordinates
    if (doTimeMes)
%          sv = tic;
    end
    [A C P_ed] = create_C_A(P, tarPtNum);
    if (doTimeMes)
%          s1prep = toc(sv)
    end
    %find aij
    % P=CA

    % Mx = 0    
%      M = zeros(2*n, 3*tarPtNum);
    up1 = -diag(U(1,:))*A';
    up2 = -diag(U(2,:))*A';
    
%     M = [A' zeros(size(A')) up1;
%         zeros(size(A')) A' up2];
    
%     for i = 1:1:n
%         for j = 1:1:tarPtNum
%             M((i-1)*2 + 1,((j-1)*3 + 1):3*j) = [A(j,i) 0 -U(1,i)*A(j,i)];
%             M(2*i,((j-1)*3 + 1):3*j) = [0 A(j,i) -U(2,i)*A(j,i)];
%         end
%     end
       
    if (doTimeMes)
%          s2prep = toc(sv)
    end

%     M1 = M'*M;
    A2 = A*A';
    AU1 = A*up1;
    AU2 = A*up2;
    U12 = up1'*up1 + up2'*up2;
    M1 = [A2 zeros(tarPtNum, tarPtNum) AU1;
          zeros(tarPtNum, tarPtNum) A2 AU2;
          AU1' AU2' U12];
    if (doTimeMes)
%          s3prep = toc(sv)
    end
     [A1, D0, V0] = svd(M1);
%     [A11, D, V0] = svd(M1'*M1);
    if (tarPtNum == 4)
        inds = [1 5 9 2 6 10 3 7 11 4 8 12];
    else
        inds = [1 4 7 2 5 8 3 6 9];
    end
    V = V0(inds, :);
    D = D0(inds, inds);
    
    
    if (doTimeMes)
%          s4prep = toc(sv)
    end
        
end