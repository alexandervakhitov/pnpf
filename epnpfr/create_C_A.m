function [A C P_ed] = create_C_A(P, tarPtNum)
    P_x_sum = sum(P(1, :));
    P_y_sum = sum(P(2, :));
    P_z_sum = sum(P(3, :));
    n = size(P, 2);    
    C = zeros(4, 3);
    if (tarPtNum == 4)
        C1 = [P_x_sum/n; P_y_sum/n; P_z_sum/n];
        C2 = [P_x_sum/n + 1; P_y_sum/n; P_z_sum/n];
        C3 = [P_x_sum/n; P_y_sum/n + 1; P_z_sum/n];
        C4 = [P_x_sum/n; P_y_sum/n; P_z_sum/n + 1];
        C = [C1 C2 C3 C4];
        E_n = ones(1, n);
        P_ed = [P; E_n];
        E_4 = [1 1 1 1];
        C_ed = [C; E_4];
        A = C_ed \ P_ed; 
        
%         pNum = size(P_ed, 2);
%         for cInd = 1:4
%             for tInd = 1:4
%                 L(tInd, 1:pNum) = P_ed(tInd, :);
%                 b(tInd) = C_ed(tInd, cInd);            
%             end
%             for pInd = 1:pNum
%                 lrow = zeros(1, pNum+4);
%                 lrow(pInd) = 2;
%                 lrow(pNum+1: pNum+4) = - P_ed(:, pInd);
%                 L(pInd+4, 1:pNum+4) = lrow;
%                 b(pInd+4) = 0;
%             end
%             c1cg = L\b';
%             q1 = c1cg(1:pNum);
%             Q(cInd, :) = q1;
%         end
        %p = a' * C
        %R*p + t = R*(a'*C) + t
        %let c = P*q
        %R*c+t = R*P*q + t = [R t]*[P; 1]*q = P_m * q 
        %P,q known from beginning; R,t,f found => (uc, vc) = known
        %c = P*q
        %\sum q_i^2 -> min
        %\sum q_i^2 + g' * (c-P*q)
        %\nabla g_j : c_j = P_j*q  <=>  P_j*q = c_j
        %\nabla q_i : 2 q_i - g'*P_i = 0  
    else
        C1 = [P_x_sum/n; P_y_sum/n; P_z_sum/n];
        Pc = P - C1*ones(1, size(P, 2));
        [U,S,V] = svd(Pc');
        C2 = C1 + V(:, 1);
        C3 = C1 + V(:, 2);        
        C(1:3, 1:3) = [C1 C2 C3];
        C1 = C(1:3, 1:3);
        E_n = ones(1, n);
        P_ed = [P; E_n];        
        C_ed = [C1; ones(1, 3)];
        A = C_ed \ P_ed;         
    end    
end
