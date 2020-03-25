function [err, projerr] = calcProjErr(R1, t1, P, f1, k1, Uc)

    P1 = [R1 t1];
    if (size(P1, 1) > 0)
        err = 0;
        projerr = zeros(size(Uc, 2), 1);
        for i = 1:size(Uc, 2)        
            proj = P1 * [P(:, i); 1];
            proj = proj / proj(3);
            proj = distort1(proj, k1, f1);
%             proj = distortFitz(proj(1:2), k1, f1);
            err = err + norm(proj(1:2) - Uc(:, i));
            projerr(i) = norm(proj(1:2) - Uc(:, i));
        end
        err = err / size(Uc, 2);      
    else
        err = -1;
        projerr = [];
    end
end

