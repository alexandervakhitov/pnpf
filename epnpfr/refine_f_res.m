function res = refine_f_res(x, L, sol1)        
    fGuess = x;
    sol = [sol1; fGuess^2*sol1];
    res = L*sol;
end