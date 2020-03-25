function v = rotMatToVect(m)
    a = rodrigues(m);
    v = a(1:3,1);
end