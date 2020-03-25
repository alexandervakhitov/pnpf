function res = refine_pos_foc_res2(x, pts3d)
    ptnum = size(pts3d, 2);
    res = zeros(ptnum*2, 1);    
    f = x(7);
    [R t] = camops.decodeEucCamera(x(1:6));
    P = [R t];
    d = 0;
    for i = 1:ptnum
        pt3d = pts3d(:, i);
        proj = P * [pt3d; 1];
        proj = proj / proj(3);
        proj = distort1(proj, d, f); 
        res(2*i-1: 2*i) = proj(1:2);
    end
end