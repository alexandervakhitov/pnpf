function [Rest, test, fest, dest, avgerr] = refine_pos_dist(Rest, test, fest, dest, projs, pts3d, u, v, dist, isFast, fMin, fMax, isFastBA, haveMex)
    noConstr = 0;
    %good and slow
    if (isFast == 0)
        options=[1E-01, 1E-15, 1E-15, 1E-15, 1e-6];
    else
    %fast
        if (isFastBA == 1)
            options=[1E-03, 1E-5, 1E-5, 1E-5, 1e-6];
        else
            options=[1E-01, 1E-10, 1E-10, 1E-10, 1e-6];
        end        
        noConstr = 1;
    end
    itmax = 200;
    imvect = projs(:);
    
    [U,S,V] = svd(Rest);    
    Rest = U*V';
    
    x = zeros(6, 1);    
    
    if (det(Rest) < 0)
        Rest = -1*Rest;
        test = -1 * test;
    end
    if (dist==1)
        x(1:6) = camops.encodeCamera(Rest, test);
    %     [ret, popt, info, covar]=levmar('refine_pos_res_d', x, imvect, itmax, options, pts3d, u, v, fest, dest);    
    %     x = zeros(7, 1);
    %     x(1:6) = popt(:);
        x(7) = fest;
        x(8) = dest;
        %[ret, popt, info, covar]=levmar('refine_pos_foc_d_res', x, imvect, itmax, options, pts3d, u, v);    
        lb = [-1e10*ones(6,1); 50; -0.5];
        ub = [1e10*ones(6,1); maxF; 0.5];
        [ret, popt, info, covar]=levmar('refine_pos_foc_d1_res', x, imvect, itmax, options, 'bc', lb, ub, pts3d, u, v);    
        [Rest, test] = camops.decodeEucCamera(popt(1:6));
        fest = popt(7);
        dest = popt(8);
        avgerr = sqrt(info(2) / size(pts3d, 2));
    else
        if (dist == 0)
            x(1:6) = encodeCamera(Rest, test);
            x(7) = fest;
            
            if (haveMex)
                resFunName = 'refine_pos_foc_res2_mex';
            else
                resFunName = 'refine_pos_foc_res2';
            end
                
            if (noConstr == 1)                
                [ret, popt, info, covar1]=levmar(resFunName, x, imvect, itmax, options, pts3d);    
            else
                lb = -1e10*ones(1,7);
                lb(7) = fMin;
                ub = 1e10*ones(1,7);
                ub(7) = fMax;
                [ret, popt, info, covar1]=levmar(resFunName, real(x), imvect, itmax, options, 'bc', lb, ub, pts3d);                    
            end
            [Rest, test] = decodeEucCamera(popt(1:6));
            fest = popt(7);
            avgerr = sqrt(info(2) / size(pts3d, 2));        
        else
            if dist == 2
                x(1:6) = camops.encodeCamera(Rest, test);
            %     [ret, popt, info, covar]=levmar('refine_pos_res_d', x, imvect, itmax, options, pts3d, u, v, fest, dest);    
            %     x = zeros(7, 1);
            %     x(1:6) = popt(:);
                x(7) = fest;            
                [ret, popt, info, covar]=levmar('refine_pos_foc_d_fix_res', x, imvect, itmax, options, pts3d, dest, fest, u, v);    
                [Rest, test] = camops.decodeEucCamera(popt(1:6));
                %fest = popt(7);            
                avgerr = sqrt(info(2) / size(pts3d, 2));
            else 
                x(1:6) = camops.encodeCamera(Rest, test);
                [ret, popt, info, covar]=levmar('refine_pos_res_d', x, imvect, itmax, options, pts3d, u, v, fest, dest);    
                [Rest, test] = camops.decodeEucCamera(popt(1:6));
                avgerr = sqrt(info(2) / size(pts3d, 2));                
            end
        end
    end
end