clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('DLT'));
addpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))

warning off;

% experimental parameters
%nls= 0:0.5:5;
nls= 1:5:20;
npt= 8;
num= 5;

% compared methods
A= zeros(size(nls));
B= zeros(num,size(nls));

addpath('C:\Users\Александр Вахитов\Documents\materials\pnp\');
addpath('C:\Users\Александр Вахитов\Documents\materials\pnp\GPnPf-Toolbox\UPnP');
addpath('C:\Users\Александр Вахитов\Documents\materials\pnp\GPnPf-Toolbox\UPnP\pose_functions');
name= {'DLT', 'UPnPf', 'UPnPf+GN', 'GPnPf',   'GPnPf+GN', 'RPnP', 'My'};
f= {   @(a,b) refined_dist_method(a, b, @DLT) , @(a,b) refined_dist_method(a, b, @upnp_interface), ...
    @(a,b)refined_dist_method(a, b, @upnp_GN_interface), @(a,b)refined_dist_method(a, b, @GPnP_f), @(a,b)refined_dist_method(a, b, @GPnP_f_GN), ...
    @(a,b)refined_dist_method(a, b, @GPnP_f_GN), @(x,y)pnpfmyDist1(x,y,4)};
marker= { 'x', 'o', 'd', '>', 's', '+', '+'};
color= {'r','g','b','m','k','c', 'r'};
markerfacecolor=  {'r','g','n','m','n','n', 'n'};
linestyle= {'-','-','-','-','-','-','-'};

inds = [3 5 7];
f = f(inds);
marker = marker(inds);
color = color(inds);
markerfacecolor = markerfacecolor(inds);
linestyle = linestyle(inds);
name = name(inds);
method_list= struct('name', name, 'f', f, 'mean_r', A, 'mean_t', A, 'mean_foc', A, 'mean_reproj', A, ...
    'med_r', A, 'med_t', A, 'med_foc', A, 'med_reproj', A, 'r', B, 't', B,...
    'foc', B, 'reproj', B, 'marker', marker, 'color', color, ...
    'markerfacecolor', markerfacecolor, 'linestyle', linestyle);

% experiments
for i= 1:length(nls)
    
    nl= nls(i);
    fprintf('nl = %.1f: ',nl);
    
    index_fail = [];
    for j= 1:num
        
        % camera's parameters
        width= 640;
        height= 480;
        f= rand(1)*1800+200; %random focal length in [200,2000]
        
        % generate 3d coordinates in camera space
        Xc= [xrand(1,npt,[-2 2]); xrand(1,npt,[-2 2]); xrand(1,npt,[4 8])];
        t= mean(Xc,2);
        R= rodrigues(randn(3,1));
        XXw= inv(R)*(Xc-repmat(t,1,npt));
        
        % projection
        kt = 3e-1;
        xxn = projectionDistD1(Xc, f, kt, 0, 0, nl);        
        %xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        %xxn= xx+randn(2,npt)*nl;

        % pose estimation
        for k= 1:length(method_list)
             if strcmp(method_list(k).name, 'RPnP')
                 [f1,R1,t1,k1]= method_list(k).f(XXw,xxn,diag([f,f,1]));
             else
                try
                    [f1,R1,t1,k1]= method_list(k).f(XXw,xxn);
                catch
                    fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                    index_fail = [index_fail, j];
                    break;
                end
             end
            
            %no solution
            while size(t1,2) < 1
                 [f1,R1,t1,k1]= method_list(k).f(XXw,xxn);
            end
            
            %choose the solution with smallest error 
            index_best = 1;
            error = inf;
            for jjj = 1:size(R1,3)
                tempy = cal_pose_err([R1(:,:,jjj) t1(:,jjj)],[R t]);
                if sum(tempy) < error
                    y = tempy;
                    error = sum(tempy);
                    index_best = jjj;
                end
            end

            method_list(k).r(j,i)= y(1);
            method_list(k).t(j,i)= y(2);
            method_list(k).foc(j,i)= abs(f1(index_best)-f)/f*100;
            method_list(k).k1(j,i)= abs(k1-kt)/abs(kt)*100;
            reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
            reproj = reproj./repmat(reproj(3,:),3,1);
            [err, projerr] = calcProjErrD1(R1(:,:,index_best), t1(:,index_best), XXw, f1(index_best), k1, xxn);    
            %err = xxn-f1(index_best)*reproj(1:2,:);
            %err = sqrt(sum(sum(err.*err))/npt);
            method_list(k).reproj(j,i) = err;
        end

        showpercent(j,num);
        end
        fprintf('\n');


        % save result
        for k= 1:length(method_list)
            method_list(k).med_r(i)= median(method_list(k).r(:,i));
            method_list(k).med_t(i)= median(method_list(k).t(:,i));
            method_list(k).med_foc(i)= median(method_list(k).foc(:,i));
            method_list(k).med_k1(i)= median(method_list(k).k1(:,i));
            method_list(k).med_reproj(i)= median(method_list(k).reproj(:,i));
            method_list(k).med_reproj(i)= median(method_list(k).reproj(:,i));
        end
end

close all;
yrange= [0 5];

i= 0; w= 300; h= 300;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_r','Rotation',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)',2);

yrange= [0 10];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_t','Translation',...
    'Gaussian Image Noise (pixels)','Translation Error (%)',2);
yrangef = [0 50];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrangef,method_list,'med_foc','Focal Length',...
    'Gaussian Image Noise (pixels)','Focal Length Error (%)',2);

yrangef = [0 150];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrangef,method_list,'med_k1','Distortion Coef',...
    'Gaussian Image Noise (pixels)','Dist Coef Error (%)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_reproj','Reprojection',...
    'Gaussian Image Noise (pixels)','Reprojection Error (pixels)',2);

rmpath(genpath('rpnp1.0'));
rmpath(genpath('DLT'));
rmpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))