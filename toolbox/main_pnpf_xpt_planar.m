clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('DLT'));
addpath(genpath('C:\Users\Александр Вахитов\Documents\materials\pnp\GPnPf-Toolbox\upnpPlanar_codeRelease\upnpPlanar_codeRelease'))
addpath(('C:\Users\Александр Вахитов\Documents\materials\pnp'))

warning off;

% experimental parameters
nls= 1:1:5;
npt= 6;
num= 50;

pnpOpts.fMax = 1e4;
pnpOpts.fMin = 30;
pnpOpts.errThr = 7.5;
pnpOpts.isFastBA = 0;

% compared methods
A= zeros(size(nls));
B= zeros(num,size(nls));
name= {'HOMO', 'UPnPf', 'GPnPf',   'GPnPf+GN', 'RPnP', 'EPnPR'};
f= {   @HOMO, @upnp_planar_interface, @GPnP_f, @GPnP_f_GN, @RPnP_interface, @(x,y)pnpfmy(x,y,3,1,pnpOpts)};
marker= { 'x', 'o', '>', 's', '+','<'};
color= {'r','g', 'm','k','c','b'};
markerfacecolor=  {'r','g','m','n','n','r'};
linestyle= {'-','-','-','-','-','-'};

inds = [2 4 5 6];
% inds = [1:7]
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
        %f= 800;
        f= rand(1)*1800+200;
        
        % generate 3d coordinates in camera space
        XXw= [xrand(2,npt,[-2 2]); zeros(1,npt)];
		R= rodrigues(randn(3,1));
        t= [rand-0.5;rand-0.5;rand*8+4];
        Xc= R*XXw+repmat(t,1,npt);
        
        % projection
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;


        % pose estimation
        for k= 1:length(method_list)
             time1 = tic;
             if strcmp(method_list(k).name, 'RPnP')
                 [f1,R1,t1]= method_list(k).f(XXw,xxn,diag([f,f,1]));
             else
                try
                    [f1,R1,t1]= method_list(k).f(XXw,xxn);
                catch
                    fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                    index_fail = [index_fail, j];
                    break;
                end
             end
                %no solution
            while size(t1,2) < 1
                 [f1,R1,t1]= method_list(k).f(XXw,xxn);
            end
            
            s = toc(time1 );
            if (s > 0.2)
                s
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
            method_list(k).tm(j, i) = s;
            reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
            reproj = reproj./repmat(reproj(3,:),3,1);
            err = xxn-f1(index_best)*reproj(1:2,:);
            err = sqrt(sum(sum(err.*err))/npt);
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
        method_list(k).med_reproj(i)= median(method_list(k).reproj(:,i));
        method_list(k).avg_t(i)= sum(method_list(k).tm(:,i)) / size(method_list(k).tm, 1);        
    end
end

close all;
yrange= [0 10];

i= 0; w= 300; h= 300;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_r','Rotation',...
    'Gaussian Image Noise (pixels)','Rotation Error (degrees)',2);

yrange= [0 40];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_t','Translation',...
    'Gaussian Image Noise (pixels)','Translation Error (%)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_foc','Focal Length',...
    'Gaussian Image Noise (pixels)','Focal Length Error (%)',2);

yrange= [0 15];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'med_reproj','Reprojection',...
    'Gaussian Image Noise (pixels)','Reprojection Error (pixels)',2);


yrange= [0 0.05];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(nls,yrange,method_list,'avg_t','Time',...
    'Average Runtime (sec)','Average Runtime (sec)',2);

rmpath(genpath('rpnp1.0'));
rmpath(genpath('DLT'));
rmpath(genpath('C:\zheng\work0919\PnP+f\upnpPlanar_codeRelease'))