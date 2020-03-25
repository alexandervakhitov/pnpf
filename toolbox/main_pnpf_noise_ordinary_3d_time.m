clear; clc;
addpath(genpath('rpnp1.0'));
addpath(genpath('DLT'));
addpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))
addpath(genpath('C:\Users\Александр Вахитов\Documents\materials\pnp\GPnPf-Toolbox\UPnP'));
warning off;
% experimental parameters
nl= 2;
npts= 10:100:1000;
num= 100;



pnpOpts.fMax = 2e4;
pnpOpts.fMin = 50;
pnpOpts.errThr = 7.5;
pnpOpts.isFastBA = 1;

pnpOptsSlow = pnpOpts;
pnpOptsSlow.isFastBA = 0;

% compared methods
A= zeros(size(npts));
B= zeros(num,size(npts));
name= {'DLT', 'UPnPf', 'UPnPf+GN', 'GPnPf',   'GPnPf+GN', 'RPnP', 'EPnPfRFast', 'EPnPfR'};
f= {   @DLT, @upnp_interface,  @upnp_GN_interface, @GPnP_f, @GPnP_f_GN, @RPnP_interface, @(x,y)pnpfmy(x,y,4,1,pnpOpts), @(x,y)pnpfmy(x,y,4,1,pnpOptsSlow)};
marker= { 'x', 'o', 'd', '>', 's', '+', '<', 'x'};
color= {'r','g','b','m','k','c','b','r'};
markerfacecolor=  {'r','g','n','m','n','n','r','g'};
linestyle= {'-','-','-','-','-','-','-','-'};

%inds = [3 5 6 7];
inds = [3 5 7 8];
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
for i= 1:length(npts)
    
    npt= npts(i);
    fprintf('npt = %d: ',npt);
    
  
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
        xx= [Xc(1,:)./Xc(3,:); Xc(2,:)./Xc(3,:)]*f;
        xxn= xx+randn(2,npt)*nl;

        if npt >= 6
            % pose estimation
            for k= 1:length(method_list)
                time1 = tic;
                try
                    [f1,R1,t1]= method_list(k).f(XXw,xxn);
                catch
                    fprintf(['   The solver - ',method_list(k).name,' - encounters internal errors! \n']);
                    index_fail = [index_fail, j];
                    break;
                end

                %no solution
                while size(t1,2) < 1
                    [f1,R1,t1]= method_list(k).f(XXw,xxn);                    
                end                                
                timeMes = toc(time1);
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
                method_list(k).tm(j, i) = timeMes;
                reproj = R1(:,:,index_best)*XXw+t1(:,index_best)*ones(1,npt);
                reproj = reproj./repmat(reproj(3,:),3,1);
                err = xxn-f1(index_best)*reproj(1:2,:);
                err = sqrt(sum(sum(err.*err))/npt);
                method_list(k).reproj(j,i) = err;
            end
        end
    
        showpercent(j,num);
    end
    fprintf('\n');
    
    % save result
    for k= 1:length(method_list)
       
        method_list(k).mean_r(i)= mean(method_list(k).r(:,i));
        method_list(k).mean_t(i)= mean(method_list(k).t(:,i));
        method_list(k).mean_foc(i)= mean(method_list(k).foc(:,i));
        method_list(k).mean_reproj(i)= mean(method_list(k).reproj(:,i));
        
        method_list(k).med_r(i)= median(method_list(k).r(:,i));
        method_list(k).med_t(i)= median(method_list(k).t(:,i));
        method_list(k).med_foc(i)= median(method_list(k).foc(:,i));
        method_list(k).med_reproj(i)= median(method_list(k).reproj(:,i));        
        method_list(k).avg_t(i)= sum(method_list(k).tm(:,i)) / size(method_list(k).tm, 1);
    end
end

%draw the boxplot of rotation error
close all;
yrange = [0 10];
i= 0; w= 300; h= 300;

i = 0;

yrange= [0 0.05];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'avg_t','Time',...
    'Average Runtime (sec)','Average Runtime (sec)',2);
yrange= [0 5];
figure('color','w','position',[w*(i+1),100,w,h]);
hold('all');
for i = 1:4
    p(i)=plot(npts, method_list(i).mean_reproj,'marker',method_list(i).marker,...
        'color',method_list(i).color,...
        'markerfacecolor',method_list(i).markerfacecolor,...
        'displayname',method_list(i).name, ...
        'LineWidth',2,'MarkerSize',8,'LineStyle',method_list(i).linestyle)
end
title('Mean Reprojection Error','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Reprojection Error (pixels)','FontSize',11);
legend(p);


rmpath(genpath('rpnp1.0'));
rmpath(genpath('DLT'));
rmpath(genpath('C:\zheng\work0919\PnP+f\UPnP'))

% close all;
% yrange= [0 15];
% 
% i= 0; w= 300; h= 300;
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_r','Mean Rotation Error',...
%     'Number of Points','Rotation Error (degrees)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_r','Median Rotation Error',...
%     'Number of Points','Rotation Error (degrees)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_t','Mean Translation Error',...
%     'Number of Points','Translation Error (%)');
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_t','Median Translation Error',...
%     'Number of Points','Translation Error (%)');
% 
% i= 0; w= 300; h= 300;
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_foc','Mean Focal Length Error',...
%     'Number of Points','Focal Length Error (%)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_foc','Median Focal Length Error',...
%     'Number of Points','Focal Length Error (%)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'mean_reproj','Mean Reproj Error',...
%     'Number of Points','Reproj Error (pixels)');
% 
% figure('color','w','position',[w*i,400,w,h]);i=i+1;
% xdrawgraph(npts,yrange,method_list,'med_reproj','Median Reproj Error',...
%     'Number of Points','Reproj Error (pixels)');
% 
