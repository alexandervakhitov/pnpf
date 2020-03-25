clear; clc;
addpath(genpath('toolbox\rpnp1.0'));
addpath(genpath('toolbox\DLT'));
addpath(genpath('toolbox\UPnP'))
addpath(genpath('toolbox'));
addpath(genpath('epnpfr'));
addpath(genpath('levmar\matlab'))
warning off;
% experimental parameters
nls = 4:1:4;
npts= 6;
num= 500;

% compared methods
A= zeros(size(npts));
B= zeros(num,size(npts));
name= {'DLT', 'UPnPf', 'UPnPf+GN', 'GPnPf',   'GPnPf+GN', 'RPnP', 'EPnPfR'};
f= {   @DLT, @upnp_interface,  @upnp_GN_interface, @GPnP_f, @GPnP_f_GN, @RPnP_interface, @epnpfr};
marker= { 'x', 'o', 'd', '>', 's', '+', '<'};
color= {'r','g','b','m','k','c','b'};
markerfacecolor=  {'r','g','n','m','n','n','r'};
linestyle= {'-','-','-','-','-','-','-'};

% inds = [5 7];
inds = [3 5 6 7];
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

mesArr  = [];
lstab_graph = [];
mctab_graph = [];
% experiments
for i= 1:length(nls)
    
    resArrNoReg = [];
    resArrReg = [];
    npt= npts;
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
        xxn= xx+randn(2,npt)*nls(i);
        N=2;
        
        [resObjNoReg resObjReg]= epnpfr_regtest(XXw,xxn, R, t, f, N);                    
        resArrReg = [resArrReg; resObjReg];
        resArrNoReg = [resArrNoReg; resObjNoReg];
    end
    lstab = [];
    ranktab = [];
    lstabr = [];
    ranktabr = [];
    mctab = [];
    for p = 1:length(resArrReg)
%         lstab = [lstab resArrNoReg(p).redLsArr];
        ranktab = [ranktab; resArrNoReg(p).rankArr];
        size(resArrReg(p).redLsArr)
        lstab = [lstab resArrReg(p).redLsArr];
        ranktabr = [ranktabr; resArrReg(p).rankArr];
        mctab = [mctab resArrReg(p).mcArr];
%         [mv, mind] = min(b);
%         if (mind > 1)
%             mind;
%         end
%         minmc = [minmc; min(resArrReg(p).mcArr)];
    end
    lstab_fin = zeros(size(lstab, 1), 1);
    mctab_fin = zeros(size(lstab, 1), 1);
    for p = 1:size(lstab, 1)
        lstab_fin(p) = sum(lstab(p, :))/size(lstab, 2);
        mctab_fin(p) = sum(mctab(p, :))/size(mctab, 2);
    end
    
    
    
    mes = [];
    for it = 1:size(lstabr, 1)
        mes = [mes; sum(lstabr(it, :)) / size(lstabr, 2)];
    end
    mesArr = [mesArr mes];
    
    lstab_graph = [lstab_graph lstab_fin];
    mctab_graph = [mctab_graph mctab_fin];
end

figure(1);
colors = 'rgbcm';
hold on;
for p = 1:length(nls)
    plot(1:length(nls), lstab_graph(p, :), colors(p));
end
hold off;
close all;

figure(1);
colors = 'rgbcm';
hold on;
for p = 1:length(nls)
    plot(1:length(nls), mctab_graph(p, :), colors(p));
end
hold off;
yrange= [0 5];

i= 0; w= 300; h= 300;

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_r','Rotation',...
    'Number of Points','Rotation Error (degrees)',2);

yrange= [0 10];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_t','Translation',...
    'Number of Points','Translation Error (%)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_foc','Focal Length',...
    'Number of Points','Focal Length Error (%)',2);

figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'med_reproj','Reprojection',...
    'Number of Points','Reprojection Error (pixels)',2);
yrange= [0 0.05];
figure('color','w','position',[w*i,100,w,h]);i=i+1;
xdrawgraph(npts,yrange,method_list,'avg_t','Time',...
    'Average Runtime (sec)','Average Runtime (sec)',2);

%draw the boxplot of rotation error
close all;
yrange = [0 10];
i= 0; w= 300; h= 300;
for i = 1:length(inds)
    methodInd = inds(i);
    figure('color','w','position',[w*i,100,w,h]);
    boxplot(method_list(i).reproj,npts);
    ylim(yrange); set(gca,'xtick',npts);
    title(method_list(i).name,'FontSize',12,'FontName','Arial');
    xlabel('Number of Points','FontSize',11);
    ylabel('Reprojection Error (pixels)','FontSize',11);
end

figure('color','w','position',[w*(i+1),100,w,h]);
hold('all');
for i = 1:length(method_list)
    plot(6:15, method_list(i).mean_reproj(1:end),'marker',method_list(i).marker,...
        'color',method_list(i).color,...
        'markerfacecolor',method_list(i).markerfacecolor,...
        'displayname',method_list(i).name, ...
        'LineWidth',2,'MarkerSize',8,'LineStyle',method_list(i).linestyle)
end
title('Mean Reprojection Error','FontSize',12,'FontName','Arial');
xlabel('Number of Points','FontSize',11);
ylabel('Reprojection Error (degrees)','FontSize',11);

% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% boxplot(method_list(3).r,npts);
% ylim(yrange); set(gca,'xtick',npts);
% title('UPnPf+GN','FontSize',12,'FontName','Arial');
% xlabel('Number of Points','FontSize',11);
% ylabel('Rotation Error (degrees)','FontSize',11);
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% boxplot(method_list(4).r,npts);
% ylim(yrange); set(gca,'xtick',npts);
% title('GPnPf','FontSize',12,'FontName','Arial');
% xlabel('Number of Points','FontSize',11);
% ylabel('Rotation Error (degrees)','FontSize',11);
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% boxplot(method_list(5).r,npts);
% ylim(yrange); set(gca,'xtick',npts);
% title('GPnPf+GN','FontSize',12,'FontName','Arial');
% xlabel('Number of Points','FontSize',11);
% ylabel('Rotation Error (degrees)','FontSize',11);
% 
% figure('color','w','position',[w*i,100,w,h]);i=i+1;
% boxplot(method_list(7).r,npts);
% ylim(yrange); set(gca,'xtick',npts);
% title('EPNPfReg','FontSize',12,'FontName','Arial');
% xlabel('Number of Points','FontSize',11);
% ylabel('Rotation Error (degrees)','FontSize',11);
% 
% 
vars = zeros(length(npts), 2);
for i = 1:length(npts)
    for k = 1:2
        vars(i,k) = var(method_list(k).reproj(:,i));
    end
end

figure('color','w','position',[w,100,w,h]);
hold('all');
markerlist = {'o-', '*-'};
for i = 1:2
    p(i) = plot(npts, vars(:,i), markerlist {i}, 'displayname',method_list(i).name);
end
title('Variance of Errors','FontSize',12,'FontName','Arial');
xlabel('Gaussian Image Noise (pixels)','FontSize',11);
ylabel('Variance (pixels)','FontSize',11);
legend(p,1);


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
