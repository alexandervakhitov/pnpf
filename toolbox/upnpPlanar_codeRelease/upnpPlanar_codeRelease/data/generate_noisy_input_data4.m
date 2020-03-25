function [A,point,Rt]=generate_noisy_input_data4(n,std_noise,percentage_outliers,draw_plot)

%clear all; close all; load generate_noisy_input_data; std_noise=6;%0.005;
%clear all; close all; n=10; std_noise=5;
if nargin<4 draw_plot=''; end

%true 3d position of the points--------------
xini=1; xend=2;%bounds
yini=1; yend=2;
zini=4; zend=8;

%intrinsic camera parameters------------------
fx=800; %focal length in pixels
fy=800;
f=0.05; %50 mm of focal length

m=fx/f;
m=1; f=1; %<<<<<<----------------------------

u0=320; v0=240;
width=640; height=480;
A=[fx/m 0 u0/m 0; 0 fy/m v0/m 0; 0 0 1 0];
std_noise=std_noise*(1/m);

i=1;
num_outliers=ceil(n*percentage_outliers);  %num_outliers=4;
while i<=n
   point(i).Xcam(1,1)=random(xini,xend,1);
   point(i).Xcam(2,1)=random(yini,yend,1);
   point(i).Xcam(3,1)=random(zini,zend,1);

   %project points into the image plane
   point(i).Ximg_true=project_3d_2d(A,point(i).Xcam);
   point(i).Ximg_true(3)=f;

   %coordinates in pixels
   if i>n-num_outliers %generate outlier
       urand=random(1,width,1);
       vrand=random(1,height,1);
       point(i).Ximg_pix_true=[urand,vrand]';
   else
   point(i).Ximg_pix_true=point(i).Ximg_true(1:2)*m;
   end

  %check if the point is into the limits of the image
  uv=point(i).Ximg_pix_true;
  if uv(1)>0 & uv(1)<width & uv(2)>0 & uv(2)<height
      i=i+1;
  end




 %  fprintf('%.1f,%.1f,%d\n',uv(1),uv(2),i);
end

%add noise
for i=1:n
   noise=randn(1,2)*std_noise;
   point(i).Ximg(1,1)=point(i).Ximg_true(1)+noise(1);
   point(i).Ximg(2,1)=point(i).Ximg_true(2)+noise(2);
   point(i).Ximg(3,1)=f;

   %noisy coordinates in pixels
   point(i).Ximg_pix=point(i).Ximg(1:2)*m;





end

%the observed data will not be in the camera coordinate system. We put the center of the
%world system on the centroid of the data. We also assume that the data is rotated with
%respect to the camera coordenate system.
centroid=[0 0 0]';
for i=1:n
   centroid=centroid+point(i).Xcam;
end
centroid=centroid/n;



x_rotation=random(0,45,1)*pi/180;
y_rotation=random(0,45,1)*pi/180;
z_rotation=random(0,45,1)*pi/180;
tx=0; ty=0; tz=0;
Rt=return_Rt_matrix(x_rotation,y_rotation,z_rotation,tx,ty,tz);
for i=1:n
   point(i).Xworld=transform_3d(Rt,point(i).Xcam-centroid);

%     point(i).Xworld(1,1)=point(i).Xcam(1,1)+rand(1,1)*0.25;
%     point(i).Xworld(2,1)=point(i).Xcam(2,1)+rand(1,1)*0.25;
%     point(i).Xworld(3,1)=point(i).Xcam(3,1)+rand(1,1)*0.25;

end


%plot noisy points in the image plane
if ~strcmp(draw_plot,'donotplot')
   figure; hold on;
   for i=1:n
       plot(point(i).Ximg_pix_true(1),point(i).Ximg_pix_true(2),'.','color',[1 0 0]);
       %txt=sprintf('%d',i);
       %text(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),txt);
       plot(point(i).Ximg_pix(1),point(i).Ximg_pix(2),'o','color',[0 1 0],'markersize',5);
   end
   title('Noise in image plane','fontsize',20);
   grid on;
   xlim([0,width]);
   ylim([0,height]);
end






%draw 3d points
% figure; hold on;
% representation_offset=10;
% plot3(0,0,0,'.','color',[0 0 0],'markersize',20);
% for i=1:n
%    plot3(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),'.','color',[1 0 0],'markersize',12);
%    txt=sprintf('%d',i);
%    text(point(i).Xcam(1),point(i).Xcam(2),point(i).Xcam(3),txt);
%    plot3(point(i).Xworld(1)+representation_offset,point(i).Xworld(2),point(i).Xworld(3),'.','color',[0.8 0.8 0],'markersize',12);
%    text(point(i).Xworld(1)+representation_offset,point(i).Xworld(2),point(i).Xworld(3),txt);
%    plot3(point(i).Ximg(1),point(i).Ximg(2),point(i).Ximg(3),'.','color',[0 0 1],'markersize',12);
% end
% grid on;
%

function [x]=random(xini,xend,n)

a=xend-xini;
b=xini;


xu=rand(n,1);
x=xu*a+b;