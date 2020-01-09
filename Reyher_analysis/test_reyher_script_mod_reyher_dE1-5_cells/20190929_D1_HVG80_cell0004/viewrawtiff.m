function viewrawtiff
clear all, close all
[fn,path]=uigetfile('*.tif','Select raw-tiff to view');
cd(path)

z=imread(fn); % z is type uint16 (0 to 65535)
% image(z)
% 
% imwrite(z(:,:,3),'test','jpg','Bitdepth',16)
% find color channel used
a=max(max(z));
%a=a(:); % make column
[dummy,channel]=max(a);
%
channel=3
z(:,:,2)=0*z(:,:,2);
% subtract offset
% set up filter function, Gaussian filter with radius rf
rf=1;
[xf,yf]= meshgrid(-rf:rf);
ff =-xf.*xf-yf.*yf;
ff =exp(ff/3/rf);
ff =ff/sum(sum(ff));
Z=double(z(:,:,3));
Z=filter2(ff,Z,'valid');
offset=min(min(Z));
Z=Z-offset;

% gamma correction and scale to 65500
gamma=3;
Z=Z.^(1/gamma);
Zmax=max(max(Z));
Z=Z/Zmax*65500;
surf(Z,'edgecolor','none')
zf(:,:,1)=uint16(zeros(size(Z)));
zf(:,:,3)=uint16(zeros(size(Z)));
zf(:,:,2)=uint16(Z);


% show image as RGB

image(zf), axis equal
title(fn,'Interpreter','none')
screensize=get(0,'ScreenSize');set(gcf,'Position',screensize);
saveas(gcf,[fn '.jpg'])

