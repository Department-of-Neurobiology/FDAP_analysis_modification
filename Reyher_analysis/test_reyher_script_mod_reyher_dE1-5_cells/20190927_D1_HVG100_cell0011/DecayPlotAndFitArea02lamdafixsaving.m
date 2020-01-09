function decayplotandfitarea
% lambda fixed
clear all, close all
[fn,path]=uigetfile('*.mat','Select any mat-file in folder to be treated.');
cd(path)
if ~exist('All.mat'); msgbox('No all.mat - returned.'); return, end

L=load('All.mat');
tP1binx=L.tP1binx;
tI1binAvg=L.tI1binAvg;
NumberOfFrames=L.NumberOfFrames;

%%%%%%%%%%%%
% Version 02 - offset not fitted but set to zero 
% See lines  156,157 and 242,243
pausetime=1;
NPixels = 256;
xScale  = 100/NPixels;       % Mikrons / Pixel
yScale=1;                    % Seconds / Frame
w=5;                         % region of activation in a fixed sample in mikrons
checkT1=10;                  % check out times in Secs
checkT2=50;
StartF=1;
D=2;                          % Estimate of D for fit
lambda=0.001;                 % fixed
%%%

NofF=NumberOfFrames-StartF+1; % number of frames
t=[0:(NofF-1)]' * yScale;         % time (column)vector
fighndl=figure(2);
plot(tP1binx(StartF,:)*xScale,tI1binAvg(StartF,:))
xlabel('Position on main axis 1 in mikrons')
XLimits=get(gca,'Xlim');set(gca,'XMinorTick','on');

YLimits=showRotCell(xScale,XLimits);
TEXT='Select tip edge, then body edge (2 clicks)';
texthndl=text('Position',[.2,.2],'String',TEXT,...
        'FontSize',10,'Units','normalized');
drawnow
[xc,yc]=ginput(2);
xTip=xc(1)
xBody=xc(2)
abs(xTip-xBody)



figure(fighndl);
pict=StartF;
CBD=struct('pict',pict,'NofF',NofF,'tP1binx',tP1binx,'tI1binAvg',tI1binAvg,'xScale',xScale); % Callbacks der Knöpfe benötigen diese Daten
guidata(fighndl,CBD); % CallBackData in Datenbereich der figure
evalin('base','CBD=guidata(gcf)'); %CBD im Base-Bereich für Callbacks erzeugen

bup = uicontrol('Parent',fighndl, ...
    'Units','Normalized', ...
    'Callback',['CBD.pict=CBD.pict+1;CBD.pict=min(CBD.pict,CBD.NofF);',...
    'plot(CBD.tP1binx(CBD.pict,:)*CBD.xScale,CBD.tI1binAvg(CBD.pict,:));title(num2str(CBD.pict));'], ...
    'Position',[0.51 0 0.1 0.05], ...
    'String','>', ...
    'Tag','Pushbutton1');

bdwn= uicontrol('Parent',fighndl, ...
    'Units','Normalized', ...
    'Callback',['CBD.pict=CBD.pict-1;',...
    'CBD.pict=max(CBD.pict,1);plot(CBD.tP1binx(CBD.pict,:)*CBD.xScale,CBD.tI1binAvg(CBD.pict,:));title(num2str(CBD.pict));'], ...
    'Position',[0.39 0 0.1 0.05], ...
    'String','<', ...
    'Tag','Pushbutton2');

bst = uicontrol('Style','PushButton',...
    'Parent',fighndl, ...
    'String','Cont.', ...
    'Callback',['[hbutton,hfigure]=gcbo; uiresume(hfigure);'], ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','StopButton');

uiwait(fighndl);

TEXT='Select center of excitation peak by 1 left button mouse click.';
[ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
from=ileft;to=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
xchannels=to-from+1;    % number of x-bins in w
hold on
plot(a,b,'r')
figure(1)
hold on
plot([a(1) a(1)],YLimits,'m',[a(xchannels) a(xchannels)],YLimits,'m') %plot vertical lines at "center +/- w/2"
% figure(2)
% TEXT='Select left side region by 1 left button mouse click.';
% [ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
% fromL=ileft;toL=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
% xchannelsL=toL-fromL+1;    % number of x-bins in w
% hold on
% plot(a,b,'r')
% figure(1)
% plot([a(1) a(1)],YLimits,'m',[a(xchannelsL) a(xchannelsL)],YLimits,'m')
% figure(2)
% TEXT='Select right side region by 1 left button mouse click.';
% [ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
% fromR=ileft;toR=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
% xchannelsR=toR-fromR+1;    % number of x-bins in w
% plot(a,b,'r')
% figure(1)
% plot([a(1) a(1)],YLimits,'m',[a(xchannelsR) a(xchannelsR)],YLimits,'m')
figure(2)
hold off
pause(pausetime)

%%%%%%%%%%%%% extract sigma0 from start frame ----------------
figure(3) 
Xs=tP1binx(StartF,:)*xScale; Ys=tI1binAvg(StartF,:);
% Normalize
dx=abs(Xs(1)-Xs(2));
AreaSF=sum(Ys)*dx;
Ys=Ys/AreaSF;
%
plot(Xs,Ys,'.') %Plot full intensity distribution for the start frame in blue and w-region in red
xlabel('Position on main axis 1 in mikrons')
hold on
plot(Xs(from:to),Ys(from:to),'r.')
    % plot(tP1binx(StartF,fromL:toL)*xScale,tI1binAvg(StartF,fromL:toL),'r')
    % plot(tP1binx(StartF,fromR:toR)*xScale,tI1binAvg(StartF,fromR:toR),'r')
% Fitting start frame Ys(Xs) 
% parameter estimates
[height,iheight]=max(Ys);
position=Xs(iheight);
wValue=0.6*height; iwValue=find(Ys(iheight:length(Ys))< wValue,1)+iheight;
width=Xs(iwValue)- position;
offset=Ys(1);
Pstart=[width, position, offset]
% show start estimates
hL=plot(Xs,Gsimu(Pstart,Xs),'g');
hT=title('Start curve. Paused for 5s.');
pause(pausetime)
delete(hT); delete(hL);
% do the fit
Pout=fminsearch(@(P) chi(P,Xs,Ys),Pstart, optimset('Display','iter'));
Pout
sigma0=Pout(1);
offsetSF=Pout(3);
ActPos=Pout(2); % Distance from tip in µm
% show the result
hL=plot(Xs,Gsimu(Pout,Xs),'g');
hT=title(['Fitted curve, Sigma0 = ', num2str(sigma0),' . Paused for 5s.']);
pause(pausetime)
hold off %figure(3)
saveas(3,'StartFramelambdafix.jpg')
%%%%%%%%%%%%% Making decay curve and fitting for D ----------------
Ic=tI1binAvg(StartF:NofF,from:to)/AreaSF; % Projected intensity in central region w, normalized    
CentralArea=sum(Ic,2)*dx;
% Averaged intensity in region w 
Iavg = CentralArea/abs(Xs(from)-Xs(to)); % 
figure(4)
plot(t,Iavg,'.') % plot Iavg(t)
hold on
title(['Center Int., avg over w ,' path]) % 
xlabel('Time in secs')
iT1=checkT1 /yScale +1;
iT2=checkT2 /yScale +1;
text(0.5,0.8,['Center Int. at ' num2str(checkT1) ' sec: ' num2str(Iavg(iT1))],'units','normalized')
text(0.5,0.7,['Center Int. at ' num2str(checkT2) ' sec: ' num2str(Iavg(iT2))],'units','normalized')
text(0.5,0.6,[cd],'units','normalized')
% fitting Iavg(t)
% parmeter estimates, partly from start frame (see above)
Pstart=[sigma0 D offsetSF ];   %  
HP=[ActPos, xTip, xBody, lambda]; %constant parameters
% show start estimates
hL=plot(t,Dsimu(Pstart,t,HP),'g');
hT=title('Start curve. Paused for 5s.');
pause(pausetime)
delete(hT); delete(hL);
% do the fit
Pout=fminsearch(@(P) Dchi(P,t,Iavg,HP),Pstart, optimset('Display','iter'));
D=Pout(2);
Pout
% show the result
hL = plot(t,Dsimu(Pout,t,HP),'g');
hL = plot(t,DsimuX(Pout,t,HP),'r--');
hT = title(['I_C, D = ', num2str(D),' .']);
axis([0 max(t)*1.05 0 max(Iavg)*1.05]);
hold off           %figure(4)
saveas(4,'NewFitlambdafix.fig')
saveas(4,'NewFitlambdafix.jpg')

%--------------------------------------------------------------------------
% Subfunctions

function [ileft,iright,x,y]=getlimits(fighndl,w,TEXT)
% function [ileft,iright,x,y]=getlimits(fighndl)
% Lets you select the central position of the excitation peak on a line plot in figure with
% handle fighndl. From this position, 2 positions  w/2 to the left and w/2 to the rigth are calculated.
% Returns indices ileft and iright of the data at these
% positions, as well as the x- and y-values within that range.
set(gcf,'pointer','fullcrosshair');
texthndl=text('Position',[.2,.2],'String',TEXT,...
        'FontSize',10,'Units','normalized');
drawnow;
[xc,yc]=ginput(1);
linehndl=findobj(fighndl,'type','line','-and','color','blue');
x=get(linehndl,'xdata');
y=get(linehndl,'ydata');
idx=find(x > xc-w/2);
ileft=min(idx);
idx=find(x < xc+w/2);
iright=max(idx);
x=x(ileft:iright);
y=y(ileft:iright);
delete(texthndl);

% Subfunction
function YLimits=showRotCell(xScale,XLimits)
CLineSet=[30 50 100 300 500 1000 2000];
rotcellh=openfig('RotCell.fig');
HG=findobj(rotcellh,'type','hggroup')
X=get(HG,'xdata');
Y=get(HG,'ydata');
Z=get(HG,'zdata');
Xs=X*xScale;Ys=Y*xScale;
contour(Xs,Ys,Z, CLineSet-2);
raxh=gca;
set(raxh,'YlimMode','auto');
YL=get(raxh,'Ylim'); YW = YL(2) - YL(1);
YLimits=[YL(1)-0.2*YW YL(2)+0.2*YW];
set(raxh,'Ylim',YLimits);
set(raxh,'Xlim',XLimits)
drawnow


% subfunctions
function out=chi(P,x,y)
S=Gsimu(P,x);
out=sum(sum((y-S).^2));

function S=Gsimu(P,x)
S = 1/(sqrt(2*pi)*P(1))* exp(-(x-P(2)).^2 / (2*P(1)^2)) + P(3);

% subfunctions
function out=Dchi(P,x,y,HP)
S=Dsimu(P,x,HP);
out=sum(sum((y-S).^2));

function S=Dsimu(P,x,HP)
sigma0=P(1);lambda=HP(4);
ActPos=HP(1); xTip=HP(2); xBody=HP(3);
sigma=sqrt(2*P(2)*x + sigma0^2);
S =exp(-lambda*x) / sqrt(2*pi)./ sigma .*(1+exp(-4*(ActPos-xTip)^2 ./ (2*sigma.^2))...   
                                   - exp(-4*(ActPos-xBody)^2 ./ (2*sigma.^2))) +P(3) ;

function S=DsimuX(P,x,HP)
sigma0=P(1);lambda=HP(4);
ActPos=HP(1); xTip=HP(2); xBody=HP(3);
sigma=sqrt(2*P(2)*x + sigma0^2);
S = exp(-lambda*x)/sqrt(2*pi) ./ sigma  + P(3);          % 


