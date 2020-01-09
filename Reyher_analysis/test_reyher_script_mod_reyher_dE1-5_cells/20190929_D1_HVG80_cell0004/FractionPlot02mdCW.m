function fractionPlot
clear all, close all
[fn,path]=uigetfile('*.mat','Select any mat-file in folder to be treated.');
cd(path)
if ~exist('All.mat'); msgbox('No all.mat - returned.'); return, end

load All.mat
NPixels=256;
xScale=100/NPixels;       % Mikrons / Pixel
yScale=1;                 % Seconds / Frame
w=5;                      % region of activation in a fixed sample in mikrons
checkT1=10;               % check out times in Secs
checkT2=50;

StartF=1;

NofF=NumberOfFrames-StartF+1; % number of frames
fighndl=figure(2);
plot(tP1binx(StartF,:)*xScale,tI1binAvg(StartF,:))
xlabel('Position on main axis 1 in mikrons')
XLimits=get(gca,'Xlim');set(gca,'XMinorTick','on');

YLimits=showRotCell(xScale,XLimits);

pict=StartF;
bup = uicontrol('Parent',fighndl, ...
    'Units','Normalized', ...
    'Callback',['pict=pict+1;pict=min(pict,NofF);',...
    'plot(tP1binx(pict,:)*xScale,tI1binAvg(pict,:));title(num2str(pict));'], ...
    'Position',[0.51 0 0.1 0.05], ...
    'String','>', ...
    'Tag','Pushbutton1');

bdwn= uicontrol('Parent',fighndl, ...
    'Units','Normalized', ...
    'Callback',['pict=pict-1;',...
    'pict=max(pict,1);plot(tP1binx(pict,:)*xScale,tI1binAvg(pict,:));title(num2str(pict));'], ...
    'Position',[0.39 0 0.1 0.05], ...
    'String','<', ...
    'Tag','Pushbutton2');

bst = uicontrol('Style','PushButton',...
    'Parent',fighndl, ...
    'String','Cont.', ...
    'Callback',['guidata(gcbo,0);'], ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','StopButton');

guidata(fighndl,1);
while guidata(fighndl), pause(0.2), end

TEXT='Select center of excitation peak by 1 left button mouse click.';
[ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
from=ileft;to=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
xchannels=to-from+1;    % number of x-bins in w
hold on
plot(a,b,'r')
figure(1)
hold on
plot([a(1) a(1)],YLimits,'m',[a(xchannels) a(xchannels)],YLimits,'m')

figure(2)
TEXT='Select left side region by 1 left button mouse click.';
[ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
fromL=ileft;toL=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
xchannelsL=toL-fromL+1;    % number of x-bins in w
hold on
plot(a,b,'r')
figure(1)
plot([a(1) a(1)],YLimits,'m',[a(xchannelsL) a(xchannelsL)],YLimits,'m')
figure(2)
TEXT='Select right side region by 1 left button mouse click.';
[ileft,iright,a,b]=getlimits(fighndl,w,TEXT);
fromR=ileft;toR=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
xchannelsR=toR-fromR+1;    % number of x-bins in w
plot(a,b,'r')
figure(1)
plot([a(1) a(1)],YLimits,'m',[a(xchannelsR) a(xchannelsR)],YLimits,'m')
figure(2)
hold off
saveas(figure(1),'Selected_Central-Left-Right.fig')
saveas(figure(1),'Selected_Central-Left-Right.jpg')

pause


figure(2) % Plot full intensity distribution for the start frame in blue and w-region in red
plot(tP1binx(StartF,:)*xScale,tI1binAvg(StartF,:))
hold on
plot(tP1binx(StartF,from:to)*xScale,tI1binAvg(StartF,from:to),'r')
plot(tP1binx(StartF,fromL:toL)*xScale,tI1binAvg(StartF,fromL:toL),'r')
plot(tP1binx(StartF,fromR:toR)*xScale,tI1binAvg(StartF,fromR:toR),'r')
xlabel('Position on main axis 1 in mikrons')
hold off
saveas(figure(2),'peak.fig')
saveas(figure(2),'peak.jpg')
%%%%%%%%%%
FullAreaStartF=sum(tI1binAvg(StartF,:)); % Full area of the projected intensity distribution for the 1st frame
medianStartF=sum(tI1binAvg(StartF,:).*tP1binx(StartF,:)) / FullAreaStartF * xScale; %Median position in my of 1st frame



%%%%%%%%%%
NofF=NumberOfFrames-StartF+1; % number of frames
x=tP1binx(StartF:NofF,from:to)*xScale; % distance runs along rows (all columns have identical values)
y=tP1biny(StartF:NofF,from:to)*yScale; % time runs along columns (all rows have identical values)
y=y-y(1,1);                       %y(1,1) contains "yscale" - so this line yields zero start time
z=tI1binAvg(StartF:NofF,from:to); % projected intensity in central region w
A1=sum(z(1,:),2);
zL=tI1binAvg(StartF:NofF,fromL:toL);
zR=tI1binAvg(StartF:NofF,fromR:toR);

t=y(:,1);                     % time vector
% mean amplitude in central region   -------------------------
Amplitude=sum(z,2)/A1;           % area in central red channels relative to 1 frame
AmplitudeL=sum(zL,2)/A1;         % area in left red channels relative to center of 1 frame
AmplitudeR=sum(zR,2)/A1;         % area in right red channels relative tocenter of 1 frame

figure(3)
plot(t,Amplitude,'.')         % decay plot
hold on
plot(t,AmplitudeL,'o')         % decay plot
plot(t,AmplitudeR,'+')         % decay plot
xlabel('Time in secs')
hold off
title(['Ampl./Ampl.(CenterFrame1) (.=center, o=left, +=right) ,' path])
iT1=checkT1 /yScale +1;
iT2=checkT2 /yScale +1;
text(0.6,0.90,['(o L) at ' num2str(checkT1) ' sec: ' num2str(AmplitudeL(iT1))],'units','normalized')
text(0.6,0.85,['(o L) at ' num2str(checkT2) ' sec: ' num2str(AmplitudeL(iT2))],'units','normalized')
text(0.6,0.80,['(. C) at ' num2str(checkT1) ' sec: ' num2str(Amplitude(iT1))],'units','normalized')
text(0.6,0.75,['(. C) at ' num2str(checkT2) ' sec: ' num2str(Amplitude(iT2))],'units','normalized')
text(0.6,0.70,['(+ R) at ' num2str(checkT1) ' sec: ' num2str(AmplitudeR(iT1))],'units','normalized')
text(0.6,0.65,['(+ R) at ' num2str(checkT2) ' sec: ' num2str(AmplitudeR(iT2))],'units','normalized')
saveas(figure(3),'Central-Left-Right.fig')
saveas(figure(3),'Central-Left-Right.jpg')


% End mean amplitude in central region ------------------------

% fraction of central area with respect to full area of start frame --
figure(4)
%[MaxA,where]=max(Amplitude);
fracCentralArea=sum(z,2)/FullAreaStartF;
plot(t,fracCentralArea)
title(['CA/FAStartF-Decay ,' path]) % CA=Central Area divided by FA=Full Area of the StartFrame
xlabel('Time in secs')
saveas(figure(4),'decay')
iT1=checkT1 /yScale +1;
iT2=checkT2 /yScale +1;
text(0.6,0.8,['CA/FAStartF at ' num2str(checkT1) ' sec: ' num2str(fracCentralArea(iT1))],'units','normalized')
text(0.6,0.7,['CA/FAStartF at ' num2str(checkT2) ' sec: ' num2str(fracCentralArea(iT2))],'units','normalized')
%End fraction of central area with respect to full area of start frame

% median ------------------------------
figure(5)
FullArea = sum(tI1binAvg,2);
median=tI1binAvg * tP1binx(StartF,:)' ./ FullArea * xScale - medianStartF; 
        % Intensitätsverteilung I(x) in Zeilen von tI1binAvg 
        % Spaltenvektor = Orte x
plot(t,median)        
title(['Median ' path])
ylabel('Shift from StartF in my')
xlabel('Time in secs')
% End median -------------------------


% Subfunction

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