function DecayFit03
%%
% Version 03 
% lambda fixed
% 3 Fitparameter: Sigma0, D, offset
% Zeifacher Fit für den Fall, dass das Anfangs-D zu schlecht ist.
% Fit mit "echtem" Chiquadrat, dafür: Fehler der Messpunkte aus Polynomfit
% Offset wird im letzten Plot als rote Linie angezeigt
clear all, close all, screensize=get(0,'ScreenSize');
%dummy='All.mat'; %literally a dummy, not used anywhere
cd(pwd);
path=pwd;
%[dummy,path]=uigetfile('*.mat','Select any mat-file in folder to be treated.');
cd(path)
if ~exist('All.mat'); msgbox('No all.mat - returned.'); return, end

L=load('All.mat');
tP1binx=L.tP1binx;
tI1binAvg=L.tI1binAvg;
NumberOfFrames=L.NumberOfFrames;

pausetime=1;
NPixels = 256;
xScale  = 100/NPixels;       % Mikrons / Pixel
yScale=1;                    % Seconds / Frame
w=5;                         % width of region of activation in mikrons
checkT1=10;                  % check-out times in Secs
checkT2=50;
StartF=1;
D=1;                          % Estimate of D for fit
lambda=0.001;                 % fixed

%model = 1;
%%
%button = questdlg('Soma und Tip im Frame? "Cancel" stops processing');
%if strcmp(button,'Cancel'), return, end
%if strcmp(button,'No')
model=0;
%end 

%%
NofF=NumberOfFrames-StartF+1; % number of frames
t=[0:(NofF-1)]' * yScale;         % time (column)vector
fig2hndl=figure(2);
plot(tP1binx(StartF,:)*xScale,tI1binAvg(StartF,:))
xlabel('Position on main axis 1 in mikrons')
XLimits=get(gca,'Xlim');set(gca,'XMinorTick','on');
pause(pausetime)
YLimits=showRotCell(xScale,XLimits); rotcellFigHndl=gcf;
pause(pausetime/2)
if model==1
    TEXT='Select tip edge, then body edge (2 clicks)';
    texthndl=text('Position',[.2,.2],'String',TEXT,...
                  'FontSize',10,'Units','normalized');
    axis('equal');drawnow
    [xc,yc]=ginput(2);
    xTip=xc(1), xBody=xc(2), abs(xTip-xBody)
end
% pict=StartF;
% CBD=struct('pict',pict,'NofF',NofF,'tP1binx',tP1binx,'tI1binAvg',tI1binAvg,'xScale',xScale); % Callbacks der Knöpfe benötigen diese Daten
% guidata(fig2hndl,CBD); % CallBackData in Datenbereich der figure
% evalin('base','CBD=guidata(gcf);'); % CBD im Base-Bereich für Callbacks erzeugen
% bup = uicontrol('Parent',fig2hndl, ...
%     'Units','Normalized', ...
%     'Callback',['CBD.pict=CBD.pict+1;CBD.pict=min(CBD.pict,CBD.NofF);',...
%     'plot(CBD.tP1binx(CBD.pict,:)*CBD.xScale,CBD.tI1binAvg(CBD.pict,:));title(num2str(CBD.pict));'], ...
%     'Position',[0.51 0 0.1 0.05], ...
%     'String','>', ...
%     'Tag','Pushbutton1');
% bdwn= uicontrol('Parent',fig2hndl, ...
%     'Units','Normalized', ...
%     'Callback',['CBD.pict=CBD.pict-1;',...
%     'CBD.pict=max(CBD.pict,1);plot(CBD.tP1binx(CBD.pict,:)*CBD.xScale,CBD.tI1binAvg(CBD.pict,:));title(num2str(CBD.pict));'], ...
%     'Position',[0.39 0 0.1 0.05], ...
%     'String','<', ...
%     'Tag','Pushbutton2');
% bst = uicontrol('Style','PushButton',...
%     'Parent',fig2hndl, ...
%     'String','Cont.', ...
%     'Callback',['[hbutton,hfigure]=gcbo; uiresume(hfigure);'], ...
%     'Units','normalized', ...
%     'Position',[0 0 0.1 0.1], ...
%     'Tag','StopButton');
% uiwait(fig2hndl);
figure(fig2hndl)
TEXT='Select center of excitation peak by 1 left button mouse click.';
[ileft,iright,a,b]=getlimits(fig2hndl,w,TEXT);
from=ileft;to=iright;   % from = index of data point center-w/2 ; to = index of data point center+w/2
xchannels=to-from+1;    % number of x-bins in w
hold on, plot(a,b,'r'), pause(pausetime)

%% Zeige Grenzen in der RotCell
figure(rotcellFigHndl)  
hold on
plot([a(1) a(1)],YLimits,'m',[a(xchannels) a(xchannels)],YLimits,'m') %plot vertical lines at "center +/- w/2"
pause(pausetime)

%%  extract sigma0 from start frame by gauss fit
figure(3) 
Xs=tP1binx(StartF,:)*xScale; Ys=tI1binAvg(StartF,:);
% Normalize
dx=abs(Xs(1)-Xs(2));
AreaSF=sum(Ys)*dx;
Ys=Ys/AreaSF;
plot(Xs,Ys,'.') %Plot full intensity distribution for the start frame in blue and w-region in red
xlabel('Position on main axis 1 in mikrons')
hold on
plot(Xs(from:to),Ys(from:to),'r.'), pause(pausetime)
    
%% ******** Fitting start frame Ys(Xs) ************************
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
Pout=fminsearch(@(P) chi(P,Xs,Ys),Pstart, optimset('Display','iter'))
sigma0=Pout(1);
offsetSF=Pout(3);
ActPos=Pout(2); % Distance from tip in µm
   % show the result
hL=plot(Xs,Gsimu(Pout,Xs),'g');
hT=title(['Fitted curve, Sigma0 = ', num2str(sigma0),' . Paused for 5s.']);
pause(pausetime)
%hold off %figure(3)
saveas(3,'StartFramelambdafix.jpg')

%% zeichne Iavg(t) mit Fehlern, Anzeige der Werte zu den Check-out-Zeiten
Ic=tI1binAvg(StartF:NofF,from:to)/AreaSF; % Projected intensity in central region w, normalized    
CentralArea=sum(Ic,2)*dx; 
Iavg = CentralArea/abs(Xs(from)-Xs(to)); % Averaged intensity in region w   
fehler=Polynomfit(t,Iavg); % Polynomfit zur Bestimmung der Fehler der Messpunkte
figure(4), errorbar(t,Iavg,fehler,'.'), hold on % plot Iavg(t)
title(['Center Int., avg over w ,' path]), xlabel('Time in secs')
iT1=checkT1 /yScale +1; iT2=checkT2 /yScale +1;
text(0.2,0.85,['CA/FAStratF at ' num2str(checkT1) ' sec: ' num2str(CentralArea(iT1))],'units','normalized')
text(0.2,0.8, ['CA/FAStratF at ' num2str(checkT2) ' sec: ' num2str(CentralArea(iT2))],'units','normalized')
text(0.2,0.75,cd,'units','normalized','interpreter','none')


%% ********** Fit von Iavg(t)
Pstart=[sigma0 D offsetSF ];   % parameter estimates, partly from start frame (see above)  
if model==1                    % Vollständiges Modell mit Images, braucht Positionen von Tip, Body und Aktivierungsstelle
    HP=[ActPos, xTip, xBody, lambda, fehler']; %constant parameters
    hL=plot(t,Dsimu(Pstart,t,HP),'g'); % show start estimates
    hT=title('Start curve. Paused for 1s.');pause(pausetime)
    % -- do the fit --
    delete(hT); delete(hL); % delete start estimates
    Pout=fminsearch(@(P) Dchi(P,t,Iavg,HP),Pstart, optimset('Display','iter'))
    Pstart=Pout; % nochmal, falls Startwert zu weit vom besten Wert entfernt war
    Pout=fminsearch(@(P) Dchi(P,t,Iavg,HP),Pstart, optimset('Display','iter'))
    D=Pout(2); chisquare=Dchi(Pout,t,Iavg,HP);
else                           % einfaches Modell
    HP=[0,0,0, lambda, fehler']; %constant parameters
    hL=plot(t,DsimuX(Pstart,t,HP),'g'); % show start estimates
    hT=title('Start curve. Paused for 1s.'); pause(pausetime)
    % -- do the fit --
    delete(hT); delete(hL); % delete start estimates
    Pout=fminsearch(@(P) DchiX(P,t,Iavg,HP),Pstart, optimset('Display','iter'))
    Pstart=Pout; % nochmal, falls Startwert zu weit vom besten Wert entfernt war
    Pout=fminsearch(@(P) DchiX(P,t,Iavg,HP),Pstart, optimset('Display','iter'))
    D=Pout(2); chisquare=DchiX(Pout,t,Iavg,HP);
end


Soll=length(t)-length(Pout);

%% show the result
figure(4)
plot([0 max(t)], [Pout(3) Pout(3)],'r') % Offset-Linie
if model==1
    F=Dsimu(Pout,t,HP); hL = plot(t,F,'g');
end
F = DsimuX(Pout,t,HP);
hL = plot(t,F,'r--');   %Kurve ohne Images
hT = title(['Icenter(t), D = ', num2str(D),' .']);
axis([0 max(t)*1.05 Pout(3)-abs(Pout(3))*0.1 max(max(Iavg),max(F))*1.1]);
text(0.2,0.9,['Chi^2: ' num2str(chisquare) ' Soll: ' num2str(Soll)],'units','normalized')
hold off           %figure(4)
set(4,'Position',screensize);
saveas(4,'NewFitlambdafix.fig'),saveas(4,'NewFitlambdafix.jpg')
disp(num2str(D))
dlmwrite('..\Deff_auto.csv',num2str(D),'delimiter','','-append')
[~, ParentFolderName] = fileparts(pwd);
outvar = [ParentFolderName, ';', num2str(D)];
dlmwrite('..\name_Deff_auto.csv',outvar,'delimiter','','-append')
%% --------------------------------------------------------------------------
% Subfunctions

function dev=Polynomfit(x,y)
zentr=mean(x); % Fit mit Polynom 4. Ordnung zur Ermittlung der Fehler der Messpunkte
skala=std(x);   %Skalieren und Zentrieren
z=(x-zentr)/skala;
[p,S]=polyfit(z,y,4);
[Y,dev]=polyval(p,z,S); % dev enthält die Fehler (im statistischen Sinn, s. Matlab-Hilfe)
figure(44), errorbar(z,y,dev,'.'); hold on, % zur visuellen Kontrolle
plot(z,Y,'g-'); 
title('Fit mit Polynom 4. Ordnung zur Fehlerermittlung');
pause(1)
    
function [ileft,iright,x,y]=getlimits(fig2hndl,w,TEXT)
% function [ileft,iright,x,y]=getlimits(fig2hndl)
% Lets you select the central position of the excitation peak on a line plot in figure with
% handle fig2hndl. From this position, 2 positions  w/2 to the left and w/2 to the rigth are calculated.
% Returns indices ileft and iright of the data at these
% positions, as well as the x- and y-values within that range.
set(gcf,'pointer','fullcrosshair');
texthndl=text('Position',[.2,.2],'String',TEXT,...
        'FontSize',10,'Units','normalized');
drawnow;
[xc,yc]=ginput(1);
linehndl=findobj(fig2hndl,'type','line','-and','color','blue');
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
CLineSet=[5 10 30 60 90 200];
if exist('CLinSet.mat','file'), load('CLinSet.mat'), end
rotcellh=openfig('RotCell.fig');
HG=findobj(rotcellh,'type','hggroup');
X=get(HG,'xdata');
Y=get(HG,'ydata');
Z=get(HG,'zdata');
Xs=X*xScale;Ys=Y*xScale;
contour(Xs,Ys,Z, CLineSet);
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
S=Dsimu(P,x,HP(1:4));
fehler=HP(5:end);
out=sum(((y-S)./fehler').^2);
strafe=abs(P(3))-0.05*y(1); % offset auf ca. 5% des ersten Iavg-Wertes begrenzen
strafe=2+tanh(10*strafe/(0.05*y(1)));
out=out*strafe;

function S=Dsimu(P,x,HP)
sigma0=P(1);lambda=HP(4); %vollständiges Modell mit Images
ActPos=HP(1); xTip=HP(2); xBody=HP(3);
sigma=sqrt(2*P(2)*x + sigma0^2);
S =exp(-lambda*x) / sqrt(2*pi)./ sigma .*(1+exp(-4*(ActPos-xTip)^2 ./ (2*sigma.^2))...   
                                   - exp(-4*(ActPos-xBody)^2 ./ (2*sigma.^2))) +P(3) ;

function out=DchiX(P,x,y,HP)
S=DsimuX(P,x,HP(1:4));
fehler=HP(5:end);
out=sum(((y-S)./fehler').^2);
strafe=abs(P(3))-0.05*y(1); % offset auf ca. 5% des ersten Iavg-Wertes begrenzen
strafe=2+tanh(10*strafe/(0.05*y(1)));
out=out*strafe;
                               
                               
function S=DsimuX(P,x,HP)
sigma0=P(1);lambda=HP(4); %Test: Ohne Images
sigma=sqrt(2*P(2)*x + sigma0^2);
S = exp(-lambda*x)/sqrt(2*pi) ./ sigma  + P(3);          % 


