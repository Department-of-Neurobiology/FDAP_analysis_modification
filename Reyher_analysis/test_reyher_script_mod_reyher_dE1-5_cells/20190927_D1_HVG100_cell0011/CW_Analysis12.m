function CW_Analysis12
%% Version 11b, addition / changes:
% 1)    If dark frame is present, the intensity distribution on the main axis of the dark frame is subtracted
%       from those of subsequent frames
% 2)    The mask M, describing the ROI, is stored in all.mat
% 3)    The dark frame also stored
% 4)    Interaktive Einstellung der Intensitätsschwelle (Definiert Zellrand)
% 5)    ROI mit 6 Eckpunkten
% 6)    Restrukturierung (Mehr Kode in subfunctions)
% 7)    Aliasing-Problem abgemildert durch Anpassung für optimale binwidth1 = Pout(2) und bin-Phase = Pout(1) bei Frame 1.
close all, clear all,screensize=get(0,'ScreenSize');
load('mymaps')
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%% to be adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DarkFramePresent=1; % 0 = not present, 1 = present
[filen,cpath]=uigetfile('*.tif','Select the LAST tif-file of series.');
cd(cpath);
drawnow;
[fn1, reminder]=strtok(filen,'_');fn1=[fn1 '_'];
NOP=strtok(reminder,'.');
digits=length(NOP)-1;
NumberOfFrames=str2num(NOP(2:digits+1))+1-DarkFramePresent;
fn3='.tif';                       % Extension
idx=strfind(cpath,'\');n=length(idx);
ID=cpath(idx(n-1)+1:idx(n)-1);                      % Identification. For figure captions
colorindex=2;                     % Color used in tif-file, 1 = red, 2 = green, 3 = blue
rf = 2;                           % radius of Gaussian filter (larger values smooth more but blur details)
%% Gaussfiltern
CLineSet=[5 10 15 30 40 50 100 300 500 1000 2000]; % provisorische Kontourlinien
ff=GaussFilter(rf);
button = questdlg('Apply or re-apply filtering? "Cancel" stops processing');
if strcmp(button,'Cancel'), return, end
if strcmp(button,'Yes')
    filterung(fn1,fn3,ID,digits,DarkFramePresent,NumberOfFrames,CLineSet,ff,colorindex);
end 

%% ----------------- Reference Intensity -------------------------
button = questdlg('Process or re-process reference area? "Cancel" stops processing');
if ~exist(['fi_' num2str(NumberOfFrames-1) '.mat'])
    msgbox('No filtered data available - stop processing')
    return
end %if ~exist
if strcmp(button,'Cancel'), return, end
if strcmp(button,'Yes')
   Referenz(digits,DarkFramePresent,NumberOfFrames,CLineSet)
end %if strcmp

%% ------------------ Main Analysis Start -------------------------------
[intlimit,CLineSet,M,Ax,sze]=FindIntLimROI(digits,DarkFramePresent,greenmap); 
                                                              % Interaktive Einstellung der Intensitäts-Schwelle intlimit,
                                                              % durch welche der Zellrand definiert wird. Weiterhin interaktive Festlegung
                                                              % der ROI. ROI gekennzeichnet durch Maske M. M ist gleich 2 innerhalb des ROI, 
                                                              % und 0 außerhalb.
save('CLinSet.mat','CLineSet');
ymax=sze(1); xmax=sze(2); % Bildhöhe und -breite

% Allocation
CGcell=zeros(2,NumberOfFrames+DarkFramePresent);       % storage array for Center of Gravity coordinates
IntTotalCell=zeros(1,NumberOfFrames+DarkFramePresent); % total intensity in the region of interest for Intensity > intlimit
U=cell(NumberOfFrames+DarkFramePresent,3); % initialize cell array to hold projection data
npmax=0;
% End Allocation

% REFERENCE AREA business
IntMean=ones(1,NumberOfFrames+DarkFramePresent);          % dummy if no reference area exists
if exist('MeanIntensityRefArea.mat')
    load('MeanIntensityRefArea.mat');       % load array IntMean
    figure(2), plot(IntMean), title(['Mean Intensity per pixel of reference area versus frame number -' ID])
    msgbox('Paused for printing'), pause
end %if exist
MeanRefI=sum(IntMean)/NumberOfFrames;

% figure preparation
fHandle=figure(4); set(fHandle,'Renderer','Painters'); guidata(fHandle,1);
b = uicontrol('Style','PushButton',...
    'Parent',fHandle, ...
    'String','STOP', ...
    'Callback',['guidata(gcbo,0);'], ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','StopButton');

%% ***********************    main loop  over pict  ***************************
%  *******************************************************************
format=['%0' num2str(digits) 'u'];
for pict=1:NumberOfFrames+DarkFramePresent
    fn2=sprintf(format,pict-1+DarkFramePresent);
    if pict == NumberOfFrames+1
        fn2='dark';
    end
    load(['fi_' fn2 '.mat']);     % load filtered data
    fz=fz*MeanRefI/IntMean(pict); % Correcting for constant intensity of the reference area
    %%%%%%%%%%%%%%%%%%%%%%%%% fz (= filtered z) containes the smoothed data with minimum = 0, array size is (1 to nsx, 1 to nsy) %%%%%%%%%%%%%%%%%%%%%%%%%
    z=(M>1).*(fz+eps);            % select ROI using  mask M; z is zero outside ROI and at least eps inside (eps is 2.2e-16 or so)
    
    if pict == 1
        %newsize=size(z); newymax=newsize(1); newxmax=newsize(2); % überflüssig in Schleife - z hat konstante Grösse
        [yInROI,xInROI]=find(z > 0);                       % find x and y coordinates where z > 0
        NinROI=length(yInROI);                             % Number of pixels in ROI
    end %if pict
    
    STEP=(z>intlimit).*z;                              % getting a 2D step function for z > intlimit
    [yInWindow,xInWindow]=find(STEP > 0);              % find x and y coordinates where z > intlimit
    NinWindow=length(yInWindow);
    % theses pixels are in the window treated in the following
    xCG=0;yCG=0;IntTotal=0;  % Reset of summing variables
    % --------------------
    %
    for j=1:NinWindow
        IntPixel_j=z(yInWindow(j),xInWindow(j));
        IntTotal=IntTotal+IntPixel_j;    % Total intensity of ROI
        xCG=xCG+xInWindow(j)*IntPixel_j; % find the center of intensity for points with  fz > intlimit
        yCG=yCG+yInWindow(j)*IntPixel_j;
    end
    xCG=xCG/IntTotal;                      % Center of intensity coordinates for actual frame (with number 'pict')
    yCG=yCG/IntTotal;
    CGcell(1,pict)=xCG;CGcell(2,pict)=yCG; % Keep center of intensity coordinates in array CGCell
    %    ------------------------------------
    IntTotalCell(pict)=IntTotal/NinWindow; % Keep total intensity in the region of interest for Intensity > intlimit (per pixel)
    %    --------------------------------------
    figure(4), subplot(1,2,1)
    contour(z,CLineSet);      % make a  contour plot
    axis(Ax);  
    hold on      %''''''''''''''''''''''''''''''''''''' hold on
    plot([xCG xCG], [1 ymax],'-k') % Plotting center of intensity as a cross into contour plot
    plot([1 xmax], [yCG yCG],'-k')
    
    % ------------------ Momentum of Intensity
    T=zeros(2,2);
    for j=1:NinWindow                                       % loop over all pixels in current window (defined by "z>intlimit")
        IntPixel_j=z(yInWindow(j),xInWindow(j));
        T(1,1)=T(1,1)+yInWindow(j)^2*IntPixel_j;            % SUM(y^2*Intensity), diagonal element
        T(2,2)=T(2,2)+xInWindow(j)^2*IntPixel_j;            % SUM(x^2*Intensity), diagonal element
        T(1,2)=T(1,2)-yInWindow(j)*xInWindow(j)*IntPixel_j; % SUM(x*y*Intensity), off-diagonal element
    end
    T(2,1)=T(1,2);                 % its a symmetric tensor
    Ts(1,1)=T(1,1)-IntTotal*yCG^2; % now tensor relative to center of intensity (Steiner's theorem)
    Ts(2,2)=T(2,2)-IntTotal*xCG^2;
    Ts(1,2)=T(1,2)+IntTotal*yCG*xCG;
    Ts(2,1)=Ts(1,2);
    
    if pict==1                    % keep main axes of 1st frame for all frames
        [V,dummy]=eig(Ts);           % eig returns eigenvectors in V, and eigenvalues in dummy
        v1=V(:,1);                   % Eigenvector 1 = directions of main axes with smaller eigenvalue ( plotted in red below)
        v2=V(:,2);                   % Eigenvector 2 = directions of main axes with bigger eigenvalue ( plotted in blue)
        sx=sign(v1(1)); sy=sign(v1(2));        % rotating and/or reflecting eigenvectors
        v1=v1*sy; % points now always upwards
        sx=sign(v2(1)); sy=sign(v2(2));
        v2=v2*sy; % points now always upwards
        RMA=[v1(:)' ; v2(:)'];        % Rotationsmatrix
        
        pMarker=600:-50: -600;         % Marker lines of position along main axis 1
        ML=4000;
        if abs(v1(1)) >= 0.707
            xM=pMarker/v1(1);
            for S=1:length(xM)
                plot([xM(S) xM(S)+v2(1)*ML],[0 0+v2(2)*ML],'k')
            end
        else
            v2=-v2;
            yM=pMarker/v1(2);
            for S=1:length(yM)
                plot([0 0+v2(1)*ML],[yM(S) yM(S)+v2(2)*ML],'k')
            end
        end %if
        RROI(RMA,z,CLineSet) % Display and/or save rotated cell
    end %if pict==1
    
    %-----------------------------------------
    for S=-200:400:200           % plot main axes of momentum of intensity tensor
        plot([xCG xCG+v1(1)*S], [yCG yCG+v1(2)*S],'-r')
        plot([xCG xCG+v2(1)*S], [yCG yCG+v2(2)*S],'-b')
    end
    %-----------------------------------------
    if abs(v1(1)) >= 0.707   % Marker lines of position along main axis 1
        xM=pMarker/v1(1);
        for S=1:length(xM)
            plot([xM(S) xM(S)+v2(1)*ML],[0 0+v2(2)*ML],'k')
        end
    else
        yM=pMarker/v1(2);
        for S=1:length(yM)
            plot([0 0+v2(1)*ML],[yM(S) yM(S)+v2(2)*ML],'k')
        end
    end %if
    axis(Ax);
    if pict == NumberOfFrames+1
        title([ID ' ROI only,  Dark Frame' ])
    else
        title([ID ' ROI only,  Frame: ' num2str(pict)])
    end
    axis('square')
    %     %TEST
    %         for j=1:NinROI
    %             plot(xInROI(j),yInROI(j),'.y')
    %         end
    %         pause
    %     %END TEST
    
    hold off %''''''''''''''''''''''''''''''''''''''''''''''''' hold off 
    %-------------------------------------------------
    % Calculating projections on main axes
    if pict == 1                                % find binwidth and scale range for both axes for all frames
        P1=zeros(NinROI,1); P2=P1;              % all values are in "pixel units"
        for j=1:NinROI
            p1=v1.*[xInROI(j);yInROI(j)]; P1(j)=sum(p1); % projection of all pixel vectors in ROI on main axis 1 (Skalarprodukt)
            p2=v2.*[xInROI(j);yInROI(j)]; P2(j)=sum(p2); % projection on main axis 2
        end
        low1=min(P1); high1=max(P1); 
        binwidth1=(high1-low1)/100;         % 100 bins (channels) along main axis
        %low2=min(P2); high2=max(P2); % axis 2 is not considered, so far
        %binwidth2=(high2-low2)/100;         
    end %if pict==1
    IntPixel=zeros(NinWindow,1); P1=IntPixel; P2=P1;
    for j=1:NinWindow
        IntPixel(j)=z(yInWindow(j),xInWindow(j));          % Intensität 
        p1=v1.*[xInWindow(j);yInWindow(j)]; P1(j)=sum(p1); % projection on main axis 1
        p2=v2.*[xInWindow(j);yInWindow(j)]; P2(j)=sum(p2); % projection on main axis 2
    end
    if pict ==1                     % Finde optimale Binwidth
        Pstart=[0 binwidth1];
        Pout=fminsearch(@(P) optbin(P,P1,IntPixel,low1,high1),Pstart,optimset('Display','iter'))
        low1=low1-Pout(1); binwidth1=Pout(2);
    end
    
    %***********************************
    [P1bin,I1binAVG,I1binSUM]=put2bins(P1,IntPixel,low1,high1,binwidth1,0.0); % subfunction put2bins puts the actual value P1 into the appropriate bin
    %***********************************
    % [P2bin,I2binAVG,I2binSUM]=put2bins(P2,IntPixel,low2,high2,binwidth2,0.0); % P1bin(n)=low1+(n-1/2)*binwidth
    % axis 2 is not considered, so far      
    np=length(P1bin);
    npmax=max([np npmax]);
    U{pict,1}=P1bin;               % hold the projection data in cell array U
    U{pict,2}=I1binSUM;
    U{pict,3}=I1binAVG;
    figure(4);subplot(1,2,2); plot(P1bin,I1binAVG); axis('square')
    %   pause % to get the possibility of interruption by CTRL-C
    drawnow
    %     printing=1:50:400;             % print every 50th frame
    %     if sum(pict==printing)
    %         disp('Paused for printing. Press any key to resume processing.')
    %         pause
    %     end%if sum
    if ~guidata(fHandle), return, end         % Stop button pressed.
end %for pict=... main loop
%**************************************************************************
%********************       END LOOP over pict        *****************************
%**************************************************************************
%**************************************************************************
%%
save('u.mat','U','NumberOfFrames','npmax','ff')

%% ---- generating a 3d intensity surface z(x,y), with z= intensity, x=position along main axis, y= frame number
% allocating arrays
    tI1binSum=zeros(NumberOfFrames,npmax); % Total intensity in bins
    tI1binAvg=tI1binSum;   % Average intensity in bins
    tP1binx=tI1binSum;     % "x"-coordinate along main axis of bins
    tP1biny=tI1binSum;     % "y"-coordinate = frame number
% end allocation

button = questdlg('Do you want to correct for constant mean intensity? "Cancel" stops processing');
if strcmp(button,'Cancel'), return, end
swtch=0;
if strcmp(button,'Yes')
    AvgTotalInt=sum(IntTotalCell)/NumberOfFrames;
    swtch=1;
end %if strcmp(button,'Yes')

I1binSUM_darkframe=zeros(1,npmax);
I1binAVG_darkframe=zeros(1,npmax);
if DarkFramePresent
    I1binSUM_darkframe=U{NumberOfFrames+DarkFramePresent,2};
    I1binAVG_darkframe=U{NumberOfFrames+DarkFramePresent,3};
end

for pict=1:NumberOfFrames
    P1bin   =U{pict,1};
    I1binSUM=U{pict,2}-I1binSUM_darkframe;   % Abzug des darkframes
    I1binAVG=U{pict,3}-I1binAVG_darkframe;   %
    figure(33); plot(P1bin,I1binAVG_darkframe,'k',P1bin,U{pict,3},'g');
     title(['Frame ' num2str(pict) '. Black=DarkFrame. Paused(0.1s)']); pause(0.1) %## Diagnose ##
    if swtch == 1
        I1binSUM=I1binSUM*AvgTotalInt/IntTotalCell(pict);
        I1binAVG=I1binAVG*AvgTotalInt/IntTotalCell(pict);
    end %if swtch
    np = length(P1bin);
    tI1binSum(pict,1:np)=I1binSUM;
    tI1binAvg(pict,1:np)=I1binAVG;
    tP1binx(pict, 1:np) =P1bin;
    tP1biny(pict, :) =pict*ones(1,npmax);
end %for

rf=1; ff=GaussFilter(rf); % Milder Filter
tI1binSumValid=filter2(ff,tI1binSum,'valid'); % Filter 2D-Data
tI1binSum(rf+1:NumberOfFrames-rf,rf+1:npmax-rf)=tI1binSumValid;   % replace "valid" region with filtered data
fHandle=figure(5);
contourf(tP1binx,tP1biny,tI1binSum,30)
hcontourgroup=findobj(gca,'Type','hggroup');
set(hcontourgroup,'LineStyle','none');

title(['Total intensity in stripes perpendicular to main axis (DarkF subtr.) ' ID])
set(gca,'XDir','reverse');view(180,90); colormap(mymap);
set(fHandle,'Renderer','Painters')
saveas(fHandle,'TotalIntProjection.fig')
saveas(fHandle,'TotalIntProjection.jpg')

tI1binAvgValid=filter2(ff,tI1binAvg,'valid');
tI1binAvg(rf+1:NumberOfFrames-rf,rf+1:npmax-rf)=tI1binAvgValid; % replace "valid" region with filtered data
fHandle=figure(6);
contourf(tP1binx,tP1biny,tI1binAvg,30)
hcontourgroup=findobj(gca,'Type','hggroup');
set(hcontourgroup,'LineStyle','none');
set(gca,'XDir','reverse')
title(['Mean intensity per pixel in stripes perpendicular to main axis (DarkF subtr.)' ID])
view(180,90); colormap(mymap);
set(fHandle,'Renderer','Painters')
saveas(fHandle,'MeanIntProjection.fig')
saveas(fHandle,'MeanIntProjection.jpg')

XL=get(gca,'Xlim');
fHandle=openfig('RotCell.fig');
set(gca,'Xlim',XL);
set(gca,'YlimMode','auto');
YL=get(gca,'Ylim'); YW = YL(2) - YL(1);
set(gca,'Ylim',[YL(1)-0.2*YW YL(2)+0.2*YW]);
saveas(fHandle,'RotCell.fig')
saveas(fHandle,'RotCell.jpg')

fHandle=figure(7);
plot(IntTotalCell)
title(['Mean intensity per pixel in ROI ' ID] )

saveas(fHandle,'MeanIntROI-decay')
saveas(fHandle,'MeanIntROI-decay.jpg')

fig8=figure(8); % displays cross sections of Intensity(x,y) along x-axis
pict=1; plot(tP1binx(pict,:),tI1binAvg(pict,:))
ax=axis;

save('All.mat','pict','NumberOfFrames','ax','tP1binx','tP1biny','tI1binAvg','tI1binSum','-mat')
BinCoordinatesPixels=tP1binx(1,:);
save('BinCoordinatesInPixels.txt','BinCoordinatesPixels','-ASCII')
save('totalInt_inBinsAvg.txt','tI1binAvg','-ASCII')
save('totalInt_inBinsSum.txt','tI1binSum','-ASCII')
uicontrol('Parent',fig8, ...
    'Units','points', ...
    'Callback',['if ~exist(''pict'');load(''All.mat'');end;pict=pict+1;pict=min(pict,NumberOfFrames);',...
    'plot(tP1binx(pict,:),tI1binAvg(pict,:));axis(ax);title(num2str(pict));'], ...
    'Position',[370 25 20 20], ...
    'String','>', ...
    'Tag','Pushbutton1');
uicontrol('Parent',fig8, ...
    'Units','points', ...
    'Callback',['if ~exist(''pict'');load(''All.mat'');end;pict=pict-1;',...
    'pict=max(pict,1);plot(tP1binx(pict,:),tI1binAvg(pict,:));axis(ax);title(num2str(pict));'], ...
    'Position',[330 25 20 20], ...
    'String','<', ...
    'Tag','Pushbutton2');
%pause

%%
%--------------------------------------------------------- subfunctions ------------------------------------------------
%
function M=setMask(fz,fHandle)
sze=size(fz);
xmax=sze(2);ymax=sze(1);
figure(fHandle)
text(0.4,-0.085,'Define ROI by 6 corner points!','units','normalized','color','r');
hold on
N=6;
% Get N points defining mask area
rx = zeros(N,1);
ry = rx;
n = 0;
% Loop, picking up the points.
%disp('Left mouse button picks 4 points.')
while 1
    [xi,yi,button] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    rx(n)=xi; ry(n)=yi;
    if n==N, break, end
end
% center of mask area
X=sum(rx)/N; Y=sum(ry)/N;
% corner points with respect to center
for k=1:N
    ux(k)=rx(k)-X;
    uy(k)=ry(k)-Y;
    [theta,u(k)]=cart2pol(ux(k),uy(k)); %polar coordinates of corner points
    if theta < 0, theta=theta+2*pi; end
    t(k)=theta;
end
% sorting by decreasing theta
[dummy,index]=sort(-t);
rx=rx(index);
ry=ry(index);
t=t(index);
u=u(index);
ux=ux(index);
uy=uy(index);
% Vectors connecting corner points. Vectors are pointing clockwise around mask area
dx(1:N-1)=ux(2:N)-ux(1:N-1);
dx(N)=ux(1)-ux(N);
dy(1:N-1)=uy(2:N)-uy(1:N-1);
dy(N)=uy(1)-uy(N);

% Define mask M
M=zeros(ymax,xmax);
% Loop to set the border line in the mask
step=1e-3;
for k=1:N
    S=norm([dx(k),dy(k)]);
    for kk=0:step:1
        Rx=rx(k)+kk*dx(k);
        Ry=ry(k)+kk*dy(k);
        I=round(Rx);
        J=round(Ry);
        M(J,I)=2;
    end
end
for k=1:ymax
    [dummy, index] = find(M(k,:)==2);
    from = min(index); to = max(index);
    M(k,from:to)=2;
end

%
%--------------------------------------------------------- subfunction put2bins ------------------------------------------------
%
function [Bbin,AbinAVG,AbinSUM] = put2bins(B,A,low,high,binWidth,pad)
% function [Bbin,AbinAVG,AbinSUM] = put2bins(B,A,binWidth,pad)
% B ist unabhängige, A ist abhängige Variable, also A(k)=f(B(k)). B-Werte sind nicht äquidistant.
% put2bins ordnet jeden B-Wert einem Intervall (=bin) der Breite binWidth auf einer äquidistanten Skala zu.
% Diese Skala läuft von min(B) bis max(B) in Schritten von binWidth und
% wird in Bbin zurückgegeben. Die zu einem Intervall j gehörenden A-Werte
% werden aufsummiert nach AbinSUM(j) gelegt.
% AbinAVG ist AbinSUM/(Anzahl der A-Werte im Intervall).
% pad gibt Breite eines Einbettbereiches E an, in den Bbin
% gelegt wird. Bbin ist [E B E] mit E=pad*(max(B)-min(B)).
% pad darf gleich Null sein.
Bmin=low;
nB=floor((B-Bmin)/binWidth) +1;       % nB ordnet jedem B-Wert eine bin- (oder Kanal-) Nummer zu
nBstop=floor((high-Bmin)/binWidth) +1;  % Die Intervalle laufen von 1 bis nBstop. Nicht jedes Bin kommt vor,
                                        % d.h. nB enthält i.A. nicht jede Zahl von  1 bis nBstop
kB=1:nBstop;                    % Index aller Bins
AbinAVG=kB*0;
AbinSUM=AbinAVG;
for k=kB                           % Schleife vektorisieren = ?
    binmembers=find(nB==k);
    AbinSUM(k)=sum(A(binmembers));   % Summe von Ereignissen
    if ~isempty(binmembers)           % Mittelwert über Bin-Mitglieder
        AbinAVG(k)=AbinSUM(k)/length(binmembers);
    end %if
end %for
Bbin=Bmin + (kB-1/2)*binWidth;     % Der Wert eines Intervalls soll in seiner Mitte liegen
% padding
nE=floor(nBstop*pad);
if nE>0
    ne=1:nE;
    Bmax=max(Bbin);
    Bbin=[Bmin-(nE-ne+1/2)*binWidth Bbin Bmax+ne*binWidth];
    AbinAVG=[zeros(1,nE) AbinAVG zeros(1,nE)];
    AbinSUM=[zeros(1,nE) AbinSUM zeros(1,nE)];
end

function out=optbin(P,B,A,Bmin,high)  % suche optimales binwidth
Bmin=Bmin-P(1); binWidth=P(2);
nB=floor((B-Bmin)/binWidth) +1;       % nB ordnet jedem B-Wert eine bin- (oder Kanal-) Nummer zu
nBstop=floor((high-Bmin)/binWidth) +1;  % Die Intervalle laufen von 1 bis nBstop. Nicht jedes Bin kommt vor,
                                        % d.h. nB enthält i.A. nicht jede Zahl von  1 bis nBstop
kB=1:nBstop;                    % Index aller Bins
AbinAVG=kB*0;
AbinSUM=AbinAVG; %out=0;
for k=kB                           % Schleife vektorisieren = ?
    binmembers=find(nB==k);
    AbinSUM(k)=sum(A(binmembers));   % Summe von Ereignissen
    if ~isempty(binmembers)           % Mittelwert über Bin-Mitglieder
        AbinAVG(k)=AbinSUM(k)/length(binmembers);
    end %if
    %out=out+AbinAVG(k)*(-1)^k;
end %for
out=sum(abs(diff(AbinAVG)));
if binWidth > 1, out=out*((binWidth-1)/0.05)^2; end % Strafe für zu große binwidth



%--------------------------------------------------------- subfunction RotateROI ------------------------------------------------

function RROI(RMA,z,CLineSet)
[a,b]=size(z);
[X,Y]=meshgrid(1:a,1:b);
XP=zeros(b,a);YP=XP;
for i=1:b
    for j=1:a
        V= RMA * [X(i,j);Y(i,j)];
        XP(i,j)=V(1);
        YP(i,j)=V(2);
    end
end
fHandle=figure(20);
contour(XP,YP,z, CLineSet)
saveas(fHandle,'RotCell.fig')
close(fHandle);

%----------------------------------------- subfunction Gaussfilter 
function  ff=GaussFilter(rf)
    % set up filter function, Gaussian filter with radius rf
    [xf,yf]= meshgrid(-rf:rf);
    ff = -xf.*xf-yf.*yf;
    ff = exp(ff/3/rf);
    ff = ff/sum(sum(ff));
%---------------------------------- subfunction filterung
function filterung(fn1,fn3,ID,digits,DarkFramePresent,NumberOfFrames,CLineSet,ff,colorindex)
% filtering, subtraction of minimum (offset)
fHandle=figure(1);
guidata(fHandle,1);
b = uicontrol('Style','PushButton',...
    'Parent',fHandle, ...
    'String','STOP', ...
    'Callback',['guidata(gcbo,0);'], ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','StopButton');
format=['%0' num2str(digits) 'u'];
for pict=1:NumberOfFrames+DarkFramePresent
    fn2=sprintf(format,pict-1);
    fn=[fn1 fn2 fn3];
    z=imread(fn);
    sze=size(z);    xmax=sze(2);         ymax=sze(1);
    Z=z(ymax:-1:1,1:xmax,colorindex); % in viewer programs y=0 is top, in matlab  y=0 is bottom
    clear z; % free memory from the very large array z
    IdxSat=find(Z>4094);
    SatPix=length(IdxSat);
    Z=double(Z); % raw data are now in double precision array "Z"
    %-----------------------------
    fz=filter2(ff,Z,'vaild');% now smoothing Z with f
    a=min(min(fz));fz=(fz-a); % now getting new minimum and subtract it
    if DarkFramePresent && pict==1
        fn2='dark';
    else
        fn2=sprintf(format,pict-1);
    end %if
    save(['fi_' fn2 '.mat'],'fz','ID','-mat'); % save filtered data
    contour(fz,CLineSet); title(['Filtered data. Frame ' fn2]); % display filtered data
    text(0.4,-0.085,['SaturatedPixels ' num2str(SatPix)],'units','normalized','color','r');
    if ~guidata(fHandle), return, end         % Stop button pressed.
    drawnow
end %for pict
%---------------------------------- subfunction Referenz
function Referenz(digits,DarkFramePresent,NumberOfFrames,CLineSet)
nullen='000';
load(['fi_' nullen(1:digits-1) num2str(DarkFramePresent) '.mat']); % load filtered data for 1st frame (after darkframe)
sze=size(fz); xmax=sze(2);ymax=sze(1);         % New array size. Filtering reduces array size by some pixels.
%----------------------------------------------------------------------------------
fHandle=figure(1);
contour (fz,CLineSet); title(['Reference area ' ID])
msgbox('Set points for reference area')
M=setMask(fz,fHandle);                  % Define reference area, M=2 innerhalb ROI
hold on
[dummy,ContourHandle]=contour(M,[1 1]);     % Display reference area by contour
set(ContourHandle,'LineColor','r')
hold off
set(fHandle,'Renderer','Painters')
msgbox('Processing Reference Area to yield reference Intensity')
IntMean=zeros(1,NumberOfFrames); % storage array for mean intensity in reference area
guidata(fHandle,1);
b = uicontrol('Style','PushButton',...
    'Parent',fHandle, ...
    'String','STOP', ...
    'Callback',['guidata(gcbo,0);'], ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','StopButton');
for pict=1:NumberOfFrames
    format=['%0' num2str(digits) 'u'];
    fn2=sprintf(format,pict-1);
    load(['fi_' fn2 '.mat']);     % load filtered data "fz" and "ID"
    z=(M>1).*fz;                  % only ref area
    STEP=(z>intlimit).*z;         % getting a 2D step function for z > intlimit
    [yInWindow,xInWindow]=find(STEP > 0); % find x and y coordinates where fz > intlimit
    IM=0; NInWindow=length(xInWindow);
    for k=1:NInWindow
        IM=IM+fz(yInWindow(k),xInWindow(k));
    end
    IntMean(pict)=IM/NInWindow;        % Mean Intensity
    contour(z,CLineSet);  % display filtered data again
    hold on                                   % and
    [dummy,ContourHandle]=contour(M,[1 1]);  title(['Reference area. Frame ' fn2]);   % show reference area
    set(ContourHandle,'LineColor','r')
    hold off
    if ~guidata(fHandle), return, end         % Stop button pressed.
    drawnow
end %for pict
save('MeanIntensityRefArea.mat','IntMean','-mat')


function [intlimit,CLineSet,M,Ax,sze]=FindIntLimROI(digits,DarkFramePresent,greenmap)
nullen='000'; screensize=get(0,'ScreenSize');
load(['fi_' nullen(1:digits-1) num2str(DarkFramePresent) '.mat']); % load filtered data for 1st frame
sze=size(fz); % filtering has changed array size
fzmax=max(max(fz));
fHandle=figure(3); colormap(greenmap); 
image(fz); set(gca,'Ydir','normal'); axis('equal'); set(fHandle,'Position',screensize); hold on
CLineSet=[0:6].^2*(fzmax*0.9)/6^2;intlimit=3;
contour(fz,CLineSet+intlimit,'k-')
guidata(fHandle,[1 intlimit CLineSet]);
hok=uicontrol('Style','PushButton',...
    'Parent',fHandle, ...
    'String','OK?', ...
    'Callback','data=guidata(gcbo);data(1)=0;guidata(gcbo,data);', ...
    'Units','normalized', ...
    'Position',[0 0 0.1 0.1], ...
    'Tag','OK?');
hHoeher=uicontrol('Parent',fHandle, ...
    'Units','normalized', ...
    'Callback',['data=guidata(gcbo);data(2)=data(2)+1;'....
                'fz=get(findobj(''type'',''image''),''Cdata'');delete(findobj(''type'',''hggroup''));'....
                'contour(fz,data(3:9)+data(2),''k-'');guidata(gcbo,data);'],...
    'Position',[0 0.6 0.1 0.1], ...
    'String','>', ...
    'Tag','Hoher');
hNiedriger=uicontrol('Parent',fHandle, ...
    'Units','normalized', ...
    'Callback',['data=guidata(gcbo);data(2)=data(2)-1;data(2)=max(data(2),1);'...
                'fz=get(findobj(''type'',''image''),''Cdata'');delete(findobj(''type'',''hggroup''));'....
                'contour(fz,data(3:9)+data(2),''k-'');guidata(gcbo,data);'],...
    'Position',[0 0.4 0.1 0.1], ...
    'String','<', ...
    'Tag','Niedriger');
title('Intensitätsschwelle vergrößern mit ">", verkleinern mit "<"')
data=guidata(fHandle); 
while data(1),data=guidata(fHandle);pause(0.2),  end     % OK-button pressed.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
delete(hok,hHoeher,hNiedriger);title(['ROI ', ID])
intlimit=data(2); CLineSet=data(3:9)+intlimit;

M=setMask(fz,fHandle);                  % Define ROI by mask M, M=2 innerhalb ROI
[dummy,ContourHandle]=contour(M,[1 1]);     % Display ROI durch Konturlinie
set(ContourHandle,'LineColor','r')
hold off %'''''''''''''''''''''''''''
axis equal, saveas(fHandle,'ROI.fig'); saveas(fHandle,'ROI.jpg')
Ax=axis;
msghndl=msgbox('Paused for printing; Next comes main anlysis');
set(msghndl,'position',[200 300 170 50]);
while ishandle(msghndl), pause(0.1), end
