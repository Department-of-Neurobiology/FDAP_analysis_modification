function CW_Analysis
%Version 9 fused together with DecayFit03 and modified for automated analysis
%The buttons and unused steps are blocked by comments, the differences in the script: 
%1)		Always takes first file named "cell_112.tif" and changes working directory to that with the script file
%2)		Doesn't "Process or re-process reference area" by skipping the yes and cancel options
%3)		Automatically corrects for mean intensity
%4)		For DecayFit03 choosing the file was done for changing the working directory - is automatic for pwd now
%5)		Skips "Soma and tip in frame", directly assumes we only selected process in the area
%6)		The final parameter is saved in a atble in the upper directory
%Nataliya Trushina, 09.01.2020

% Version 9, addition / changes: 
% 1)    If dark frame is prsent, the intensity distribution on the main axis of the dark frame is subtracted 
%       from those of subsequent frames
% 2)    The mask M, describing the ROI, is stored in all.mat
% 3)    The dark frame also stored
close all
clear all
load('mymap')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% to be adjusted %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
DarkFramePresent=1; % 0 = not present, 1 = present
filen='cell_112.tif'; %is always called that when using bash renaming procedure
%[filen,cpath]=uigetfile('*.tif','Select the LAST tif-file of series.');
cd(pwd);
cpath=pwd;
drawnow;
[fn1, reminder]=strtok(filen,'_');fn1=[fn1 '_'];
NOP=strtok(reminder,'.');
digits=length(NOP)-1;
NumberOfFrames=str2num(NOP(2:digits+1))+1-DarkFramePresent; 
%fn1='CellA_';  % First part of the filename including "_". May not start with a number, may not have a "blank".
fn3='.tif';                       % Extension
%NumberOfFrames=238 ;
idx=strfind(cpath,'\');n=length(idx);
ID=cpath(idx(n-1)+1:idx(n)-1);                      % Identification. For figure captions
intlimit=3;                       % lowest intensity limit to be selected just above noise level AFTER FILTERING
colorindex=1;                     % Color used in tif-file, 1 = red, 2 = green, 3 = blue
rf = 2;                           % radius of Gaussian filter (larger values smooth more but blur details)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
CLineSet=[40 50 100 300 500 1000 2000];
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Filtering %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% set up filter function, Gaussian filter with radius rf
[xf,yf]= meshgrid(-rf:rf);
ff =-xf.*xf-yf.*yf;
ff =exp(ff/3/rf);
ff =ff/sum(sum(ff));
%----------------------------------------------------------------------------------
%button = questdlg('Apply or re-apply filtering? "Cancel" stops processing');
%if strcmp(button,'Cancel'), return, end
%if strcmp(button,'Yes')
    % filtering, subtraction of minimum (offset), mean intensity in ref-area
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
    %z=double(z); % raw data are now in double precision array "z"
    sze=size(z);
    xmax=sze(2); 
    ymax=sze(1);
    Z=z(ymax:-1:1,1:xmax,colorindex); % in viewer programs y=0 is top, in matlab  y=0 is bottom
    clear z; pack; % free memory from the very large array z
    IdxSat=find(Z>4094);
    SatPix=length(IdxSat);
    Z=double(Z); % raw data are now in double precision array "Z"
    %-----------------------------
    fz=filter2(ff,Z,'vaild');% now smoothing Z with f
    a=min(min(fz));fz=(fz-a); % now getting new minimum and subtract it
    if DarkFramePresent & pict==1
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
%end %if strcmp
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Reference Intensity %%%%%%%%%%%%%%%%%%%%%%%%%%%%
%button = questdlg('Process or re-process reference area? "Cancel" stops processing');
if ~exist(['fi_00' num2str(DarkFramePresent) '.mat'])
    msgbox('No filtered data available - stop processing')
    return
end %if ~exist

%{
if strcmp(button,'Cancel'), return, end
if strcmp(button,'Yes')
    print('We are here.')
    load(['fi_00' num2str(DarkFramePresent) '.mat']); % load filtered data for 1st frame (after darkframe)
    sze=size(fz);
    xmax=sze(2);ymax=sze(1);         % New array size. Filtering reduces array size by some pixels.
    %----------------------------------------------------------------------------------
    fHandle=figure(1);
    contour (fz,CLineSet); title(['Reference area ' ID])
    msgbox('Set points for reference area')
    M=setMask(fz,fHandle);                  % Define reference area
    hold on
    [dummy,ContourHandle]=contour(M,[1 1]);% Display reference area
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
        fn2=sprintf('%03u',pict-1);
        load(['fi_' fn2 '.mat']);     % load filtered data "fz" and "ID"
        z=(M>1).*fz;                  % only ref area
        STEP=(z>intlimit).*z;         % getting a 2D step function for z > intlimit
        [yInWindow,xInWindow]=find(STEP > 0);% find x and y coordinates where fz > intlimit
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
end %if strcmp
%}

if exist('MeanIntensityRefArea.mat')
    load('MeanIntensityRefArea.mat')
    figure(2)
    plot(IntMean), title(['Mean Intensity per pixel of reference area versus frame number -' ID])
    %msgbox('Paused for printing')
    %pause
end %if exist

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Main Analysis %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
load(['fi_00' num2str(DarkFramePresent) '.mat']); % load filtered data for 1st frame
sze=size(fz); ymax=sze(1); xmax=sze(2);% filtering has changed array size
fHandle=figure(3);
contour(fz,CLineSet); title(['ROI ', ID])
M=setMask(fz,fHandle);                  % Define ROI by mask M
hold on %''''''''''''''''''''''''''''
[dummy,ContourHandle]=contour(M,[1 1]); % Display ROI
set(ContourHandle,'LineColor','r')
hold off %'''''''''''''''''''''''''''
axis equal
saveas(fHandle,'ROI.fig')
saveas(fHandle,'ROI.jpg')
Ax=axis;
%msghndl=msgbox('Paused for printing; Next comes main anlysis');
%set(msghndl,'position',[200 300 170 50]);
%while ishandle(msghndl), pause(0.1), end


% Allocation
    CGcell=zeros(2,NumberOfFrames+DarkFramePresent);       % storage array for Center of Gravity coordinates
    IntTotalCell=zeros(1,NumberOfFrames+DarkFramePresent); % total intensity in the region of interest for Intensity > intlimit
    U=cell(NumberOfFrames+DarkFramePresent,3); % initialize cell array to hold projection data
    npmax=0;
% End Allocation

% rEFERENCE AREA business
    IntMean=ones(1,NumberOfFrames+DarkFramePresent);          % dummy if no reference area exists
    RefMeanI=1;                                % dummy if no reference area exists
    if exist('MeanIntensityRefArea.mat')
     load('MeanIntensityRefArea.mat');       % load array IntMean
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

%*******************************************************************
%*******************************************************************
%***********************    main loop    ***************************
%*******************************************************************
%*******************************************************************


for pict=1:NumberOfFrames+DarkFramePresent   % main loop
    fn2=sprintf('%03u',pict-1+DarkFramePresent);
    if pict == NumberOfFrames+1
        fn2='dark';
    end
    load(['fi_' fn2 '.mat']);     % load filtered data
    fz=fz*MeanRefI/IntMean(pict); % Correcting for constant intensity of the reference area
    %%%%%%%%%%%%%%%%%%%%%%%%% fz (= filtered z) containes the smoothed data with minimum = 0, array size is (1 to nsx, 1 to nsy) %%%%%%%%%%%%%%%%%%%%%%%%%
    z=(M>1).*(fz+eps);            % select ROI using  mask M; z is zero outside ROI and at least eps inside (eps is 2.2e-16 or so)
   
    if pict == 1
     newsize=size(z); newymax=newsize(1); newxmax=newsize(2); % überflüssig in Schleife - fz hat konstante Grösse
     [yInROI,xInROI]=find(z > 0);                       % find x and y coordinates where z > 0 
     NinROI=length(yInROI);                             % Number of pixels in ROI
    end %if pict
    
    STEP=(z>intlimit).*z;                              % getting a 2D step function for fz > intlimit
    [yInWindow,xInWindow]=find(STEP > 0);              % find x and y coordinates where z > intlimit
    % theses pixels are in the window treated in the following
    xCG=0;yCG=0;NinWindow=length(yInWindow);IntTotal=0;% Reset of summing variables
    % --------------------
    %
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
    figure(4)
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
        RMA=[v1(:)' ; v2(:)'];
        
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
           RROI(RMA,z,CLineSet) % Display rotated cell
    end %if
    
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
    %TEST
    %     for j=1:NinROI
    %         plot(xInROI(j),yInROI(j),'.y')
    %     end
    %     pause
    %END TEST
    
    hold off %''''''''''''''''''''''''''''''''''''''''''''''''' hold off
    % plot rotated ROI
    
    
    %-------------------------------------------------
    % Calculating projections on main axes
    if pict == 1                                % find binwidth and scale range for both axes for all frames
        P1=zeros(NinROI,1); P2=P1;              % all values are in "pixel units"
        for j=1:NinROI
            p1=v1.*[xInROI(j);yInROI(j)]; P1(j)=sum(p1); % projection of all pixel vectors in ROI on main axis 1
            p2=v2.*[xInROI(j);yInROI(j)]; P2(j)=sum(p2); % projection on main axis 2
        end
        low1=min(P1); high1=max(P1); low2=min(P2); high2=max(P2); 
        binwidth1=(high1-low1)/100;         % 100 bins (channels) along main axis
        binwidth2=(high2-low2)/100;         % fixed scale for all frames
    end %if pict
    IntPixel=zeros(NinWindow,1); P1=IntPixel; P2=P1;
    for j=1:NinWindow
        IntPixel(j)=z(yInWindow(j),xInWindow(j)); 
        p1=v1.*[xInWindow(j);yInWindow(j)]; P1(j)=sum(p1); % projection on main axis 1
        p2=v2.*[xInWindow(j);yInWindow(j)]; P2(j)=sum(p2); % projection on main axis 2
    end
    [P1bin,I1binAVG,I1binSUM]=put2bins(P1,IntPixel,low1,high1,binwidth1,0.0); % subfunction put2bins puts the actual value P1 into the appropriate bin
    [P2bin,I2binAVG,I2binSUM]=put2bins(P2,IntPixel,low2,high2,binwidth2,0.0); % P1bin(n)=low1+(n-1/2)*binwidth
    % form here on, axis 2 is not considered, so far
    np=length(P1bin);
    npmax=max([np npmax]);
    U{pict,1}=P1bin;               % hold the projetion data in cell array U
    U{pict,2}=I1binSUM;
    U{pict,3}=I1binAVG;
   %figure; plot(P1bin,I1binAVG);
   %pause(0.1) % to get the possibility of interruption by CTRL-C
    drawnow
    printing=1:50:400;             % print every 50th frame
    if sum(pict==printing)
        %disp('Paused for printing. Press any key to resume processing.')
        %pause
    end%if sum
    if ~guidata(fHandle), return, end         % Stop button pressed. 
end %for pict=... main loop
save('u.mat','U','NumberOfFrames','npmax','ff')
%**************************************************************************
%********************       END LOOP          *****************************
%**************************************************************************
%**************************************************************************


% ---- generating a 3d intensity surface z(x,y), with z= intensity, x=position along main axis, y= frame number
% allocating arrays
tI1binSum=zeros(NumberOfFrames,npmax); % Total intensity in bins
tI1binAvg=tI1binSum;   % Average intensity in bins
tP1binx=tI1binSum;     % "x"-coordinate along main axis of bins
tP1biny=tI1binSum;     % "y"-coordinate = frame number
% end allocation

%button = questdlg('Do you want to correct for constant mean intensity? "Cancel" stops processing');
%if strcmp(button,'Cancel'), return, end
%swtch=0;
%if strcmp(button,'Yes')
AvgTotalInt=sum(IntTotalCell)/NumberOfFrames;
swtch=1;
%end %if strcmp(button,'Yes')

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
                title(['Frame ' num2str(pict) '. Black=DarkFrame. Paused..']); %pause %## Diagnose ##
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

a=figure(8); % displays cross sections of Intensity(x,y) along x-axis
pict=1;
plot(tP1binx(pict,:),tI1binAvg(pict,:))
ax=axis;

save('All.mat','pict','NumberOfFrames','ax','tP1binx','tP1biny','tI1binAvg','tI1binSum','-mat')
BinCoordinatesPixels=tP1binx(1,:);
save('BinCoordinatesInPixels.txt','BinCoordinatesPixels','-ASCII')
save('totalInt_inBinsAvg.txt','tI1binAvg','-ASCII')
save('totalInt_inBinsSum.txt','tI1binSum','-ASCII')
b = uicontrol('Parent',a, ...
    'Units','points', ...
    'Callback',['if ~exist(''pict'');load(''All.mat'');end;pict=pict+1;pict=min(pict,NumberOfFrames);',...
    'plot(tP1binx(pict,:),tI1binAvg(pict,:));axis(ax);title(num2str(pict));'], ...
    'Position',[370 25 20 20], ...
    'String','>', ...
    'Tag','Pushbutton1');
b = uicontrol('Parent',a, ...
    'Units','points', ...
    'Callback',['if ~exist(''pict'');load(''All.mat'');end;pict=pict-1;',...
    'pict=max(pict,1);plot(tP1binx(pict,:),tI1binAvg(pict,:));axis(ax);title(num2str(pict));'], ...
    'Position',[330 25 20 20], ...
    'String','<', ...
    'Tag','Pushbutton2');

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
%Save final parameter into a table
disp(num2str(D))
dlmwrite('..\Deff_auto.csv',num2str(D),'delimiter','','-append')
[~, ParentFolderName] = fileparts(pwd);
outvar = [ParentFolderName, ';', num2str(D)];
dlmwrite('..\name_Deff_auto.csv',outvar,'delimiter','','-append')
%Check the fit
msgbox('Paused for checking the fit. Press any key to continue.')
disp('Press any key to close session.')
pause
clear all
close all

%% --------------------------------------------------------------------------
% Subfunctions


%pause
%
%
%--------------------------------------------------------- subfunctions ------------------------------------------------
%
%
function M=setMask(fz,fHandle)
sze=size(fz);
xmax=sze(2);ymax=sze(1);
figure(fHandle)
text(0.4,-0.085,'Define ROI!','units','normalized','color','r');
hold on
% Get 4 points defining mask area
rx = zeros(4,1);
ry = rx;
n = 0;
% Loop, picking up the points.
%disp('Left mouse button picks 4 points.')
while 1
    [xi,yi,button] = ginput(1);
    plot(xi,yi,'ro')
    n = n+1;
    rx(n)=xi; ry(n)=yi;
    if n==4, break, end
end
% center of mask area
X=sum(rx)/4; Y=sum(ry)/4;
% corner points with respect to center
for k=1:4
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
dx(1:3)=ux(2:4)-ux(1:3);
dx(4)=ux(1)-ux(4);
dy(1:3)=uy(2:4)-uy(1:3);
dy(4)=uy(1)-uy(4);

% Define mask M
M=zeros(ymax,xmax);
% Loop to set the border line ins the mask
step=1e-3;
for k=1:4
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
hold off
%
%
%--------------------------------------------------------- subfunctions ------------------------------------------------
%
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
nB=floor((B-Bmin)/binWidth) +1;       % nB ordnet jedem B-Wert eine bin- (oder Intervall-) Nummer zu
nBstop=floor((high-Bmin)/binWidth) +1;;  % Die Intervalle laufen von 1 bis nBstop. Nicht jedes Intervall kommt vor,
% d.h. nB enthält i.A. nicht jede Zahl von
% 1 bis nBstop
kB=1:nBstop;                    % Index aller Bins

AbinAVG=kB*0;
AbinSUM=AbinAVG;
for k=kB                           % Schleife vektorisieren = ?
    binmembers=find(nB==k);
    AbinSUM(k)=sum(A(binmembers));   % Summe von Eregignissen
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


%---------------------------------------

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
contour(XP,YP,z, CLineSet-2)
saveas(fHandle,'RotCell.fig')
close(fHandle);

%-----------------------------------------


%DecayFit03.m subfunctions below

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

