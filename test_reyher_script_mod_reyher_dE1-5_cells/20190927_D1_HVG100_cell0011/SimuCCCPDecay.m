function simudecay
sigma0=2.6;
D=0.8;
NumberOfFrames=110;  %
Noise=20 * 0.01 % Noise in Percent
L=47;     % Process length in µm
ActPos=20; % Distance from tip in µm
A=350;   % Signal amplitude
Noise=A*Noise; % Noise absolute
offset = 0;
% simulated data
N=500;           % Spatial channels
h=L/(N-1);       % Space stepwidth 
disp(['Spacestep (mikrons)  = ' num2str(h)])
x=(0:N-1)*h;     %"True" coordinate in mikrons (from 0 to L)
FrameTime=1:NumberOfFrames; % Time channels, Startframe at 1s after acitvation
Signal=zeros(N,NumberOfFrames);
for NFrame= 1:NumberOfFrames
    time=FrameTime(NFrame);
    Sig=gauss(x,ActPos,L,A,sigma0,time,D);
    Sig=Sig+(rand(size(Sig))-0.5)*Noise + offset; %add noise
    Signal(:,NFrame)= Sig;% Analytical, method of images
end
% set up filter function, Gaussian filter with radius rf
rf=2;
[xf,yf]= meshgrid(-rf:rf);
ff =-xf.*xf-yf.*yf;
ff =exp(ff/3/rf);
ff =ff/sum(sum(ff));
fSignal=filter2(ff,Signal,'valid');  % now smoothing as in CW-Analysis09
Signal(rf+1:N-rf,rf+1:NumberOfFrames-rf)=fSignal; % replace "valid" region

% plot simulated Color Coded Contour Plot  
figure(1), clf
[T,X]=meshgrid(-FrameTime,x);
[dummy,hc]=contourf(X,T,Signal,30);
%
set(hc,'LineColor','none');



% Simulated "Decay plots"
% at the activation spot according to Method of Images
ActPosN=round(ActPos/h);
n_simu=Signal(ActPosN,:);

fig4=figure(4), clf
plot(FrameTime,n_simu,'.')
hold on
text(0.5,0.95,'Start values ','units','normalized')
text(0.5,0.90,['Sigma0=' num2str(sigma0)],'units','normalized')
text(0.5,0.85,['D=' num2str(D) ' kBleach=' num2str(0)],'units','normalized')

% fitting n_simu(t)
D=0.5; A=120; offset=0;
Pstart=[A, D, offset]; 
HP=[sigma0, ActPos, L]; %constant parameters
% show start estimates
hL=plot(FrameTime,Dsimu(Pstart,FrameTime,HP),'g');
hT=title('Start curve. Paused for 5s.');
pause(2)
delete(hT); delete(hL);
% do the fit
Pout=fminsearch(@(P) Dchi(P,FrameTime,n_simu,HP),Pstart, optimset('Display','iter'));
D=Pout(2)
Pout
% show the result
hL = plot(FrameTime,Dsimu(Pout,FrameTime,HP),'g');
hT = title(['I_{center}, D_{fit}=', num2str(D)]);
axis([0 max(FrameTime)*1.05 0 max(n_simu)*1.05]);
hold off           %figure(4)


% subfunctions
function out=Dchi(P,x,y,HP)
S=Dsimu(P,x,HP);
out=sum(sum((y-S).^2));

function S=Dsimu(P,x,HP)
sigma0=HP(1); ActPos=HP(2); L=HP(3);
sigma=sqrt(2*P(2)*x + sigma0^2);
S = P(1)*sigma0 ./ sigma .*(1+exp(-4*ActPos^2 ./ (2*sigma.^2))-...     % left image
                                    exp(-4*(ActPos-L)^2 ./ (2*sigma.^2))) + P(3);

function G=gauss(x,ActPos,L,A,sigma0,time,D)
sigma = sqrt(2*D*time +   sigma0^2);
G= A*sigma0*1/sigma * exp(-(x-ActPos).^2 / (2*sigma^2));
G=G + A*sigma0*1/sigma*exp(-(x+ActPos).^2 / (2*sigma^2)); % left image
G=G - A*sigma0*1/sigma*exp(-(x-(2*L-ActPos)).^2 / (2*sigma^2)); % right image


