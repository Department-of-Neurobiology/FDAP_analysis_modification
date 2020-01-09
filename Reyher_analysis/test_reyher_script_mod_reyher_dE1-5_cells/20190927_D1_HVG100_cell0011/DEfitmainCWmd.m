function DEfitmain
% Script for fitting decay plots
% 
clear all, close all
[fn,path]=uigetfile('*.fig','Select decay figure');
cd(path)
open(fn)
% Figure with "Decay Plot" is open


fighndl=gcf;
lh=findobj(fighndl,'Type','line');
t=get(lh,'xdata'); % Extraction of data from figure
t=t-t(1);
I=get(lh,'ydata');
figure
plot(t,I,'.');

% Start values:
% Start Amplitude, exp decay constant, Diffusion D, offset
Pin=[0.6 0.002 0.05 0]; 

S=desimu(Pin, t); % show start guess
hold on
plot(t,S,'g')
hold off
%return
pause

Pout=fminsearch(@(P) DEchi(P,t,I),Pin, optimset('Display','iter'));

figure(3)
S=desimu(Pout, t);
plot(t,I,'.');
hold on
plot(t,S,'r')
hold off
xlabel(Pout)
saveas(figure(3),'fit.fig')
saveas(figure(3),'fit.jpg')


Pout

% subfunctions
function out=DEchi(P,x,y)
S=desimu(P,x);
out=sum(sum((y-S).^2));


function S=desimu(P,x)
StartWidth=1/P(1);
S = exp(-P(2)*x) ./sqrt(2*P(3)*x + StartWidth^2) + P(4);

