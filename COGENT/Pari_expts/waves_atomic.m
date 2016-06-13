clear
clc

%Square wave
%Triangle wave
%Ricker wavelet

N=100;
t=linspace(0,10,N);

x1=.4*sign(sin(.1*t));
x2=.3*sign(sin(.2*t-pi/3));
x3=.6*sign(sin(.3*t-.7*pi));

y=x1+x2+x3;
y=y';


% COGENT Formulation
Phi=eye(N,N);
tau=1.3;
Ainit=sign(sin(t-.2*pi))';
%selfun=@(gradf) selfun_square;
selfun = @(gradf) selfun_square(gradf);
[x, At, iter, obj, time, back_count]=CoGEnT(y, Phi, tau, Ainit, selfun,'tol',.01,'verbose',1,'dropcount',100) 

hold on
plot(y,'r')