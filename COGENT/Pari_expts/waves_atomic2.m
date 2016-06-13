clear
clc

%Square wave
%Triangle wave
%Ricker wavelet

N=100;
t=linspace(0,10,N);

% x1=.4*sign(sin(.1*t));
% x2=.3*sign(sin(.2*t-pi/3));
% x3=.6*sign(sin(.3*t-.7*pi));

x1=.4*sawtooth(.8*t)/norm(sawtooth(.8*t));
x2=.3*sawtooth(.7*t-pi/3)/norm(sawtooth(.7*t-pi/3));
x3=.6*sawtooth(.3*t-.7*pi)/norm(sawtooth(.3*t-.7*pi));

y=x1+x2+x3;
y=y';

plot(y)
%figure
%pause

% COGENT Formulation
Phi=eye(N,N);
tau=1.3;
Ainit=sawtooth(.4*t-.3*pi)'/norm(sawtooth(.4*t-.3*pi));
%selfun=@(gradf) selfun_square;
selfun = @(gradf) selfun_sawtooth(gradf);
[x, At, ~]=CoGEnT(y, Phi, tau, Ainit, selfun,'tol',.001,'verbose',1); 

hold on; plot(y,'r .')