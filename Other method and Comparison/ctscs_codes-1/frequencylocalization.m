
clear;
clc;
javaaddpath java;
n=32; % # of equispaced samples
k=3;  % # of frequencies
m=10;  % # of random samples

amps = exp(2*pi*1i*rand(k,1));%.*(randn(k,1).^2+.5); %random amplitudes with random phase

%generate signal with random frequencies but separated by 2/n
[signal,true_c,true_poles] = moment_vector(n,k,'random',amps,2); 

I = randperm(n);
I = sort(I(1:m));
Ic = setdiff(1:n,I);
xVals = signal(I); %get m random samples

fprintf('Running SDP...')
T0 = clock;
%recover missing samples and identify the frequencies from the random
%samples using semidefinite programming
[x_sdp,poles_sdp,coeffs_sdp] = ctscs_sdpt3(xVals,I,n); 
fprintf(' done. %.2f s\n',etime(clock,T0));

%%
cvx_begin
variable xh(n) complex;
variable u0;
variable u(n-1) complex;
variable T(n,n) complex hermitian;
variable t;
minimize ((trace(T)/n + t)/2)
subject to
T == toeplitz([u0;conj(u)],[u0;conj(u)]');
[T xh;
    xh' t] == semidefinite(n+1,n+1);
xh(I) == xVals;
cvx_end

%%
cvx_begin
variable q(n) complex;
variable H(n,n) hermitian;
maximize (real(signal'*q))
subject to
[H q;
    q' 1] == semidefinite(n+1,n+1);
trace(H) == 1;
for j = 1:n-1
    sum(diag(H,j)) == 0;
end
q(Ic) == 0;
cvx_end

%%
s = real(signal'*q);

%%
cvx_begin
variable q(n) complex;
variable H(n,n) hermitian;
minimize norm(q,2)
subject to
s == (real(signal'*q));
[H q;
    q' 1] == semidefinite(n+1,n+1);
trace(H) == 1;
for j = 1:n-1
    sum(diag(H,j)) == 0;
end
q(Ic) == 0;
cvx_end
%%
N = 2^13;
ff = (0:N-1)/N;
ff = ff(:);
nn = 0:n-1;
nn = nn(:);
A = exp(1i*2*pi*nn*ff');

Q = A'*q;


set(0,'DefaultTextFontName','Times','DefaultTextFontSize',36,...
    'DefaultAxesFontName','Times','DefaultAxesFontSize',36,...
    'DefaultLineLineWidth',3,'DefaultLineMarkerSize',6)

stem(true_poles/(2*pi),abs(true_c),'b--','LineWidth',3); 
hold on
plot(ff,abs(Q),'r')
hold off
%[legh,objh,outh,outm] = legend('true poles','dual polynomial');
title('Frequency Localization')
grid on;



