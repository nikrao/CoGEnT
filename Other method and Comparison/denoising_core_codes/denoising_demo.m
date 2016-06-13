%%%%%%
% This demo generates n equispaced samples of a signal composed of k
% complex sinusoids corrupted by Gaussian noise. Different denoising
% algorithms are compared.
%
% Please refer to
% B. Bhaskar, G. Tang, and B. Recht, "Atomic norm denoising with
% applications to line spectral estimation" for more details
% To run the codes, one needs sdpt3 (available at
% http://www.math.nus.edu.sg/~mattohkc/sdpt3.html)
% and cvx (available at http://cvxr.com/cvx/)
%
% Zoom in to see the frequencies.
%%%%%%
%%
clear;
clc;
javaaddpath java;
n=64; % # of equispaced samples
k=2;  % # of frequencies
SNR = 10; %SNR in dB


amps = exp(2*pi*1i*rand(k,1)); %unit amplitudes with random phase

%generate signal with random frequencies but separated by 2/n
[signal,true_c,true_poles] = moment_vector(n,k,'random',amps,2);
noise_std = norm(signal)/sqrt(n)*10^(-SNR/20);
observed = signal + noise_std*(randn(n,1) + 1i*randn(n,1))/sqrt(2);

e  = @(x) norm(signal(:)-x(:))/norm(signal); %inline function to evaluate mse
s  = @(x) svd(T(x));

figure(1);
plot(1:n,abs(signal),'b');grid on;xlim([0 n+1]);
figure(2);
stem(true_poles/(2*pi),abs(true_c),'bo');grid on;

%%
%AST-ADMM
T0 = clock;
[x_debias,~,fs,cs,~,~,~] = ast_denoise(observed, noise_std);
mseval = e(x_debias);
fprintf('AST-ADMM       MSE = %.4f, TiME = %.2f s\n',mseval,etime(clock,T0));
figure(1);hold on;plot(1:n,abs(x_debias),'r');hold off;
figure(2);hold on;stem(fs,ones(size(fs)),'rs');hold off;

%%
%AST-SDPT3
T0 = clock;
[x_debias,Tu,fs,cs,pts,polyval,~] = ast_denoise_sdpt3(observed, noise_std);
mseval = e(x_debias);
fprintf('AST-SDPT3      MSE = %.4f, TiME = %.2f s\n',mseval,etime(clock,T0));
figure(1);hold on;plot(1:n,abs(x_debias),'g');hold off;
figure(2);hold on;stem(fs,ones(size(fs)),'gx');hold off;

%%
%CADZOW
T0 = clock;
[x_debias,X,fs,cs,x] = cadzow_denoise(observed, k);
fs = fs/(2*pi);
mseval = e(x_debias);
fprintf('CADZOW         MSE = %.4f, TiME = %.2f s\n',mseval,etime(clock,T0));
figure(1);hold on;plot(1:n,abs(x_debias),'m');hold off;
figure(2);hold on;stem(fs,ones(size(fs)),'m*');hold off;

%%
%LASSO
grid_size = 2^12;
tau = sqrt(log(n)+log(4*pi*log(n)))*noise_std*1.5;
T0 = clock;
[x_lasso,c_lasso,~,~] = moment_sparsa(observed,tau*sqrt(n),grid_size);
idx = find(c_lasso);
fs_lasso = (idx-1)/grid_size;
est_basis = exp( 1i*(0:(n-1))'*2*pi*fs_lasso(:)' )/sqrt(n);
c_lasso_debias = est_basis\observed;
x_lasso_debias = est_basis*c_lasso_debias;
mseval = e(x_lasso_debias);
fprintf('LASSO          MSE = %.4f, TiME = %.2f s\n',mseval,etime(clock,T0));
figure(1);hold on;plot(1:n,abs(x_lasso_debias),'k');hold off;
figure(2);hold on;stem(fs_lasso,ones(size(fs_lasso)),'kd');hold off;


%%
%MUSIC
T0 = clock;
[x_debias,fs,cs] = music_denoise(observed, k);
mseval = e(x_debias);
fprintf('MUSIC          MSE = %.4f, TiME = %.2f s\n',mseval,etime(clock,T0));
figure(1);hold on;plot(1:n,abs(x_debias),'y');hold off;
figure(2);hold on;stem(fs,ones(size(fs)),'y+');hold off;

figure(1);legend('true','ast-admm','ast-sdpt3','cadzow','lasso','music');
figure(2);legend('true','ast-admm','ast-sdpt3','cadzow','lasso','music');
ylim([0 1.1]);




