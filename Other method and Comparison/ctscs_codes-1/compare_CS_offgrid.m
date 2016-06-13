% compare CS_offgrid methods

clear;
clc;
close all;
lengths = 2^11;
ntests = 10;
for ind = 1:ntests
%% generate signal
N =lengths;  % signal length
n =ceil(N/4);   % number of random measurements obtained
sig = 0.01;
S=randsample(N,n);
S=S-1;
act_freq = rand(5,1);
a0=rand(length(act_freq),1);     % random weights for the actual frequencies
Atrue = [];
% true  signal
x=zeros(n,1);
for jj = 1:length(act_freq)
    
    x = x + 1/sqrt(n)*a0(jj)*exp(1i*2*pi*act_freq(jj)*S);
    Atrue = [Atrue exp(1i*2*pi*act_freq(jj)*S)];
    
end
xtrue = x; % noiseless signal
x = x+sig*randn(n,1);  % observations


%% CoGEnT
tau = norm(a0,1) + sig;   % regularization parameter (clairvoyant)
epsilon = 1e-3;
maxiter = 2*n;
tic
[poles_cogent, coeffs_cogent] = CoGEnT_offgrid(x,S,maxiter,epsilon, tau);
time_cogent(ind) = toc;
L = length(poles_cogent);
x_cogent = 0;
for ii = 1:L
    x_cogent = x_cogent + 1/sqrt(n)*coeffs_cogent(ii)*exp(1i*2*pi*poles_cogent(ii)*([S])); 
end
% err_cogent = norm(x_cogent - xtrue)^2/N;
fprintf('cogent done \n')

%% DAST [Tang et al'13]
dastfreq = act_freq*2*pi;
signal = sum( diag(a0)*exp(1i * dastfreq * (0:N-1)),1).';
xObs = signal(S+1);
% tic;
% [x_dast,poles_dast,coeffs_dast] = ctscs_dast(xObs,S'+1,N,64); 
% time_dast(ind) = toc;
% err_dast = norm(xtrue-x_dast(S))^2/N;
fprintf('\n dast done \n')

%% SDP implementation
% tic;
% [x_sdp,poles_sdp,coeffs_sdp] = ctscs_sdpt3(xObs,S'+1,N);
% time_sdp(ind) = toc;
% % err_sdp = norm(xtrue-x_sdp(S+1))^2/N;
% fprintf('\n sdp done \n')


end

mean(time_cogent)
% mean(time_sdp)
std(time_cogent)
% std(time_sdp)