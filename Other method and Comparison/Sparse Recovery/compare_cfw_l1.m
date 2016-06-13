%%% compare different methods for L1 recovery in the presence of noise

clear;
clc
close all;
warning('off')


% parameters for all algos
rng(0);
n = 10000;  % signal dimension
k = 3000;   % number of measurements
s = 300;   % signal sparsity
ntests = 5; % number of tests to average over
noiserange = linspace(0,1,5); % range of noise values
tol = 1e-6;
maxiter = 1000;

% parmeters for cogent and fw
% A = [eye(n) -eye(n)];  % atomic set
Ainit = [1; zeros(n-1,1)];
selfun = @(gradf) find_next_atom_l1(gradf);
gp_forward = 0;
eta  =0.5;

nind = 0;
for noise = noiserange
    
    nind = 1+nind;
    fprintf('testing noise = %f \n',noise);
    tol = noise^2/sqrt(k);
    for test = 1:ntests
        
        x = 2*rand(n,1)-1;
        r = randsample(n,n-s);
        x(r) = 0; % this is the true signal
        
        Phi = randn(k,n)/sqrt(k);
        
        y = Phi*x + noise*randn(k,1); % observations
        
        fob_param = norm(x,1);
        
        % cosamp
        tic;
        [Sest,~] = cosamp(Phi,y,s,tol,maxiter,x);
        timec(test) = toc;
        err2c(test) = norm(Sest - x)^2/n;
        err1c(test) = sum((Sest~=0)~=(x~=0))/n;
        
        % cogent
        [xhat, ~, ~, TT] = CoGEnT_sparse(y, Phi, fob_param, Ainit, selfun,...
            'backward',1,...
            'maxiter',maxiter,...
            'tol',tol,...
            'gpiter',10);
        timef(test) = TT(end);
        err2f(test) = norm(xhat - x)^2/n;
        err1f(test) = sum((xhat~=0)~=(x~=0))/n;
        
        
        % fwolfe
        [xhat, ~, ~, TT] = CoGEnT_sparse(y, Phi, fob_param, Ainit, selfun,...
            'backward',0,...
            'maxiter',maxiter,...
            'tol',tol,...
            'gp_forward',0,...
            'gpiter',10);
        
        timew(test) = TT(end);
        err2w(test) = norm(xhat - x)^2/n;
        err1w(test) = sum((xhat~=0)~=(x~=0))/n;
        
        
        % subspace pursuit
        tic;
        Sest = CSRec_SP(s,Phi,y,maxiter);
        times(test) = toc;
        err2s(test) = norm(Sest.x_hat - x)^2/n;
        err1s(test) = sum((Sest.x_hat~=0)~=(x~=0))/n;
        
        
        % grades
        tic;
        [x i] = grades(y, Phi, s, tol, maxiter);
        timeg(test) = toc;
        err2g(test) = norm(xhat - x)^2/n;
        err1g(test) = sum((xhat~=0)~=(x~=0))/n;
        
        fprintf('.');
        
    end
    
    fprintf('\n this noise done \n');
    COGENTl2(nind) = mean(err2f);
    COGENTh(nind)  = mean(err1f);
    COSAl2(nind) = mean(err2c);
    COSAh(nind)  = mean(err1c);
    FWOLl2(nind) = mean(err2w);
    FWOLh(nind)  = mean(err1w);
    SPURl2(nind) = mean(err2s);
    SPURh(nind)  = mean(err1s);
    GRADESl2(nind) = mean(err2g);
    GRADESh(nind) = mean(err1g);
    
    TIMEcogent(nind) = mean(timef);
    TIMEfwolfe(nind) = mean(timew);
    TIMEspursu(nind) = mean(times);
    TIMEcosamp(nind) = mean(timec);
    TIMEgrades(nind) = mean(timeg);
    
    save l1_expt_time_noise...
        COGENTl2 COGENTh COSAl2 COSAh FWOLl2 FWOLh SPURl2 GRADESl2 GRADESh...
        SPURh noiserange TIMEcogent TIMEfwolfe TIMEspursu TIMEcosamp TIMEgrades...
    
end
    
        
    
        
        
        
