% OSCAR using CoGEnT

clear;
clc;
close all;

% Measurements and true signal
n = 100;
p = 5000;
Phi = randn(n,p);
beta = zeros(p,1);
a = 7*ones(20,1)+ randn(20,1);
b = 9*ones(20,1)+ randn(20,1);
% c = -8*ones(20,1)+ randn(20,1);
beta(1:20) = a; beta(301:320)=b; 
noise = 0.1;
y = Phi*beta + noise*randn(n,1);

% cogent parameters
maxiter = 500;
eta = 0.5;
gptol = 1e-3;
gpiter = 50;
gp_forward = 1;


fprintf('... Measurements Obtained. Begin CoGEnT ... \n')

Ainit = [1; zeros(p-1,1)];
tauset = 250:25:450;
cset = 2.^(-10:2:10);

for tind = 1:length(tauset)
    tau = tauset(tind);
    for cind = 1:length(cset);
        c = cset(cind);
        
        
        selfun =@(gradf) find_next_atom_oscar(gradf,c);
        [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'gp_forward',gp_forward,...
            'stopcrit',1,...
            'debias',1);
        
        fprintf('.');
        
        errt(tind,cind) = norm(x-beta)^2/p;
        
    end
    fprintf('\n');
end

fprintf('done \n')


[t c]  = find(errt == min(min(errt)));

tau = tauset(t(1));
c = cset(c(1));
selfun =@(gradf) find_next_atom_oscar(gradf,c);
        [x, At, objcog, timecog] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'gp_forward',gp_forward,...
            'stopcrit',1,...
            'debias',1);
 
 stem(beta);
 hold on
 plot(x,'r .');
 legend('Original','Recovered')
 fprintf('\n MSE = %f \n',norm(x-beta)^2/p);
 fprintf('\n Mean Hamming Error = %f \n',sum((x~=0).*(beta~=0))/p);
x_cogent = x;
err_cogent = norm(x_cogent- beta)^2/p

%% frank wolfe


Ainit = [1; zeros(p-1,1)];

for tind = 1:length(tauset)
    tau = tauset(tind);
    for cind = 1:length(cset)
        c = cset(cind);
        
        selfun =@(gradf) find_next_atom_oscar(gradf,c);
        [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'gp_forward',0,...
            'stopcrit',1,...
            'debias',1,...
            'backward',0);
        
        fprintf('.');
        err(tind,cind) = norm(x - beta)^2/p;
        
    end
    fprintf('\n');
end

fprintf('rerunning cond gradient with best params \n');
[t c] = find(err == min(min(err)));
t = t(1); c = c(1);
tau = tauset(t);
c = cset(c);
selfun =@(gradf) find_next_atom_oscar(gradf,c);
[x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun, ...
    'maxiter',maxiter,...
    'eta', eta,...
    'gptol', gptol,...
    'gpiter', gpiter,...
    'gp_forward',0,...
    'stopcrit',1,...
    'debias',1,...
    'backward',0);

time_cgrad = time(end);
err_cgrad =  norm(x - beta)^2/p;
x_cgrad = x;

fprintf('\n cond gradient done \n')



%% adlas


lam_range = 2.^(-20:-6);

options.verbosity = 0;
options.iterations = maxiter;

for lind = 1:length(lam_range);
    lambda = lam_range(lind);
    
    for cind= 1:length(cset)
        c = cset(cind);
        
%         lambda_vec = lambda*(ones(p,1) + c*(p*ones(p,1) - (1:p)'));

        % DWSL1
        [x_hat,nIter,obj, timeSteps, errorSteps] = SolveFISTA_oscar(Phi,y,c,...
            'maxiteration',maxiter,'lambda',lambda);
        
        fprintf('.');
        nnz(x_hat)
        err(lind,cind) = norm(x_hat - beta)^2;

    end
    fprintf('\n');
end

fprintf('rerunning fista with best params \n');
[t c] = find(err == min(min(err)));
t = t(1); c = c(1);
lambda = lam_range(t);
c = cset(c);

% lambda_vec = lambda*(ones(p,1) + c*(p*ones(p,1) - (1:p)'));

% DWSL1
[x_hat,nIter,obj, timeSteps, errorSteps] = SolveFISTA_oscar(Phi,y,c,...
            'maxiteration',maxiter,'lambda',lambda);


err_prox =  norm(x_hat - beta)^2/p;
x_prox = x_hat;

fprintf('\n fista done \n')

