clear;
clc;
close all;

p = 2000; % signal dimension
n = 600;  % number of measurements
noise = 0.05; % noise std. deviation

% sparse signal
thr = 0.95;  % sparsity level
x = randn(p,1).*(rand(p,1)>thr);
fprintf('... sparse signal created ... \n')

% measurements
Phi = randn(n,p)/sqrt(n);
y = Phi*x + noise*randn(n,1);
fprintf('... measurements obtained ... \n');

% CoGEnT parameters
tau = norm(x,1);
maxiter = 300;
selfun = @(gradf) find_next_atom_l1(gradf);
Ainit  = [1; zeros(p-1,1)];
gpiter = 10;
gptol = 1e-4;

% CoG
fprintf('... beginning CG ... \n')
[xhatcg, Atcg, objcg, timecg] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',0,...
    'gp_forward',0,...
    'maxiter',maxiter,...
    'stopcrit',1,...
    'gptol',gptol,...
    'gpiter',gpiter);


% CoGEn
fprintf('... beginning CoGEn ... \n')
[xhatcogen, Atcogen, objcogen, timecogen] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',0,...
    'gp_forward',1,...
    'maxiter',maxiter,...
    'stopcrit',1,...
    'gptol',gptol,...
    'gpiter',gpiter);


% CoGEnT
fprintf('... beginning CoGEnt ... \n')
[xhatcogent, Atcogent, objcogent, timecogent] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',1,...
    'gp_forward',1,...
    'maxiter',maxiter,...
    'stopcrit',1,...
    'gptol',gptol,...
    'gpiter',gpiter);



%FW classical
fprintf('... beginning classical FW ... \n')
[xhatfw, Atfw, objfw, timefw] = FW_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',0,...
    'gp_forward',0,...
    'maxiter',maxiter,...
    'stopcrit',1,...
    'gptol',gptol,...
    'gpiter',gpiter,...
    'classic',1);



% FW with full corrections
fprintf('... beginning FW with full corrections ... \n')
[xhatfull, Atfull, objfull, timefull] = FW_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',0,...
    'gp_forward',0,...
    'maxiter',maxiter,...
    'stopcrit',1,...
    'gptol',gptol,...
    'gpiter',gpiter,...
    'classic',0);


%% Print everything
clc
fprintf('... FW done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atfw,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatfw)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatfw,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timefw(end));

fprintf('... FW with corrections done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atfull,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatfull)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatfull,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timefull(end));

fprintf('... CG done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atcg,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatcg)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatcg,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timecg(end));

fprintf('... CoGEn done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atcogen,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatcogen)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatcogen,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timecogen(end));

fprintf('... CoGEnt done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atcogent,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatcogent)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatcogent,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timecogent(end));

%% plot stuff 

hold on
indstart = 3;
% plot(timefw(indstart:end),objfw(indstart:end),'c','LineWidth',2);
% plot(timefull(indstart:end),objfull(indstart:end),'k','LineWidth',2);
% plot(timecg(indstart:end),objcg(indstart:end),'b','LineWidth',2);
% plot(timecogen(indstart:end),objcogen(indstart:end),'m','LineWidth',2);
% plot(timecogent(indstart:end),objcogent(indstart:end),'r','LineWidth',2);
% legend('CG-fixed','CG-fc','CG','CG-En','CoGEnT')
% xlabel('Time (sec)');
% ylabel('Objective');


plot(log(objfw(indstart:end)),'c','LineWidth',2);
plot(log(objfull(indstart:end)),'k','LineWidth',2);
plot(log(objcg(indstart:end)),'b','LineWidth',2);
plot(log(objcogen(indstart:end)),'m','LineWidth',2);
plot(log(objcogent(indstart:end)),'r','LineWidth',2);
legend('CG-fixed','CG-fc','CG','CG-En','CoGEnT')
xlabel('Iterations');
ylabel('log(Objective)');

%% save 
save ROLE_OF_ENHANCE