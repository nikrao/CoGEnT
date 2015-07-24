% TestBench to run a standard L1 recovery using the CoGEnT method

clear;
clc;
close all;
warning off

p = 2000; % signal dimension
n = 500;  % number of measurements
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
maxiter = 1000;
selfun = @(gradf) find_next_atom_l1(gradf);
Ainit  = [1; zeros(p-1,1)];
fprintf('... beginning CoGEnT ... \n')

% run CoGEnT
fprintf('********************************************************* \n')
fprintf('RUNNING CoGEnT \n')
fprintf('********************************************************* \n')
[xhat, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',1,...
    'maxiter',maxiter,...
    'gp_forward',1,...
    'stopcrit',1);
fprintf('... CoGEnT done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(At,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhat)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhat,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',time(end));


% compare to conditional gradient
fprintf('********************************************************* \n')
fprintf('COMPARING TO CONDITIONAL GRADIENT \n')
fprintf('********************************************************* \n')

[xhatcg, Atcg, objcg, timecg] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
    'backward',0,...
    'gp_forward',0,...
    'maxiter',maxiter,...
    'stopcrit',1);
fprintf('... CG done ... \n');
fprintf('ORIGINAL SPARSITY = %d \n',nnz(x));
fprintf('FINAL SPARSITY    = %d \n',size(Atcg,2));
fprintf('NORMALIZED MSE    = %f \n',norm(x-xhatcg)^2/norm(x)^2);
fprintf('AVG. L1 ERROR     = %f \n',norm(x-xhatcg,1)/p);
fprintf('TIME TAKEN        = %f seconds \n',timecg(end));



% plot results
figure
subplot(211)
stem(x);hold on; plot(xhat,'r .');
title('CoGEnT')
subplot(212)
stem(x);hold on; plot(xhatcg,'r .');
title('CG')

figure
plot(time(2:end),log(obj(2:end)),'r','LineWidth',2);
hold on
plot(timecg(2:end),log(objcg(2:end)),'b','LineWidth',2);
legend('CoGEnT','CG')
xlabel('Time (seconds)','FontSize',18)
ylabel('log(f(x))','FontSize',18)
grid on

figure
plot(log(obj(2:end)),'r','LineWidth',2);
hold on
plot(log(objcg(2:end)),'b','LineWidth',2);
legend('CoGEnT','CG')
xlabel('Iterations','FontSize',18)
ylabel('log(f(x))','FontSize',18)
grid on