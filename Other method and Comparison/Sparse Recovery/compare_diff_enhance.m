clear;
clc;
close all;

p = 500; % signal dimension
n = 200;  % number of measurements
noise = 0.05; % noise std. deviation

% sparse signal
thr = 0.9;  % sparsity level
iterset = [1:20 50:50:250];

for run = 1:5

x = randn(p,1).*(rand(p,1)>thr);
fprintf('... sparse signal created ... \n')

% measurements
Phi = rand(n,p);
Phi = sparse(Phi>0.5);
y = Phi*x + noise*randn(n,1);
fprintf('... measurements obtained ... \n');

% CoGEnT parameters
tau = norm(x,1);
maxiter = 5000;
selfun = @(gradf) find_next_atom_l1(gradf);
Ainit  = [1; zeros(p-1,1)];

gptol = 1e-5;

iterind = 0;
for iter = iterset
    iterind = 1+iterind;
    % CoGEnT
    [w, Atcogent, obj, t] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
        'backward',0,...
        'gp_forward',1,...
        'maxiter',maxiter,...
        'stopcrit',0,...
        'gptol',gptol,...
        'gpiter',iter,...
        'tol',1e-5,...
        'true_x',x);
    
    nmse(run,iterind) = norm(x-w)^2/norm(x)^2;
    T(run,iterind) = t(end);
    S(run,iterind) = nnz(w);
   
    fprintf('.')
end

fprintf('\n')

end
% save ENHANCEMENT S T nmse iterset