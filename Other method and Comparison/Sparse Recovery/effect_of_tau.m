clear;
clc;
close all;


rng(0);
p = 2000;
n = 500;
s = 100;

Phi = randn(n,p)/sqrt(n);
inds = randsample(p,s);
xtrue = sparse(zeros(p,1));
xtrue(inds) = (2*rand(s,1)-1);

y = Phi*xtrue + 0.05*randn(n,1);

fprintf('comparing to various tau \n');

gptol = 1e-3;
gpiter = 5;
maxiter = p;
tol = 1e-5;
Ainit = [1; zeros(p-1,1)];
Ainit = sparse(Ainit);
selfun = @(gradf) find_next_atom_l1(gradf);
tau = norm(xtrue,1);
tauset = tau*2.^(-2:4);

for t = 1:length(tauset)
	
	tau = tauset(t);

[x_cogent, ~, obj_cogent, time_cogent] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
            'backward',1,...
            'maxiter',maxiter,...
            'tol',tol,...
            'gpiter',gpiter,...
            'gptol',gptol);
   
OBJ{t} = obj_cogent;
TIM{t} = time_cogent;   
MSE(t) = norm(xtrue - x_cogent)/norm(xtrue)^2
        
	fprintf('t = %d of %d \n',t,length(tauset));
end