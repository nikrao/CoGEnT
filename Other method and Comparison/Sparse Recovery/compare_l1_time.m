clear;
clc;
close all;
addpath(genpath('/home/nikhilrao/Dropbox/group-regularization/Code/FOBA_atomic'));

rng(0);
p = 20000;
n = 5000;
s = 1000;

Phi = randn(n,p)/sqrt(n);
for test = 1:5
inds = randsample(p,s);
xtrue = sparse(zeros(p,1));
xtrue(inds) = (2*rand(s,1)-1);
sig = 0.05;
y = Phi*xtrue + sig*randn(n,1);

fprintf('beginning timing comparisons \n');

%% cogent
fprintf('cogent \n');

gptol = 1e-3;
gpiter = 5;
maxiter = 1000;
tol = 1e-4;
Ainit = [1; zeros(p-1,1)];
Ainit = sparse(Ainit);
selfun = @(gradf) find_next_atom_l1(gradf);
tau = norm(xtrue,1);

[x_cogent, ~, obj_cogent, time_cogent,mcogent] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
            'backward',1,...
            'maxiter',maxiter,...
            'tol',tol,...
            'gpiter',gpiter,...
            'gptol',gptol);
        
        
        
%% cond grad
fprintf('cgradient \n')
 


[x_cgrad, ~, obj_cgrad, time_cgrad,mcgrad] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
            'backward',0,...
            'maxiter',maxiter,...
            'tol',tol,...
            'gpiter',gpiter,...
            'gptol',gptol);
        
        
%% cosamp

fprintf('cosamp \n')
[x_cosamp,d,obj_cosamp,time_cosamp]=cosamp(Phi,y,s,tol,maxiter);


%% subspace pursuit
fprintf('S Pursuit \n')
        
[Rec,obj_sp,time_sp] = CSRec_SP(s,Phi,y,maxiter,tol);
x_sp = Rec.x_hat;

%% grades
fprintf('grades \n');

[x_grades ,i,obj_grades,time_grades] = grades(y, Phi, s, tol, maxiter);

%% compute stats
Egrad(test) = norm(x_grades - xtrue)^2/norm(xtrue)^2;
Espur(test) = norm(x_sp - xtrue)^2/norm(xtrue)^2;
Ecosa(test) = norm(x_cosamp - xtrue)^2/norm(xtrue)^2;
Ecgra(test) = norm(x_cgrad - xtrue)^2/norm(xtrue)^2;
Ecoge(test) = norm(x_cogent - xtrue)^2/norm(xtrue)^2;

Tgrad(test) = max(time_grades);
Tspur(test) = max(time_sp);
Tcosa(test) = max(time_cosamp);
Tcgra(test) = max(time_cgrad);
Tcoge(test) = max(time_cogent);
end
%%
%save L1_TIME Egrad Tgrad Espur Tspur Ecgra Tcgra Ecoge Tcoge Ecosa Tcosa

fprintf('all done \n')

        
        