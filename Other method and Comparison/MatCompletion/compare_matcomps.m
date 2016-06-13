% cogent matcomp comparisons 

clear;
clc;
close all;
warning('off');
%% create a low rank matrix
maxiter = 1000;
tol =1e-6;
N =  [60
         100
         150
         250
         500
         750
        1000
        1500
        2000]
for index = 1:7
n1 = N(index)
n2 = ceil(n1*4/3);
r = max(ceil(min(n1,n2)/100),2);
X0 = randn(n1,r)*randn(r,n2);
sig = 0.05;
noise = sig*randn(size(X0));
X_noise = X0 + noise;

%% create the "sampling" matrix and obtain measurements
sam_fraction = 0.25;
M = rand(n1,n2);
M = M<=sam_fraction;
J = M(:);
Phi = spalloc(n1*n2,n1*n2,n1*n2);
for ii = 1:n1*n2
   Phi(ii,ii)= J(ii);
end
clear J;
% Phi = sparse(diag(M(:)));
y = Phi*X_noise(:);

%% cogent
tau = sum(svd(X0));
selfun = @(gradf) find_next_atom_nucnorm(gradf,n1,n2);
Ainit = randn(n1,1)*randn(1,n2);
Ainit = Ainit(:)/norm(Ainit(:));
[x, Atc, iterc, objc, timec, back_count] = CoGEnT_MC(y, Phi, tau, Ainit, [n1,n2] ,selfun,...
    'maxiter',maxiter,'tol',tol,'debias',1,'gptol',1e-3,'gpiter',5,'gp_forward',0);

X_cogent = reshape(x,n1,n2);
err_cogent(index) = norm(X_cogent - X0,'fro')^2/norm(X0,'fro')^2;
time_cogent(index) = timec(iterc);
rank_cogent(index) = rank(X_cogent);
fprintf('\n cogent done \n')

%% cg
[x, At, iter, obj, time, back_count] = CoGEnT_MC(y, Phi, tau, Ainit, [n1,n2] ,selfun,...
    'maxiter',maxiter,...
    'backward',0,'tol',tol);

X_cg = reshape(x,n1,n2);
err_cg(index) = norm(X_cg - X0,'fro')^2/norm(X0,'fro')^2;
time_cg(index) = time(iter);
rank_cg(index) = rank(X_cg);
fprintf('\n cg done \n')

%% OptSpace
D = X_noise.*M;
tic;
[X S Y dist] = OptSpace(D,[],maxiter,tol);
X_optspace = X*S*Y';
t = toc;
err_optspace(index) = norm(X_optspace - X0,'fro')^2/norm(X0,'fro')^2;
time_optspace(index) = t;
fprintf('\n OptSpace done \n')

%% ALM
D = X_noise.*M ;
tic
[A iter svp] = inexact_alm_mc(D,tol,maxiter);
t = toc;
time_alm(index)=  t;

X_alm = A.U*A.V';
err_alm(index) = norm(X_alm - X0,'fro')^2/norm(X0,'fro')^2;

fprintf('\n ALM done \n')

%% SET
D = X_noise.*M ;
tic
[U,S,V,Xr_norm] = MatrixSET02(D,M,r,tol);
t = toc;
time_set(index) = t;
X_set = U*S*V';
err_set(index) = norm(X_set - X0,'fro')^2/norm(X0,'fro')^2;

fprintf('\n SET done')

save matcomp_times time_set time_alm time_optspace time_cogent time_cg...
    err_cogent err_cg err_optspace err_alm err_set
end

