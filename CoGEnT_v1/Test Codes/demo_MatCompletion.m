% test CoGEnT for matrix completion

clear ;
clc;
close all;
warning off;

n1 = 100; n2 =80; r = 2; % matrix is n1Xn2 rank = r;
sam_frac = 0.35; %sampling fraction
Xtrue = randn(n1,r)*randn(r,n2);

xtrue = Xtrue(:); %vectorize the matrix

% create the matrix mask
M = rand(n1,n2);
M = M<=sam_frac;
% Y = M.*X;
Phi = sparse(diag(M(:))); % this is the mask we will use

% noise parameter
noisevar = 0.01;

% measurements
y = Phi*xtrue + noisevar*randn(n1*n2,1);
Y = M.*Xtrue;

% CoGEnT PARAMETERS
selfun =@(gradf) find_next_atom_nucnorm(gradf,n1,n2);
uinit = randn(n1,1);vinit = randn(1,n2);
uinit = uinit/norm(uinit);
vinit = vinit/norm(vinit);
Ainit = uinit*vinit;
Ainit = Ainit(:);
maxiter = 100;
tol = 1e-10;
eta = 0.5;
gptol = 1e-5;
gpiter = 50;
tau = sum(abs(svd(Xtrue)));
svt = 10;


% CoGEnT with SVT in the Backward Step
tic
[x, ~, ~, ~,~,~] = CoGEnT_MC(y, Phi, tau, Ainit, [n1,n2],selfun,...
    'maxiter',maxiter,...
    'eta', eta,...
    'gptol', gptol,...
    'gpiter', gpiter,...
    'gp_forward',1,...
    'backward',1,...
    'verbose', 0, ...
    'tol',tol,...
    'debias',1,...
    'svt',svt);
toc
errs = norm(xtrue-x)^2/numel(x);
recs = x;


% CoGEnT without SVT in the Backward Step
tic
[x, ~, ~, ~,~,~] = CoGEnT_MC(y, Phi, tau, Ainit, [n1,n2],selfun,...
    'maxiter',maxiter,...
    'eta', eta,...
    'gptol', gptol,...
    'gpiter', gpiter,...
    'gp_forward',1,...
    'backward',1,...
    'verbose', 0, ...
    'tol',tol,...
    'debias',1);

toc
errc = norm(xtrue-x)^2/numel(x);
recc = x;

Xc = reshape(recc,n1,n2);
Xs = reshape(recs,n1,n2);

figure
stem(svd(Xtrue),'MarkerSize',12)
hold on
plot(svd(Xs),'g o','MarkerSize',12);
plot(svd(Xc),'k x','MarkerSize',12);
xlim([0.8,3*r])
legend('True Singular Values','CoGEnT-SVT','CoGEnT-noSVT')