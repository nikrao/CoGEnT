% script to test FB, FW, COSAMP and SP on multiple instances

% problem parameters
clear;clc; close all;
p = 1024;
n = 200;
sigma = 0.001;
ntests = 1000;


% error initializations
FBerr = zeros(ntests,1);
FWerr = zeros(ntests,1);
SPerr = zeros(ntests,1);
COSAMPerr = zeros(ntests,1);
FBsparsity = zeros(ntests,1);
FWsparsity = zeros(ntests,1);
SPsparsity = zeros(ntests,1);
COSAMPsparsity = zeros(ntests,1);
timesFB = zeros(ntests,1);
timesFW = zeros(ntests,1);
maxiter = 1000;
sparsity = 64;
tol = 1e-20;

% specifics for various algos
A = [eye(p) -eye(p)];  % atomic set
Ainit = A(:,end);
selfun = @(gradf) find_next_atom_l1(gradf,A);
gptol = 1e-4;
gpiter = 100;
for test = 1:ntests
    
    % generate sparse signal
    x = zeros(p,1);
    locs = randsample(p,sparsity);
    x(locs) = 2*rand(sparsity,1) - 1;
    
    % generate measurements
    Phi = randn(n,p)/sqrt(n);
    y = Phi*x + sigma*randn(n,1);
    
    % FOBA
    tau = norm(x,1)*1.01;
    [xfoba, At,iter, obj, time, backs] = FOBA_L2(y, Phi, tau, Ainit, selfun, ...
    'mbackward', 1, ...
    'maxiter',maxiter,...
    'eta', 0.5,...
    'gptol', gptol,...
    'gpiter', gpiter,...
    'verbose', 0);

    FBerr(test) = norm(x-xfoba)^2/p;
    timesFB(test) = time(end);
    FBsparsity(test) = size(At,2);
    
    % COSAMP
    [xcosamp,errs] = cosamp(Phi,y,sparsity,tol,maxiter,x);
    
    COSAMPerr(test) = norm(xcosamp-x)^2/p;
    COSAMPsparsity(test) = nnz(xcosamp);
    
    % SP
    Rec = CSRec_SP(sparsity,Phi,y,maxiter);
    xsp = Rec.x_hat;
    SPerr(test) = norm(x-xsp)^2/p;
    SPsparsity(test) = nnz(xsp);
    
    % FrankWolfe
    [xfw, At,iter, obj, time, backs] = FOBA_L2(y, Phi, tau, Ainit, selfun, ...
    'backward',0,...    
    'mbackward', 1, ...
    'maxiter',maxiter,...
    'eta', 0.5,...
    'gptol', gptol,...
    'gpiter', 1,...
    'verbose', 0);

    FWerr(test) = norm(x-xfw)^2/p;
    timesFW(test) = time(end);
    FWsparsity(test) = size(At,2);
    
    fprintf('test %d done \n',test);
    
end
save COMPARE_GREEDYS FWerr FBerr COSAMPerr SPerr timesFW timesFB...
    FWsparsity FBsparsity SPsparsity COSAMPsparsity

% make performance profiles MSE
close all
nbetas = 1000;
betas = logspace(0,5,nbetas);
E = [FBerr FWerr COSAMPerr SPerr];% error matrix

M = min(E')';
Mr = repmat(M,1,nbetas);
br = repmat(betas,ntests,1);
B = br.*Mr;  % each element contains \beta*min(that test)

% PP for FOBA
V = FBerr;
V = repmat(V,1,nbetas);
V = V<=B;
FOBAPP = sum(V)/ntests;

% PP for FrankWolfe
V = FWerr;
V = repmat(V,1,nbetas);
V = V<=B;
FWPP = sum(V)/ntests;

% PP for COSAMP
V = COSAMPerr;
V = repmat(V,1,nbetas);
V = V<=B;
COSAMPPP = sum(V)/ntests;

% PP for SP
V = SPerr;
V = repmat(V,1,nbetas);
V = V<=B;
SPPP = sum(V)/ntests;

% make performance profiles SPARSITY
E = [FBsparsity FWsparsity COSAMPsparsity SPsparsity];% error matrix

M = min(E')';
Mr = repmat(M,1,nbetas);
br = repmat(betas,ntests,1);
B = br.*Mr;  % each element contains \beta*min(that test)

% PP for FOBA
V = FBerr;
V = repmat(V,1,nbetas);
V = V<=B;
FOBAPPs = sum(V)/ntests;

% PP for FrankWolfe
V = FWerr;
V = repmat(V,1,nbetas);
V = V<=B;
FWPPs = sum(V)/ntests;

% PP for COSAMP
V = COSAMPerr;
V = repmat(V,1,nbetas);
V = V<=B;
COSAMPPPs = sum(V)/ntests;

% PP for SP
V = SPerr;
V = repmat(V,1,nbetas);
V = V<=B;
SPPPs = sum(V)/ntests;

save COMPARE_GREEDYS FWerr FBerr COSAMPerr SPerr timesFW timesFB...
    SPPP COSAMPPP FWPP FOBAPP SPPPs COSAMPPPs FWPPs FOBAPPs sigma
% plot(betas,SPP,'k')
hold on
plot(betas,COSAMPPPs,'g')
plot(betas,FOBAPPs)
plot(betas,FWPPs,'r')
legend('COSAMP','FB','FW');