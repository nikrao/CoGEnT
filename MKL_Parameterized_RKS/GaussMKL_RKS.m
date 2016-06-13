function [w,Xout] = GaussMKL_RKS(y,X,SIGS,D, algparams)

% function that starts with an initial set of kernels, and refines to give a "good" representaiton

% D = number of random features to draw

[n,p] = size(X);

Z = zeros(n,D*length(SIGS));

% initial setup
groups = spalloc(D*length(SIGS),length(SIGS),D*length(SIGS));


for s = 1:length(SIGS)
% generate a RKS matrix for each std. dev
stddev = SIGS(s);
w = randn(p, D)*stddev/sqrt(p);
inds = (s-1)*D +1 : s*D;
Z(:,inds) = exp(i*X*w);
groups(inds,s) = 1;
end
	
% ALGO PARAMETERS
T = algparams.T;
tau = algparams.tau;
gptol = 1e-3;
gpiter= 5;
Ainit = zeros(D*length(SIGS),1);
Ainit(1:D) = 1/sqrt(D);
maxiter = 20;

siglist = SIGS;


%%%%%%%%%%%%%%%
for t = 1:T

% cogent part
selfun =@(gradf) find_next_atom_group(gradf,groups);
[x, At, ~, ~,~] = CoGEnT_sparse(y, Z, tau, Ainit, selfun , 'maxiter',maxiter)


% determine what variances have been selected, using At
for a = 1:size(At,2)
atom = At(:,a);
ind = find(atom~=0);
ind = ind(1);
% ind(1) = (s-1)*D+1
ind = ((ind-1)/D) + 1;
selsigs(a) = siglist(ind);
end
	

% add new variances
[signew] = addsigs(selsigs);
siglist = [siglist signew];

% use new variances to remake the Z, and groups
l = size(groups,2);
for s = 1:length(signew)
% generate a RKS matrix for each std. dev
stddev = signew(s);
w = randn(p, D)*stddev/sqrt(p);
Z = [Z exp(i*X*w)];
groups(inds,l+s) = 1;
end


end

