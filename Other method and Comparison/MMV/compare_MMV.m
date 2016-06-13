% script for MMV comparison

clear;clc;close all

m = 300;
n = 1000;
L = 5;
r = 50;
act_inds = randsample(n,r);
SNR = 1;
noise = r/m*10^(-SNR/10);

ntests = 10;
% methods to compare:  COGENT, MFOCUSS
lambda_focuss = 2*noise^2;
tau_cogent    = r*sqrt(L);


Phi = randn(m,n);
X0 = zeros(n,L);
X0(act_inds,:) = randn(r,L);
N = noise*randn(m,L);
Y = Phi*X0 + N;

% MFOCUSS
% [X, gamma_ind, gamma_est, count] = MFOCUSS(Phi, Y, lambda_focuss,...
%     'p',1,...
%     'MAX_ITERS',1000);
% 
% emfocuss = norm(X-X0);

% COGENT

% form a vector for all
groups = sparse(L*n,n);
for ii = 1:n
    g = L*[0:L-1] + ii;
    groups(g,ii) = 1;
end

selfun = @(gradf) find_next_atom_group(gradf,groups);
AINIT= [ones(1,L); zeros(n-1,L)];
Ainit= AINIT(:);
x0 = X0(:);
phi= sparse(m*L,n*L);
for ii = 1:L
   phi((ii-1)*m+1:ii*m, (ii-1)*n+1:ii*n) = Phi ;
end
nn = N(:);
y = phi*x0 + nn;
[x, Atc, objc, timec, back_countc] = CoGEnT(y, phi, tau_cogent, Ainit, selfun ,...
    'maxiter',n);

XCOG = reshape(x,n,L);

[x, At, obj, time, back_count] = CoGEnT(y, phi, tau_cogent, Ainit, selfun ,...
    'maxiter',n,'gp_forward',0,'backward',0);

XCG = reshape(x,n,L);






