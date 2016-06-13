%%%  script for solving the system ID problem using FOBA

clear;clc;close all

n = 100 ; % number of measurements
n_atoms = 10; % number of atoms in the true signal
z = 1:n;

true_thetas = rand(n_atoms,1);
true_thetas = true_thetas*2*pi;

true_a = exp(1i*true_thetas);
true_a = true_a.*rand(n_atoms,1); 
for ii = z
    true_atoms(ii,:) = (1 - abs(true_a).^2)./(ii - true_a);  % the true atoms
end

true_c = rand(n_atoms,1);

%measurements
noisestd = 0;
y = zeros(n,1);
for jj = 1:n_atoms
    y = y + true_c(jj)*true_atoms(:,jj) + noisestd*randn(n,1);
end

% ALGORITHM PARAMETERS
gridding=0;
maxiter_gp = 200;  % as of now we are not using gradient projections
tol_gp = 1e-4;
eta = 0;
usecvx = 1;   % if ON, use CVX for the optimization step
verbose = 0;  % print messages and diagnostics
backcount = 20;  % to check convergence
tol = 1e-10 ; % to check convergence
maxiter = 500; % max iterations

drops = 0;


for iter = 2:maxiter


