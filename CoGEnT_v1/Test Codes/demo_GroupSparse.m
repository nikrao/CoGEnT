%%%% CoGEnT demo for recovery of (overlapping) group sparse vectors

clear;
clc;


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% OVERLAPPING GROUPS %%%%%
% we look to recover wavelet coefficients of 1D signals, with parent-child
% groups as defined in
% Rao et. al "Convex Approaches to Model Wavelet Sparsity Patterns" ICIP 11
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

p = 2048; %signal length
sig = MakeSignal('Piece-Polynomial',p);
% replace with signal of choice. See MakeSignal.m for more options

% normalize the signal
sig = sig-min(sig);
sig = sig/max(sig);

[xtrue, strue] = wavedec(sig,log2(p),'haar');
xtrue = xtrue';
power = log2(p);
G = cell(0);
for i = 2:2^(power-1)
    
    group1 = {[i (2*i - 1)]};
    group2 = {[i 2*i]};
    G      = [G ; group1 ;group2];
    
end
G = [G;{1}];

% measurements
m = floor(p/4);
Phi = randn(m,p)/sqrt(m);
noise = 0.01;
y = Phi*xtrue + noise*randn(m,1);

% CoGEnT inputs
maxiter = 500;
Ainit = [1 ; zeros(p-1,1)];  % first index active
groups = sparse(p,length(G));
for ii = 1:length(G)
   inds = G{ii};
   groups(inds,ii) = 1;
end
selfun =@(gradf) find_next_atom_group(gradf,groups);
tau = 300; % might need tuning

errtaus = [];
rectaus = [];

[x, ~, ~, ~] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun ,...
    'maxiter',maxiter);
err = norm(x-xtrue)^2/p;

sigrec = waverec(x',strue,'haar');

figure, 
stem(xtrue);
hold on
plot(x,'r .');
title('Recovered DWT Coefficients');
legend('True','Recovered');


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%% NON OVERLAPPING GROUPS %%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% consider a toy block sparse vector
M = 100; % number of groups
B = 10;  % group size
k = 5;   % number of active groups
xtrue = zeros(M*B,1);
active_groups = randsample(M,k);
for ii = 1:k
   inds = (active_groups(ii)-1)*B+1:active_groups(ii)*B;
   xtrue(inds) = randn(length(inds),1);
end

m = floor(M*B/4);
Phi = randn(m,M*B)/sqrt(m);
y = Phi*xtrue;

groups = sparse(M*B,M);
for ii = 1:M
    inds = (ii-1)*B+1:ii*B;
    groups(inds,ii) = 1;
end

% CoGEnT parameters
maxiter = 500;
Ainit = zeros(M*B,1);
Ainit(inds) = 1;
Ainit = Ainit/norm(Ainit);
selfun =@(gradf) find_next_atom_group(gradf,groups);
tau = 20;

[x, ~, ~, ~] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun ,...
        'maxiter',maxiter);
    
err  = norm(x-xtrue)^2/p;

figure, 
stem(xtrue);
hold on
plot(x,'r .');
title('Recovered Block Sparse Signal');
legend('True','Recovered');





