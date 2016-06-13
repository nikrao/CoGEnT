%%% GREEDY METHOD FOR GROUP TESTING OVER GRAPHS

clear; clc;
close all;

% create a random adjacency matrix
N = 10000; % number of nodes
thr = 0.98;
Adj = rand(N);
Adj = Adj>thr;
Adj = ((Adj + Adj')/2)~=0;
fprintf('Adjacency Matrix Created \n')

% form groups. Each group is "centered" at a node and contains all its
% neighbord
G = cell(N,1);
for node = 1:N
    vec = Adj(:,node);
    vec(node) = 0;
    G(node) = {find(vec~=0)'};
end
fprintf('Groups made \n')

% measurement matrix. Each matrix is a path along the graph, and is binary

% for now, let it be a random binary matrix
m = 4000;
A = (rand(m,N))>=0.75;
fprintf('Measurements Obtained \n')

% form a graph activation pattern. We pick a node at random, and then pick
% a few more connected nodes, and set all the associated nodes to be
% defective
num_defective = 1;
seed_node = randsample(N,1);
active_inds = [];
while num_defective>0
    g = G{seed_node};
    active_inds = [active_inds g];
    seed_node = randsample(g,1);
    num_defective = num_defective-1;
end
active_inds = unique(active_inds);
xstar = zeros(N,1);
xstar(active_inds) = 1;
fprintf('activation pattern made \n')

% obtain measurements
y = A*xstar + 0.2*randn(m,1);

% FOBA params
tau = 1;
selfun =@(gradf) find_next_atom_group_linf(gradf,G);
g = G{1};
Ainit = zeros(N,1);
Ainit(g) = 1;Ainit = Ainit/norm(Ainit);
maxiter = 10000;

fprintf('Beginning CoGEnT \n')
[rec_pattern, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun ,...
    'maxiter',maxiter,'gp_forward',1);
fprintf('CoGEnT done \n')



