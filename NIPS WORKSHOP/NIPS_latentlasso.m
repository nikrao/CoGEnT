clear;
clc;
close all;

% latent group lasso for high dimensional regression. For the NIPS workshop

M = 50;
B = 50;
k = floor(M/10);

% form groups
G = cell(M,1);
G(1) = {1:B};
lastind = B;

for ii = 2:M
    G(ii) = {lastind-29:lastind-30+B};
    lastind = lastind-30+B;
end
p = lastind;
n = ceil(p/2);

Phi = randn(n,p)/sqrt(n);

% true signal and observations
acts = randsample(M,k);
xstar = sparse(zeros(p,1));
nx = 0;
for grp = 1:k
    inds = G{acts(grp)};
    v = (2*rand(B,1)-1);
    xstar(inds) = xstar(inds) + v;
    nx = nx + norm(v);
end
noise = 0.1;
y = Phi*xstar + noise*randn(n,1);

Ainit = [ones(B,1)/sqrt(B); zeros(p-B,1)];
errtausfb = []; errtausfw = [];
rectausfb = []; rectausfw = [];
timesfb = []; timesfw = [];

% groups matrix
groups = sparse(zeros(p,M));
for g = 1:M
    gind = G{g};
    groups(gind,g) = 1;
end
selfun =@(gradf) find_next_atom_group(gradf,groups);
for tau = nx.*1.25
    
    [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
        'maxiter',5000,'tol', noise/sqrt(n),'gp_forward',0, 'gptol',1e-1, 'gpiter', 1);
    timesfb = [timesfb ; time(end)];
    errtausfb = [errtausfb norm(x-xstar)^2/p];
    rectausfb = [rectausfb sparse(x)];
    
    fprintf('FB \t')
    
%     [x, ~, iter, ~, time, ~] = FOBA_L2(y, Phi, tau, Ainit, selfun,...
%         'backward',0, 'tol', 1e-8,'maxiter',1000);
%     
%     timesfw = [timesfw ; time(iter)];
%     errtausfw = [errtausfw norm(x-xstar)^2/p];
%     rectausfw = [rectausfw sparse(x)];
%     
%     fprintf('FW \n')
end


% setup for latent lasso
[repindex, ginds, groupsR] = get_replication_data(G);
Phit = Phi(:,repindex);
lamset = linspace(0.1,1,5);
fprintf('matrix replicated \n')
errlams = [];
reclams = [];
times   = [];

%
% pt = uint32(sum(cellfun(@length,G)));
% groupinds = sparse(zeros(p,pt));
% % row i of groupinds has a 1 in the locations where the ith entry in the
% % ambient space is in the location in the replicated space
% for ii = 1:p
%     ind = find(repindex == ii);
%     groupinds(ii,ind) = 1;
% end

for lam = 10.^[-2:1]
    
    
    %     x = LatentLasso(y,Phit,Phi,lam,ginds,...
    %     groupsR,G,groupinds',repindex);
    [x, time] = LatentLasso(y,Phit,Phi,lam,ginds,...
        groupsR,G, 5000, noise/sqrt(n));
    times = [times; time];
    errlams = [errlams; norm(x-xstar)^2/p];
    reclams = [reclams sparse(x)];
    fprintf('LGL \n')
end
clear Phit
% setup for FB and CGD methods


% save NIPS_LATENT errtausfb errtausfw rectausfb rectausfw timesfb timesfw times errlams reclams p M



