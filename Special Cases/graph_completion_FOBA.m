% generate a toy graph, obscure data and look to interpolate it back

clear;
clc;
close all;

N = 100; % number of nodes in the graph
S = rand(N);
S = S.*(rand(N)>0.5); % inverse coveriance matrix of the graph
C = inv(S);  % covariance matrix of the variables

num_samples = 25; %number of samples per node
VALUES = real(mgd(num_samples,N,zeros(N,1),C));

% randomly blank out some data
blank = 0.5;
mask = rand(size(VALUES));
DATA = VALUES.*(mask>=blank);

% form groups
G = cell(0);
for i = 1:N
    inds = find(S(i,1:i)~=0);
    inds = setdiff(inds,i);
    if ~isempty(inds)
        gtemp = repmat(i,length(inds),1);
        gtemp = [gtemp inds'];
        gtemp = mat2cell(gtemp,ones(length(inds),1),2);
        G = [G;gtemp];
    end
    fprintf('.')
end
fprintf('\n')

%FOBA parameters
tauset = .2:.2:1;
gtemp = G{1};
Ainit = zeros(N,1);
Ainit(gtemp) = ones(length(gtemp),1);
Ainit = Ainit./norm(Ainit);
selfun =@(gradf) find_next_atom_group(gradf,G);
maxiter = 100;
eta = 0.5;
gptol = 1e-3;
gpiter = 20;
dropcount = 20;
do_fw = 1;
gp_forward = 1;
errtau = [];
for tau = tauset
    etau = 0;
    for node = 1:N
        this_data = DATA;
        y = this_data(:,node);
        this_data(:,node) = zeros(rows(DATA),1);
        
        [x, At,iter, obj, time, backs] = FOBA_L2(y, this_data, tau, Ainit, selfun, ...
            'mbackward', 1, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'dropcount',dropcount,...
            'do_fw',do_fw,...
            'gp_forward',gp_forward,...
            'verbose', 0);
        
        xfin = this_data*x;
        
        etau = etau + (norm(xfin - VALUES(:,node))^2)/numel(x);
        fprintf('.')
    end
    errtau = [errtau etau];
    fprintf('\n tau = %f \n', tau);
end

[e i] = min(errtau);

RECVALS = zeros(size(VALUES));
for tau = tauset(i)
    for node = 1:N
        this_data = DATA;
        y = this_data(:,node);
        this_data(:,node) = zeros(rows(DATA),1);
        
        [x, At,iter, obj, time, backs] = FOBA_L2(y, this_data, tau, Ainit, selfun, ...
            'mbackward', 1, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'dropcount',dropcount,...
            'do_fw',do_fw,...
            'gp_forward',gp_forward,...
            'verbose', 0);
        
        RECVALS(:,node) = this_data*x;
        fprintf('.')
    end
end
    

