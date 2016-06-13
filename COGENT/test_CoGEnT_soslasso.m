% testbench for SOSlasso multitask learning

clear;
clc;
close all;

M = 100;
B = 20;
k = 10;
tasks = 1;
overlap = true;
% form groups
G = cell(M,1);
for ii = 1:M
    if overlap
        if ii==1
            G(ii) = {(ii-1)*B+1:ii*B};
            
        else
            G(ii) = {lastind+1:lastind+B};
        end
        lastind = floor(ii*B*2/3);
    else
        G(ii) = {(ii-1)*B+1:ii*B};
    end
end
p = G{end};
p = max(p);


thresh = 0.85;
% decide what groups are active
acts = randsample(M,k);
% form groups matrix
groups = sparse(p*tasks,M);
xstar = zeros(p,tasks);
mat = p*(0:tasks-1);
mat = repmat(mat,B,1);
for ii = 1:M
    inds = G{ii}';
    if intersect(ii,acts)
        vec = 2*rand(B,tasks)-1;
        vec = vec.*(rand(B,tasks)>thresh);
        xstar(inds,:) = vec + xstar(inds,:);
    end
    inds = repmat(inds,1,tasks)+mat;
    groups(inds(:),ii) = 1;
end

% noise and measurements
m = ceil(p/3);
Phi = randn(m,p)/sqrt(m);
noise = 0.01;
y = Phi*xstar + noise*randn(m,tasks);

% COGENT params
maxiter = 200;
eta = 0.5;
tol = max(noise,1e-8);
Ainit = [ones(1,tasks); zeros(p-1,tasks)];
selfun =@(gradf) find_next_atom_soslasso_multitask(gradf,groups);
errtau = [];rectau = [];

tauset = 2.^(2:5);

for tau = tauset
    [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,...
        'tol',tol,...
        'eta',eta,...
        'debias',0,...
        'backward',1);
    
    errtau = [errtau, norm(x-xstar,'fro')^2/numel(xstar)];
    rectau = [rectau x];
    
    fprintf('.');
end
fprintf('\n');
[err,i] = min(errtau);
x = rectau(:,(i-1)*tasks+1:i*tasks);
subplot(121)
imagesc(xstar~=0); 
subplot(122)
imagesc(x~=0);

[ecogent, i] = min(errtau);
x = rectau(:,i);
indstrue = find(xstar~=0);
valstrue = xstar(indstrue);
indscogent = find(x);
valscogent = x(indscogent);
stem(indstrue,valstrue);
hold on
plot(indscogent,valscogent,'r .');
legend('True','CoGEnT')
