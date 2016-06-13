function [thetaHATgrp , time] = LatentLasso(y,A,Aold,tau,groups,...
    group_arr,group_MAT,maxiter, tol)

% this function solves the overlap group lasso method, with replication of
% columns of A as explained in Jacob et.al.

%INPUTS :
% y         = the observed vector
% A         = the measurement matrix after replication
% Aold      = original measurement matrix for debiasing
% tau       = regularization parameter for the overlap lasso penalty
% groups    = row vector indicating groups
% group_arr = matrix with each row indexing group elements in the replicate
% group_MAT = cell array having groups as rows


[~,n] = size(Aold);

% MATRIX OPERATION
% hA = @(x) Ax(x,Aold,groupsINDS); %A*x;
% hAt = @(x) Atx(x,Aold,repindex); %A'*x;
hA = @(x) A*x;
hAt = @(x) A'*x;

% REGULARIZER

psi = @(x,tau) l1l2vectorsoft(x,tau,groups,group_arr);
phi = @(x) l1l2norm(x,group_arr);


% GROUP LASSO
tic
[thetaGrp,xg,~,~,~,~,~]= ...
    SpaRSA(y,hA,tau,...
    'At',hAt,...
    'Psi',psi,...
    'Phi',phi,...
    'Debias',1,...
    'StopCriterion',4,...
    'Monotone',1,...
    'Continuation',1,...
    'MaxiterA',maxiter, ...
    'ToleranceA',tol,...
    'Verbose',0 ...
    );

% RECONSTRUCTION


% THIS IS THE DEBIASING STEP FOR GROUP LASSO WITH
% REPLICATION

% 1. find groups that have been selected
g = []; % this indicates the groups that are selected
keep = find(thetaGrp ~= 0);
%if is empty, then it causes problems, hence, let keep correspond
%to all the elements in thetaHAT
if isempty(keep)
    keep = 1:length(thetaGrp);
end

%  determine the "group numbers" to select. i.e. the row number of
%  the group_MAT matrix that corresponds to the non zero element
%  locations of the group LASSO estimate

g = groups(keep);

g = unique(g);% these are the groups that have been selected

% 2. now that we know the groups, collapse the repetitions , due to
% overlap
selgroups = [];
for j = 1:length(g)
    selgroups = [selgroups group_MAT{g(j)}];
end

selcoeff  = unique(selgroups); %these are the coefficients on which to run LS

% 3. run least squares on selcoeff
Als = Aold(:,selcoeff);
thetaLS = Als\y;

% 4. append this new theta into the estimate, by putting 0s wherever
% needed

thetaHATgrp = zeros(n,1);
thetaHATgrp(selcoeff) = thetaLS;

time = toc;

end
