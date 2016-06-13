% CODE TO PERFORM MULTITASK LEARNING USING DIRTY MODELS AND THE COGENT METHOD,
% AND COMPARING TO THE CODE IN MALSAR PACKAGE

clear;
clc;
close all;
warning off

% generate a sparse matrix to recover
T = 5; % number of tasks
index = 2;
p = 500; % dimension
kg = floor(0.01*p); % active groups
ks = floor(0.01*p); % other random coefficeints

XSTAR = zeros(p,T); 
XT1 = XSTAR;
XT2 = XSTAR;
a_g = randsample(p,kg);
XT1(a_g,:) = sign(2*rand(kg,T) - 1);
for t = 1:T
    a_s = randsample(p,ks);
    XT2(a_s,t) = sign(2*rand(ks,1) - 1);
end
XSTAR = XT1 + XT2;

% measurement matrices
noise = 0.0;
m = p;
for i = 1:T
    str = sprintf('Phi%d = eye(m,p);',i); eval(str);
    str = sprintf('w%d = noise*randn(m,1);',i); eval(str);
    str = sprintf('y%d = Phi%d * XSTAR(:,%d) + w%d;',i,i,i,i); eval(str);
end

% construct observations for FB
y = y1;
Phi = sparse(Phi1);
for i = 2:T
    str = sprintf('y = [y;y%d];',i);eval(str);
    str = sprintf('Phi = blkdiag(Phi,Phi%d);',i);eval(str);
end


% construct groups for FB vectorized
G = cell(p + T*p,1);
groups = sparse(T*p,p);

for ii = 1:p
   G(ii) = {ii+ p*(0:T-1)};
   groups(ii+ p*(0:T-1),ii) = 1;
end


% FB atom selection functions
selfun1 = @(gradf) find_next_atom_group_linf(gradf,groups);
selfun2 = @(gradf) find_next_atom_l1(gradf);
Ainit1 = zeros(T*p,1);
g = G{1};
Ainit1(g) = 1;
Ainit1 = Ainit1/norm(Ainit1);
Ainit2 = [1; zeros(T*p-1,1)];


% other FB parameters
maxiter = (kg*5+ks*5);
eta = 0.8;
gptol = 1e-5;
gpiter = 20;
dropcount = 200;
do_fw = 1;
gp_forward = 1;
backward = 1;
tauset1 = linspace(10,500,10);
tauset2 = norm(XT2(:),1); %linspace(5,100,20);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
t1 = 0;
for tau1 = tauset1
    t1 = 1+t1;
    t2 = 0;
    for tau2 = tauset2
        t2 = 1+t2;
        [x1,x2,~, ~] = CoGEnT_Demix(y, Phi, tau1,tau2,...
            Ainit1, Ainit2, selfun1, selfun2,...
            'maxiter',maxiter,...
            'eta',eta,...
            'gptol',gptol,...
            'gpiter',gpiter,...
            'dropcount',dropcount,...
            'debias',1,...
            'verbose',0);
        
        fprintf('.')
        errtau(t1,t2) = norm((x1+x2)-XSTAR(:))^2/norm(XSTAR(:))^2;
    end
    fprintf('TAU1 %d \n',t1)
end

% obtain best parameters and rerun COGENT and note time
[r c] = find(errtau == min(min(errtau))); 
tau1 = tauset1(r(1));
tau2 = tauset2(c(1));
[x1,x2,~, ~] = CoGEnT_Demix(y, Phi, tau1,tau2,...
    Ainit1, Ainit2, selfun1, selfun2,...
    'maxiter',maxiter,...
    'eta',eta,...
    'gptol',gptol,...
    'gpiter',gpiter,...
    'dropcount',dropcount,...
    'debias',1);

XFB = reshape((x1+x2),size(XSTAR));
EFB = norm(XFB-XSTAR,'fro')^2/norm(XSTAR(:))^2;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
