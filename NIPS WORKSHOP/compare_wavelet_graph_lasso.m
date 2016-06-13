% script to compare CG with SPARSA and CoGEnT

clear
clc
close all

for i = 1
% signal
signal = 'HeaviSine';
p = 1024; % signal length
n = 300; % number of measurements
noise = 0.01; % noise std. dev
Phi = randn(n,p)/sqrt(n); % measurement matrix

% groups etc
grouping = 'pc';
group_MAT = FormGroups(grouping,1,p);
groups = sparse(p,length(group_MAT));
for g = 1:length(group_MAT)
    G = group_MAT{g};
    groups(G,g) = 1;
end


% make signal and obtain measurements
xstar = MakeSignal(signal,p);
xstar = xstar - min(xstar);
xstar = xstar/max(xstar);
xstar = xstar';
W = MakeHaarBases(p,1); % (inverse) wavelet matrix
s = W'*xstar; % waveler coefficients
y = Phi*s + noise*randn(n,1);

% algorithm parameters
maxiter = 1000;
tol = 1e-6;
tauset = norm(s,1)*linspace(0.5,2,25);
lamset = linspace(0.001,5,25);
selfun = @(gradf) find_next_atom_group(gradf,groups);
Ainit = [1; zeros(p-1,1)];


% cogent and cg
SrecCOG = zeros(p,length(tauset));
SrecCG = zeros(p,length(tauset));
tind = 0;
for tau = tauset
    tind = 1+tind;
    [SrecCOG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol);
    
    [SrecCG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol,'gp_forward',0, 'backward',0);
    
    fprintf('.');
    
end
fprintf('\n');

% % sparsa
% SrecSP = zeros(p,length(lamset));
% [repindex, ginds, groupsR] = get_replication_data(group_MAT);
% Phit = Phi(:,repindex);
% lind = 0;
% for lam = lamset
%     lind = 1:lind;
%     [SrecSP(:,lind) , time] = LatentLasso(y,Phit,Phi,lam,ginds,...
%     groupsR,group_MAT,maxiter, tol);
%     
%     fprintf('.');
% end
% fprintf('\n');


% compute errors
XhatCOG = W*SrecCOG;
XhatCG = W*SrecCG;
% XhatSP = W*SrecSP;
Xcog = repmat(xstar,1,length(tauset));
% Xsp  = repmat(xstar,1,length(lamset));

ecog = sum((Xcog-XhatCOG).^2)/p;
ecg  = sum((Xcog-XhatCG).^2)/p;
% esp  = sum((Xsp-XhatSP).^2)/p;

[cog, icog] = min(ecog);
[cg, icg] = min(ecg);
% [sp, isp] = min(esp);

cogent = XhatCOG(:,icog);
cg = XhatCG(:,icg);
% sp = XhatSP(:,isp);

nnz_cogent = nnz(SrecCOG(:,icog));
nnz_cg = nnz(SrecCG(:,icg));
% nnz_sp = nnz(SrecSP(:,isp));
nnz_true = nnz(s);

save WAVHSINE cogent xstar cg nnz_cogent nnz_cg nnz_true

end

fprintf('Hsine done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1
% signal
signal = 'Blocks';
p = 1024; % signal length
n = 300; % number of measurements
noise = 0.01; % noise std. dev
Phi = randn(n,p)/sqrt(n); % measurement matrix

% groups etc
grouping = 'pc';
group_MAT = FormGroups(grouping,1,p);
groups = sparse(p,length(group_MAT));
for g = 1:length(group_MAT)
    G = group_MAT{g};
    groups(G,g) = 1;
end


% make signal and obtain measurements
xstar = MakeSignal(signal,p);
xstar = xstar - min(xstar);
xstar = xstar/max(xstar);
xstar = xstar';
W = MakeHaarBases(p,1); % (inverse) wavelet matrix
s = W'*xstar; % waveler coefficients
y = Phi*s + noise*randn(n,1);

% algorithm parameters
maxiter = 1000;
tol = 1e-6;
tauset = norm(s,1)*linspace(0.5,2,25);
lamset = linspace(0.001,5,25);
selfun = @(gradf) find_next_atom_group(gradf,groups);
Ainit = [1; zeros(p-1,1)];


% cogent and cg
SrecCOG = zeros(p,length(tauset));
SrecCG = zeros(p,length(tauset));
tind = 0;
for tau = tauset
    tind = 1+tind;
    [SrecCOG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol);
    
    [SrecCG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol,'gp_forward',0, 'backward',0);
    
    fprintf('.');
    
end
fprintf('\n');

% sparsa
% SrecSP = zeros(p,length(lamset));
% [repindex, ginds, groupsR] = get_replication_data(group_MAT);
% Phit = Phi(:,repindex);
% lind = 0;
% for lam = lamset
%     lind = 1:lind;
%     [SrecSP(:,lind) , time] = LatentLasso(y,Phit,Phi,lam,ginds,...
%     groupsR,group_MAT,maxiter, tol);
%     
%     fprintf('.');
% end
% fprintf('\n');


% compute errors
XhatCOG = W*SrecCOG;
XhatCG = W*SrecCG;
% XhatSP = W*SrecSP;
Xcog = repmat(xstar,1,length(tauset));
% Xsp  = repmat(xstar,1,length(lamset));

ecog = sum((Xcog-XhatCOG).^2)/p;
ecg  = sum((Xcog-XhatCG).^2)/p;
% esp  = sum((Xsp-XhatSP).^2)/p;

[cog, icog] = min(ecog);
[cg, icg] = min(ecg);
% [sp, isp] = min(esp);

cogent = XhatCOG(:,icog);
cg = XhatCG(:,icg);
% sp = XhatSP(:,isp);

nnz_cogent = nnz(SrecCOG(:,icog));
nnz_cg = nnz(SrecCG(:,icg));
% nnz_sp = nnz(SrecSP(:,isp));
nnz_true = nnz(s);

save WAVBLOCKS cogent xstar cg nnz_cogent nnz_cg nnz_true

end

fprintf('Blocks done')

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


for i = 1
% signal
signal = 'Piece-Regular';
p = 1024; % signal length
n = 300; % number of measurements
noise = 0.01; % noise std. dev
Phi = randn(n,p)/sqrt(n); % measurement matrix

% groups etc
grouping = 'pc';
group_MAT = FormGroups(grouping,1,p);
groups = sparse(p,length(group_MAT));
for g = 1:length(group_MAT)
    G = group_MAT{g};
    groups(G,g) = 1;
end


% make signal and obtain measurements
xstar = MakeSignal(signal,p);
xstar = xstar - min(xstar);
xstar = xstar/max(xstar);
xstar = xstar';
W = MakeHaarBases(p,1); % (inverse) wavelet matrix
s = W'*xstar; % waveler coefficients
y = Phi*s + noise*randn(n,1);

% algorithm parameters
maxiter = 1000;
tol = 1e-6;
tauset = norm(s,1)*linspace(0.5,2,25);
lamset = linspace(0.001,5,25);
selfun = @(gradf) find_next_atom_group(gradf,groups);
Ainit = [1; zeros(p-1,1)];


% cogent and cg
SrecCOG = zeros(p,length(tauset));
SrecCG = zeros(p,length(tauset));
tind = 0;
for tau = tauset
    tind = 1+tind;
    [SrecCOG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol);
    
    [SrecCG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol,'gp_forward',0, 'backward',0);
    
    fprintf('.');
    
end
fprintf('\n');

% sparsa
% SrecSP = zeros(p,length(lamset));
% [repindex, ginds, groupsR] = get_replication_data(group_MAT);
% Phit = Phi(:,repindex);
% lind = 0;
% for lam = lamset
%     lind = 1:lind;
%     [SrecSP(:,lind) , time] = LatentLasso(y,Phit,Phi,lam,ginds,...
%     groupsR,group_MAT,maxiter, tol);
%     
%     fprintf('.');
% end
% fprintf('\n');


% compute errors
XhatCOG = W*SrecCOG;
XhatCG = W*SrecCG;
% XhatSP = W*SrecSP;
Xcog = repmat(xstar,1,length(tauset));
% Xsp  = repmat(xstar,1,length(lamset));

ecog = sum((Xcog-XhatCOG).^2)/p;
ecg  = sum((Xcog-XhatCG).^2)/p;
% esp  = sum((Xsp-XhatSP).^2)/p;

[cog, icog] = min(ecog);
[cg, icg] = min(ecg);
% [sp, isp] = min(esp);

cogent = XhatCOG(:,icog);
cg = XhatCG(:,icg);
% sp = XhatSP(:,isp);

nnz_cogent = nnz(SrecCOG(:,icog));
nnz_cg = nnz(SrecCG(:,icg));
% nnz_sp = nnz(SrecSP(:,isp));
nnz_true = nnz(s);

save WAVPREG cogent xstar cg nnz_cogent nnz_cg nnz_true

end

fprintf('PREG done')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



for i = 1
% signal
signal = 'Piece-Polynomial';
p = 1024; % signal length
n = 300; % number of measurements
noise = 0.01; % noise std. dev
Phi = randn(n,p)/sqrt(n); % measurement matrix

% groups etc
grouping = 'pc';
group_MAT = FormGroups(grouping,1,p);
groups = sparse(p,length(group_MAT));
for g = 1:length(group_MAT)
    G = group_MAT{g};
    groups(G,g) = 1;
end


% make signal and obtain measurements
xstar = MakeSignal(signal,p);
xstar = xstar - min(xstar);
xstar = xstar/max(xstar);
xstar = xstar';
W = MakeHaarBases(p,1); % (inverse) wavelet matrix
s = W'*xstar; % waveler coefficients
y = Phi*s + noise*randn(n,1);

% algorithm parameters
maxiter = 1000;
tol = 1e-6;
tauset = norm(s,1)*linspace(0.5,2,25);
lamset = linspace(0.001,5,25);
selfun = @(gradf) find_next_atom_group(gradf,groups);
Ainit = [1; zeros(p-1,1)];


% cogent and cg
SrecCOG = zeros(p,length(tauset));
SrecCG = zeros(p,length(tauset));
tind = 0;
for tau = tauset
    tind = 1+tind;
    [SrecCOG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol);
    
    [SrecCG(:,tind), ~, ~, ~, ~] = CoGEnT(y, Phi, tau, Ainit, selfun,...
        'maxiter',maxiter,'tol',tol,'gp_forward',0, 'backward',0);
    
    fprintf('.');
    
end
fprintf('\n');

% sparsa
% SrecSP = zeros(p,length(lamset));
% [repindex, ginds, groupsR] = get_replication_data(group_MAT);
% Phit = Phi(:,repindex);
% lind = 0;
% for lam = lamset
%     lind = 1:lind;
%     [SrecSP(:,lind) , time] = LatentLasso(y,Phit,Phi,lam,ginds,...
%     groupsR,group_MAT,maxiter, tol);
%     
%     fprintf('.');
% end
% fprintf('\n');


% compute errors
XhatCOG = W*SrecCOG;
XhatCG = W*SrecCG;
% XhatSP = W*SrecSP;
Xcog = repmat(xstar,1,length(tauset));
% Xsp  = repmat(xstar,1,length(lamset));

ecog = sum((Xcog-XhatCOG).^2)/p;
ecg  = sum((Xcog-XhatCG).^2)/p;
% esp  = sum((Xsp-XhatSP).^2)/p;

[cog, icog] = min(ecog);
[cg, icg] = min(ecg);
% [sp, isp] = min(esp);

cogent = XhatCOG(:,icog);
cg = XhatCG(:,icg);
% sp = XhatSP(:,isp);

nnz_cogent = nnz(SrecCOG(:,icog));
nnz_cg = nnz(SrecCG(:,icg));
% nnz_sp = nnz(SrecSP(:,isp));
nnz_true = nnz(s);

save WAVPPOLY cogent xstar cg nnz_cogent nnz_cg nnz_true

end

fprintf('Ppoly done')
