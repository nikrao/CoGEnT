% video denoising using the SOSLASSO and FOBA

clear;
clc;
close all;
tic;
%% problem parameters
frame_range = 1051:1100;
noise = 0.001;


%% form groups
imsize = 64;
gsize = 100;
shift = 50;
G = cell(0);
for gind = 1:shift:imsize^2-gsize+1
    %     G{gind} = gind:gind+gsize-1;
    g = gind:gind+gsize-1;
    G = [G; {[g]}];
end
toc;
fprintf('\n groups made \n')
%% measurement matrix and replication
m = ceil(imsize^2/3);
p = imsize^2;
A = randn(m,p)/sqrt(m);
toc;
fprintf('\n measurement matrix created \n')
%% load frames and obtain noisy measurements
THETAS = [];
for i = frame_range
    if i>frame_range(1)
        prevI = I;
    end
    str = strcat('load walk',num2str(i),';');
    eval(str);
    I = im2double(I);
    I= imresize(I,[imsize,imsize]);
    I = I(:);
    %     [theta,struct] = wavedec(I,level,'haar');
    %     STRUCTS = [STRUCTS; {[struct]}];
    if i>frame_range(1)
        theta = I-prevI;
    else
        theta = I;
    end
    THETAS = [THETAS theta];
end
Y = A*THETAS(:,2:end) + noise*randn(m,size(THETAS(:,2:end),2));
toc;
fprintf('measurements obtained \n')

%% FOBA parameters
tauset = 10.^[1:5];
gamset = 2.^[-3:3];
maxiter = 1000;
eta = 0.5;
gptol = 1e-3;
gpiter = 10;
dropcount = 20;
do_fw = 0;
gp_forward = 1;
errtau = zeros(length(lamset),length(gamset));
rectau = cell(length(lamset),length(gamset));

%% FORM THE "VECTORIZED" martices and groups
y = Y(:);
Gmat = cell(length(G),1);
[p,T]  = size(THETAS(:,2:end));
addvec = [0:T-1]*p;
addvec = repmat(addvec,[gsize,1]);
for i = 1:length(G)
    g = G{i};
    g = g';
    g = repmat(g,[1,T]);
    g = g+addvec;
    g = g(:)';
    Gmat(i) = {g};
end
Phi = sparse(m*T,p*T);
for i = 1:T
    Phi((i-1)*m+1:i*m , (i-1)*p+1:i*p) = A;
end
%% START FOBA
M = length(Gmat);B = numel(g);
Ainit = [ones(B,1); zeros(M*B-B,1)];
Ainit = Ainit/norm(Ainit);
for tind = 1:length(tauset)
    tau = tauset(tind);
    
    for gind = 1:length(gamset)
        gamma = gamset(gind);
        selfun =@(gradf) find_next_atom_soslasso(gradf,Gmat,gamma);
        [x, At,iter, obj, time, backs] = FOBA_L2(Y, Phi, tau, Ainit, selfun, ...
            'mbackward', 1, ...
            'maxiter',maxiter,...
            'eta', eta,...
            'gptol', gptol,...
            'gpiter', gpiter,...
            'dropcount',dropcount,...
            'do_fw',do_fw,...
            'gp_forward',gp_forward,...
            'verbose', 0);
        
        fprintf('.')
        
        