% CoGEnT for symmetric orthogonal tensor completion

clear ;
clc;
close all;
warning off;

%using tensor toolbox. See for notes:
%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.103.8658&rep=rep1&type=pdf

n=100;

%Define tensor
A=randn(n,n);
r = 3;
[U1,S1,V1]=svd(A);
S2=diag(S1);
S=S2(1:r);
U=U1(:,1:r);

T=ktensor(S,U,U,U);

Tfull=full(T); % This is the true solution


sam_frac = 0.7; %sampling fraction

xtrue = Tfull(:); %vectorize the matrix

% create the matrix mask
M = rand(n,n,n);
M = M<=sam_frac;
Phi = sparse(diag(M(:))); % this is the mask we will use

% noise parameter
noisevar = 0;

% measurements
y = Phi*xtrue;

% CoGEnT PARAMETERS
selfun =@(gradf) find_next_atom_tensor(gradf,n);
uinit = randn(n,1);
uinit = uinit/norm(uinit);
Ainit1 = ktensor(1,uinit,uinit,uinit);
Ainit2=full(Ainit1);
Ainit=Ainit2(:);

maxiter = 20;
svt = 2;
tic;

tau=sum(S)+.001;

tic
[x, At, iter, obj, time, back_count] = CoGEnT_TC(y, Phi, tau, Ainit, [n n n] ,selfun,...
    'maxiter',maxiter,'svt',svt);

toc
err=norm(xtrue-x)^2/numel(x);


