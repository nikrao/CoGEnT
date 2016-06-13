clear
clc

% Tensor recovery using square deal
%using tensor toolbox. See for notes:
%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.103.8658&rep=rep1&type=pdf

n=10;

sampling_range = .1:.1:1;
num_tests = 10;

maxiter = 20;
tol = 1e-5;
eta = 0.5;
gptol = 1e-5;
gpiter = 50;
do_fw = 1;
gp_forward = 1;
backward = 1;
sparsify = 1;
E = [];


%Define tensor
A=randn(n,n);
r = 3;
[U1,S1,V1]=svd(A);
S2=diag(S1);
S=S2(1:r);
U=U1(:,1:r);

T=ktensor(S,U,U,U);

Tfull=full(T); % This is the true solution

xtrue = Tfull(:); %vectorize the matrix

sam_frac=.8;

% create the matrix mask
M = rand(n,n,n);
M = M<=sam_frac;
% Y = M.*X;
Phi = sparse(diag(M(:))); % this is the mask we will use

% noise parameter
noisevar = 0;

% measurements
y = Phi*xtrue;
%Y = M.*Xtrue;
Y=reshape(y,10,10,10);

%form the fibers and create matrix
% get mode 1 fiber
count=0;
for j=1:n
    for k=1:n
        count=count+1;
        A(:,count)=Y(:,j,k);
    end
end

%reshape to "square" shape? No, it will not help! X_square=A

B=sparse(A);

[Trec iter svp] = inexact_alm_mc(B',1e-4);
Trec2=Trec.U*Trec.V';

Trec3=reshape(Trec2,n,n,n);

T4=double(T);
err=T4-double(Trec3);
