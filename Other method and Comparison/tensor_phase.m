% tensor completion phase transition

clear ;
clc;
close all;
warning off;

%%%%%%%%%% Problem setup
clear
clc

%using tensor toolbox. See for notes:
%http://citeseerx.ist.psu.edu/viewdoc/download?doi=10.1.1.103.8658&rep=rep1&type=pdf

n=100;

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
for sam_frac = sampling_range
    
    errtaus = [];
    for test = 1:num_tests
        
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
        
        % create the matrix mask
        M = rand(n,n,n);
        M = M<=sam_frac;
        % Y = M.*X;
        Phi = speye(n^3);
        m = M(:);
        for index = 1:length(m)
           Phi(index,index)=m(index); 
           if mod(index,5000)==0
               fprintf('~');
           end
        end
        fprintf('\n');
        % noise parameter
        noisevar = 0;
        
        % measurements
        y = Phi*xtrue;
        %Y = M.*Xtrue;
        
        % FOBA PARAMETERS
        selfun =@(gradf) find_next_atom_tensor(gradf,n);
        uinit = randn(n,1);
        uinit = uinit/norm(uinit);
        Ainit1 = ktensor(1,uinit,uinit,uinit);
        Ainit2=full(Ainit1);
        Ainit=Ainit2(:);
        
        svt = 2;
        
        tau=sum(S)+.001;
        
        % CoGEnT with SVT in the Backward Step
        [x, At, iter, obj, time, back_count] = CoGEnT_TC(y, Phi, tau, Ainit, [n n n] ,selfun,...
            'maxiter',maxiter,'svt',svt);
        
        errtaus = [errtaus (norm(xtrue-x)^2/numel(x))<1e-4];
        
        fprintf('.');
    end
    fprintf('\n');
    E = [E ; errtaus];
end
save tensorPhase_BIG E