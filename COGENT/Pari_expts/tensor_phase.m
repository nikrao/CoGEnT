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

n=30;

sampling_range = .1:.1:1;
num_tests = 3;

maxiter = 20;
tol = 1e-5;
eta = 0.5;
gptol = 1e-5;
gpiter = 50;
do_fw = 1;
gp_forward = 1;
backward = 1;
sparsify = 1;
E = [];E1 = [];
for sam_frac = sampling_range
    
    errtaus = [];errtaus1 = [];
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
        Phi = sparse(diag(M(:))); % this is the mask we will use
        
        % noise parameter
        noisevar = 0;
        
%         measurements
        y = Phi*xtrue;
        
%         FOBA PARAMETERS
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
        
        
        M = rand(n,n,n);
        M = M<=sam_frac;
        J = M(:);
        Phi = spalloc(n^3,n^3,n^3);
        for ii = 1:n^3
            Phi(ii,ii)= J(ii);
        end
        clear J;
        
        y = Phi*xtrue; %add noise here, choose tau
        Y=reshape(y,n,n,n);
        
        A = zeros(n,n^2);
        count=0;
        for j=1:n
            for k=1:n
                count=count+1;
                A(:,count)=Y(:,j,k);
            end
        end
        clear Y
        
        tau = sum(svd(A));
        
        A = A(:);
        actind = find(A ~= 0);
        Phi = spalloc(n^3,n^3,n^3);
        for ii = 1:length(actind)
           Phi(actind(ii),actind(ii)) = 1; 
        end
        y = Phi*A; %add noise here, choose tau
       
        selfun = @(gradf) find_next_atom_nucnorm(gradf,n,n^2);
        Ainit = randn(n,1)*randn(1,n^2);
        Ainit = Ainit(:)/norm(Ainit(:));
        [x, ~, iterc, ~, timec, ~] = CoGEnT_MC(y, Phi, tau, Ainit, [n,n^2] ,selfun,...
            'maxiter',maxiter);
        
        
        errtaus1 = [errtaus1 (norm(xtrue(:)-x(:))^2/numel(x))<1e-4];
        
    end
    fprintf('\n');
    E = [E ; errtaus];
    E1 = [E1 ; errtaus1];
end
save tensorPhaseM1 E E1