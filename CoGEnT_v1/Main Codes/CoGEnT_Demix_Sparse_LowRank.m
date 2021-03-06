function [x1,x2,At1, At2, iter, obj, time, back_count] =...
    CoGEnT_Demix_Sparse_LowRank(y, Phi, tau1, tau2,Ainit1, Ainit2, selfun1, matsize,varargin)

% usage [x, iters, obj, time] = CoGEnT_Demix_lowrank(y, Phi, tau1,tau2,Ainit1,Ainit2, selfun1 , varargin)
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x1,2 0.5|| y - Phi*(x1+x2) ||^2  s.t
% ||x1||_A1 \leq tau1
% ||x2||_*  \leq tau2
% where ||.||_* is the nuclear norm of x2.
% This code uses the "recomputing the bases" backward step
% used for matrix completion problems. selful2 is the nuclear norm
% selection function by default.
%
% INPUTS:
%
% Ainit       = first atom selected
% selfun      = function handle to pick an atom. one of its inputs must be the gradient.
% y, Phi      = observarion and data matrix
% tau         = constraint on the atomic norm
%
% OPTIONAL INPUTS
% 'maxiter'   = maximum number of iterations allowed default = 10000
% 'tol'       = tolerance parameter
% 'backward'  = default = 1. If active, multiple backward steps are performed
% 'eta'       = backward parameter. Default = 0.5
% 'gptol'     = tolerance for gradient projection. Default = 1e-3
% 'gpiter'    = maximum iterations of gradient projection, Default = 10;
% 'verbose'   = default = 0. If set, gives text feedback
% 'dropcount' = default = 9999. The value will determine the number of
%               times an atom can be picked and dropped before we declare convergence
% 'gp_forward'= default = 1 perform gradient projection in forward step
% 'sparsify'  = default = 0. If set, performs a final sparsification step
%               by merging atoms and deleting zeros
% 'debias'    = default = 1. If set, performs a final debiasing step
%
% OUTPUTS :
% x           = final result
% iters       = iterations performed
% obj         = objective function value
% time        = time taken since start of the iterations
% At          = final atomic set
% back_count  = number of backward steps taken in total
%
% For more information, please refer to CoGEnT_Demix, and CoGEnT_MC
%
% %%%%%%%%%%%%
% Nikhil Rao : Nov 5 2013
% Parikshit Shah : 
% %%%%%%%%%%%%


% set optional parameter defaults
maxiter   = 10000;
tol       = 1e-5;
backward  = 1;
eta       = 0.5;
gptol     = 1e-2;
gpiter    = 10;
verbose   = 0;
dropcount = 9999;
gp_forward= 1;
sparsify  = 0;
debias    = 1;

% check for optional parameters entered by the user
if (rem(length(varargin),2)==1)
    error('Optional parameters should be in pairs');
else
    for i=1:2:(length(varargin)-1)
        switch lower(varargin{i})
            case 'maxiter'
                maxiter = varargin{i+1};
            case 'tol'
                tol = varargin{i+1};
            case 'backward'
                backward = varargin{i+1};
            case 'eta'
                eta = varargin{i+1};
            case 'gptol'
                gptol = varargin{i+1};
            case 'gpiter'
                gpiter = varargin{i+1};
            case 'verbose'
                verbose = varargin{i+1};
            case 'dropcount'
                dropcount = varargin{i+1};
            case 'gp_forward'
                gp_forward = varargin{i+1};
            case 'sparsify'
                sparsify = varargin{i+1};
            case 'debias'
                debias = varargin{i+1};
            case 'svt'
                svt = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% initialize iterate
iter = 1;
x1 = Ainit1*tau1 ;
% At1 = zeros(length(x1),maxiter);
At1(:,1) = Ainit1;
coeft1 = zeros(maxiter,1);
coeft1(1) = tau1;
history1 = zeros(size(Phi,1),maxiter);
x2 = Ainit2*tau2 ;
% At2 = zeros(length(x2),maxiter);
At2(:,1) = Ainit2;
coeft2 = zeros(maxiter,1);
coeft2(1) = tau2;
n1 = matsize(1); n2 = matsize(2);
selfun2 =@(gradf) find_next_atom_nucnorm(gradf,n1,n2);


currlength1 = iter;
currlength2 = iter;
obj = 9999*ones(1,maxiter);
obj(1) = 0.5 * norm(y - Phi*(x1+x2))^2;
back_count = 0;
if verbose
    fprintf('\n Beginning iterations \n')
end

drops = 0;

time = zeros(1,maxiter);
tic;
time(1) = toc;
while iter <= maxiter
    iter  = 1+iter;
    
    %%%%% FORWARD STEP FOR ATOMIC NORM 1 %%%%%%
    resid = Phi*(x1+x2)-y; % residual
    gradf = Phi.'*resid;  % gradient
    currlength1 = currlength1 + 1;
    
    % pick next atom
    next_atom = selfun1(gradf);
    At1=[At1, next_atom];
    % find alpha -- line search
    alpha1=linesearch2(y-x1-x2,tau1*next_atom-x1); %1 D line search alpha \in [0,1]
%     x1=x1+alpha1*(tau1*next_atom-x1);
    [coeft_n]=grad_proj_pari2(At1,y-x2,tau1);
    x1=At1*coeft_n;
    obj_fwd = 0.5 * norm((x1+x2)-y)^2;
    
    %     %%%%%% BACKWARD STEP FOR ATOMIC NORM 1  %%%%%%

                resid = Phi*(x1+x2)-y; % residual
                gradfb =  Phi'*resid;
                temp5=gradfb'*At1;
                [~, ind_rem]=min(abs(temp5));
                clear temp5;

                A_new=At1;
                A_new(:,ind_rem)=[];
                [coeft_new]=grad_proj_pari2(A_new,y-x2,tau1);
    
    
                x_rem = A_new*coeft_new;
                obj_back = 0.5* norm(y - Phi*(x_rem+x2))^2;
    
                if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
                    % backward step taken
                    currlength1 = currlength1-1;
                    At1 = A_new;
                    x1 = x_rem;
                    back_count = back_count +1;
                    if verbose
                        fprintf('#backward steps atom 1 = %d, iterations = %d \n', back_count,iter);
                    end
                    obj(iter)=obj_back;
                else
                    obj(iter) = obj_fwd;
                end
    
    
            %end
        %end
    clear next_atom;
    
    %%%%% FORWARD STEP FOR ATOMIC NORM 2 %%%%%%
    resid = Phi*(x1+x2)-y; % residual
    gradf = Phi.'*resid;  % gradient
    currlength2 = currlength2 + 1;
    
    % pick next atom
    next_atom = selfun2(gradf);
    %svd(reshape(next_atom,20,20))
    
    
    % conditional gradient step
   % alpha2=linesearch(y-x1-x2,tau2*next_atom-x2); %1 D line search alpha \in [0,1]
    %x2=x2+alpha2*(tau2*next_atom-x2);
    %obj_fwd = 0.5 * norm((x1+x2)-y)^2;
    
    At2=[At2, next_atom];

    [coeft_n]=grad_proj_pari2(At2,y-x1,tau2);
    x2=At2*coeft_n;
    
    'fwd step 2';
    if ~backward
        % update output variables
        time(iter) = toc;
        obj(iter)  = obj_fwd;
        continue
    end
    

    
    %         %%%%%% BACKWARD STEP FOR NUCLEAR NORM  %%%%%%
    %
    Tx = reshape(x2,n1,n2);
    [U,S,V] = svd(Tx,'econ');
    S = diag(S);

    cutoff = find(diff(-S) == max(diff(-S)));
    S = S(1:cutoff);
    U = U(:,1:cutoff);
    V = V';
    V = V(1:cutoff,:);


numats = length(S);
A_new = zeros(rows(U)*cols(V),numats);
for ats = 1:numats
    a = U(:,ats)*V(ats,:);
    A_new(:,ats) = a(:);
end
[coeft_new]=grad_proj_pari2(A_new,y-x1,tau2);
x_rem = A_new*coeft_new;
obj_back = 0.5* norm(y - (x1+x_rem))^2;
if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
    currlength2 = size(A_new,2);
    At2(:,1:currlength2) = A_new;
    x2 = x_rem;
    coeft2(1:currlength2) = coeft_new;
    
    obj(iter) = obj_back;

else
    obj(iter) = obj_fwd;

end
%     %%%%%%%%end backward step code

%obj(iter) = obj_fwd;


if (obj(iter) < tol)  % convergence check
    if verbose
        fprintf('\n Convergence achieved. Final objective value = %f \n',obj(iter));
    end
    break
end

clear next_atom;
end

end

