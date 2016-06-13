function [x, At, iter, obj, time, back_count] = CoGEnT_MC_v2(y, Phi, tau, uinit,vinit, varargin)
% usage  [x, At, iter, obj, time, back_count] = CoGEnT_MC_v2(y, Phi, tau, Ainit, varargin)
%
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - Phi x ||^2  s.t ||x||_* \leq tau
% INPUTS:
%
% Ainit       = first atom selected
% selfun      = function handle to pick an atom. one of its inputs must be the gradient.
% y, Phi      = observarion and data matrix
% tau         = constraint on the atomic norm
% matsize     = 2d vector of rows cols of matrix
%
% OPTIONAL INPUTS
% 'maxiter'   = maximum number of iterations allowed default = 10000
% 'tol'       = tolerance parameter
% 'backward'  = default = 1. If active, backward steps are performed
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
% 'svt'       = if present, performs svt in the backward step
%
% OUTPUTS :
% x           = final result
% iters       = iterations performed
% obj         = objective function value
% time        = time taken since start of the iterations
% At          = final atomic set
% back_count  = number of backward steps taken in total
%
% %%%%%%%%%%%%
% Nikhil Rao, Parikshit Shah and Stephen Wright
% %%%%%%%%%%%%

% set optional parameter defaults
maxiter   = 100;
tol       = 1e-8;
backward  = 1;
eta       = 0.5;
gptol     = 1e-2;
gpiter    = 10;
verbose   = 0;
gp_forward= 0;
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

n1 = length(uinit);
n2 = length(vinit);

% initialize iterate
iter = 1;
x = uinit*vinit*tau;

% preallocate Ut, Vt and coeft
Ut = zeros(n1,maxiter);
Ut(:,1) = uinit;
Vt = zeros(maxiter,n2);
Vt(1,:) = vinit;

coeft = spalloc(maxiter,maxiter,maxiter);
coeft(1,1) = tau;

currlength = iter;
obj = 9999*ones(1,maxiter);

obj(1) = 0.5 * norm(y - Phi.*x,'fro')^2;
if verbose
    fprintf('\n Beginning iterations \n')
    figure(1)
end

time = zeros(1,maxiter);
tic;
time(1) = toc;
y = sparse(y);
while iter <= maxiter
    iter  = 1+iter;
    currlength = currlength + 1;
    %%%%% FORWARD STEP %%%%%%
    
    resid = Phi.*x-y; % residual
    gradf = resid;  % gradient
    
    % pick next atom
    [u,v] = find_next_atom_nucnorm_matrix(gradf);
    
    Ut(:,currlength) = u;
    Vt(currlength,:) = v;
    
%     At(:,currlength) = next_atom; % update the atomic set
    
    % conditional gradient step
    phinext=Phi*(tau*u*v-x);
    gamma=-trace(resid'*phinext)/trace(phinext'*phinext);
    coeft=(1-gamma)*coeft;
    coeft(currlength)=gamma*tau;
    
    % GRADIENT PROJECTION FOR MATRICES
    if gp_forward==1
        [coeft(1:currlength),~]=grad_proj(Phi*At(:,1:currlength),y,coeft(1:currlength),tau,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x = Ut(:,1:currlength)*coeft(1:currlength,1:currlength)*Vt(1:currlength,:);
    
    obj_fwd = 0.5 * norm(Phi.*x-y,'fro')^2;
    
    if ~backward
        % update output variables
        time(iter) = toc;
        obj(iter)  = obj_fwd;
        continue
    end

    %%%%%% BACKWARD STEP  %%%%%%
    if backward
        
%         Tx = reshape(x,n1,n2);
        [U,S,V,~,~] = lansvd(x);
        S = diag(S);
        if exist('svt','var')
            remove = find(S<=svt); % SVT thresholding step
            S(remove) = [];
            U(:,remove) = [];
            V = V';
            V(remove,:) = [];
        else
            cutoff = find(diff(S) == max(diff(S)));
            S = S(1:cutoff);
            U = U(:,1:cutoff);
            V = V';
            V = V(1:cutoff,:);
        end

        numats = length(S);
        A_new = zeros(rows(U)*cols(V),numats);
        for ats = 1:numats
            a = U(:,ats)*V(ats,:);
            A_new(:,ats) = a(:);
        end
        [coeft_new,~]=grad_proj(Phi*A_new,y,zeros(numats,1),tau,gpiter,gptol);
        x_rem = A_new*coeft_new;
        obj_back = 0.5* norm(y - Phi*x_rem)^2;
        if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
            currlength = size(A_new,2);
            At(:,1:currlength) = A_new;
            x = x_rem;
            coeft(1:currlength) = coeft_new;
            %                 history = Hnew;
            obj(iter) = obj_back;
            time(iter) = toc;
            back_count = back_count +1;
            if verbose
                fprintf('#backward steps = %d \n', back_count);
            end
            if currlength<=1
            end
        else
            obj(iter) = obj_fwd;
            time(iter) = toc;
        end
        
    end
    
    
    if (obj(iter) < tol)  % convergence check
        if verbose
            fprintf('\n Convergence achieved. Final objective value = %f \n',obj(iter));
        end
        break
    end
    
end

% final output variables
At = At(:,1:currlength);

if sparsify
    if verbose
        fprintf('\n Entering final sparsification step \n');
    end
    % remove the zeros
    nonzers = find(coeft ~= 0);
    At = At(:,nonzers);
    coeft = coeft(nonzers);
    x = At*coeft;
    
end

% debias
if debias
    if verbose
        fprintf('\n Entering Debiasing Phase. %d atoms seleted \n',cols(At));
    end
    coeft = (Phi*At)\y;
    x = At*coeft;
end

end

