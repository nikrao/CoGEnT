function [x, At, iter, obj, time, back_count] = FOBA_L2_temp(y, Phi, tau, Ainit, selfun , varargin)
% usage [x, iters, obj, time] = FOBA_L2(y, Phi, tau, x0, selfun , varargin)
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - Phi x ||^2  s.t ||x||_A \leq tau
% INPUTS:
%
% Ainit       = first atom selected
% selfun      = function handle to pick an atom. one of its inputs must be the gradient.
% y, Phi      = observarion and data matrix
% tau         = constraint on the atomic norm
%
% OPTIONAL INPUTS
% 'maxiter'   = maximum number of iterations allowed default = 500
% 'tol'       = tolerance parameter
% 'backward'  = default = 1. If active, backward steps are performed
% 'mbackward' = default = 0. If set, perform multiple backward steps per iteration
% 'eta'       = backward parameter.
% 'gptol'     = tolerance for gradient projection
% 'gpiter'    = maximum iterations of gradient projection
% 'verbose'   = default = 0. If set, gives text feedback
% 'dropcount' = default = 9999. The value will determine the number of
%               times an atom can be picked and dropped before we declare convergence
% 'do_fw'     = default = 1 if set, do a frank-wolfe update in forward step
% 'gp_forward'= default = 1 perform gradient projection in forward step
% NOTE : atleast one of 'do_fw' or 'gp_forward' has to be active
% 'sparsify'  = default = 0. If set, performs a final sparsification step
%               by merging atoms and deleting zeros
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
% Nikhil Rao : Dec 4 2012
% Last Update: Sep 12 2013
% %%%%%%%%%%%%

% set optional parameter defaults
maxiter   = 1000;
tol       = 1e-10;
backward  = 1;
mbackward = 0;
eta       = 0.5;
gptol     = 1e-2;
gpiter    = 10;
verbose   = 0;
dropcount = 9999;
do_fw     = 1;
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
            case 'mbackward'
                mbackward = varargin{i+1};
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
            case 'do_fw'
                do_fw = varargin{i+1};
            case 'gp_forward'
                gp_forward = varargin{i+1};
            case 'sparsify'
                sparsify = varargin{i+1};
            case 'debias'
                debias = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% check if atom selection is working OK
try
    dummy = selfun(Phi.'*y);
catch
    error('Something is wrong with function handle for atom selection')
end

% check if there is some update rule for forward step
if ((do_fw||gp_forward) == 0)
    error('There is no update rule specified for forward step \n');
end

% initialize iterate
iter = 1;
x = Ainit*tau;

% At, coeft and history need to be preallocated 
At = zeros(length(x),maxiter);
coeft = zeros(maxiter,1);
history = zeros(rows(Phi),maxiter);
currlength = iter;
% At = Ainit;
% coeft(1) = tau;
% history = Phi*Ainit*coeft;
obj = 9999*ones(1,maxiter);
obj(1) = 0.5 * norm(y - Phi*x)^2;
back_count = 0;
if verbose
    fprintf('\n Beginning iterations \n')
    figure(1)
end

drops = 0;

time = zeros(1,maxiter);
tic;
time(1) = toc;
while iter <= maxiter
    iter  = 1+iter;
    currlength = currlength + 1;
    %%%%% FORWARD STEP %%%%%%
    
    resid = Phi*x-y; % residual
    
    gradf = Phi.'*resid;  % gradient
    
    % pick next atom
    next_atom = selfun(gradf);
    
    At(:,currlength) = next_atom; % update the atomic set
    
    if do_fw==1
        phinext=Phi*(tau*next_atom-x);
        gamma=-(resid'*phinext)/(phinext'*phinext);
        coeft=(1-gamma)*coeft;
        coeft(currlength)=gamma*tau;
    end
    
    if gp_forward==1
        [temp,~]=grad_proj(Phi*At(:,1:currlength),y,coeft(1:currlength),tau,gpiter,gptol);
        coeft(1:currlength) = temp;
    end
    
    % update iterate, and store the atoms and coefficients
    x = At(:,1:currlength)*coeft(1:currlength);
    history(:,currlength) = Phi*next_atom*coeft(currlength);
    
    obj_fwd = 0.5 * norm(resid)^2;
    
    if ~backward
        % update output variables
        time(iter) = toc;
        obj(iter)  = obj_fwd;
        continue
    end
    if verbose
        clf
        stem(x); pause(0.1)
    end
    %%%%%% BACKWARD STEP  %%%%%%
    if backward
        do_multiple = 1;
        while do_multiple
            gradfb =  Phi'*(Phi*x-y);
            
            [~, idx_rem] = min(-coeft(1:currlength).*(At(:,1:currlength)'*gradfb)+0.5*sum(history(:,1:currlength).^2,1)');
            if idx_rem == currlength
                drops = drops + 1;
            end
            
            % remove an atom and find coefficients
            A_new = At(:,1:currlength);
            A_new(:,idx_rem) = [];
            coeft_new = coeft(1:currlength);
            coeft_new(idx_rem) = [];
            Hnew = history(:,1:currlength);
            Hnew(:,idx_rem) = [];
%             if verbose
%                 fprintf('Entering Backward GP: iter = %d, obj = %f \n', iter, obj_fwd);
%                 
%             end
            [coeft_new,~]=grad_proj(Phi*A_new,y,coeft_new,tau,gpiter,gptol);
            
            if i==2
                obj(iter) = obj_fwd;
                time(iter) = toc;
                continue
            end
            
            x_rem = A_new*coeft_new;
            obj_back = 0.5* norm(y - Phi*x_rem)^2;
            
            if ~mbackward
                do_multiple = 0;  % take multiple backward steps only if specified
            end
            
            if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta)) 
                % backward step taken
                currlength = currlength-1;
                At(:,1:currlength) = A_new;
                x = x_rem;
                coeft(1:currlength) = coeft_new;
                history(:,1:currlength) = Hnew;
                obj(iter) = obj_back;
                time(iter) = toc;
                back_count = back_count +1;
                if verbose
                    fprintf('#backward steps = %d \n', back_count);
                end
                if currlength<=1
                    do_multiple = 0;
                end
            else
                obj(iter) = obj_fwd;
                time(iter) = toc;
                do_multiple = 0;
            end
            
            
        end
    end
    
    
    if drops >= dropcount
        if verbose
            fprintf('\n Dropcount achieved. Final objective value = %f \n',obj(iter));
        end
        break
    end
    
    if (obj(iter)) < tol  % convergence check
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
        fprintf('\n Entering Debiasing Phase \n');
    end
    coeft = (Phi*At)\y;
    x = At*coeft;
end

end

