function [x, At, iter, obj, time, back_count] = FOBA_L2(y, Phi, tau, Ainit, selfun , varargin)
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
maxiter   = 5000;
tol       = 1e-10;
backward  = 1;
eta       = 0.5;
gptol     = 1e-3;
gpiter    = 50;
verbose   = 0;
dropcount = 9999;
gp_forward= 1;
sparsify  = 0;
away      = 0;
debias    = 0;

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
            case 'away'
                away = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end



% initialize iterate
iter = 1;
x = Ainit*tau;
At = Ainit;
coeft(1) = tau;
history = Phi*Ainit*coeft;
history = Phi*Ainit;
obj = 9999*ones(1,maxiter);
obj(1) = 0.5 * norm(y - Phi*x)^2;
back_count = 0;



if verbose
    fprintf('\n Beginning FOBA iterations \n')
    figure(1)
end

drops = 0;

if away == 1
   gp_forward =0;
   backward = 0;
end

tic;
time = zeros(1,maxiter);
time(1) = toc;

while iter < maxiter
    iter  = 1+iter
    
    %%%%% FORWARD STEP %%%%%%
    
    resid = Phi*x-y; % residual
    gradf = Phi.'*resid;  % gradient
    
    % pick next atom
    next_atom = selfun(gradf);
    if away==1
        % perform away step if away = 1
        if iter>2
            away_atom = selfun(-gradf);
            dt = next_atom - x;
            da = x - away_atom;
            if dt'*gradf > da'*gradf
                next_atom = away_atom;
            end
        end
    end
    
    At = [At next_atom]; % update the atomic set
    
    % gradient projection
%     if verbose
%         fprintf('Entering Forward GP: iter = %d, obj = %f \n', iter, obj(iter));
%     end
    
    w =y - tau*Phi*next_atom;
    r = -resid;rw = r-w;
    gamma = (r'*rw)/(rw'*rw);
    coeft = (1-gamma)*coeft;
    coeft = [coeft; gamma*tau];
    
%     if do_fw==1
%         phinext=Phi*(tau*next_atom-x);
%         gamma=-(resid'*phinext)/(phinext'*phinext);
%         x=x+gamma*(tau*next_atom-x);
%         coeft=(1-gamma)*coeft;
%         coeft=[coeft; gamma*tau];
%     else
%         coeft=[coeft;0];
%     end
    
    if gp_forward==1
        [coeft,err]=grad_proj(Phi*At,y,coeft,tau,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x = At*coeft;
%     history = [history Phi*next_atom*coeft(end)];
    history = [history Phi*next_atom];
    
    obj_fwd = 0.5 * norm(resid)^2;
    
    if ~backward
        % update output variables
        time(iter) = toc;
        obj(iter)  = obj_fwd;
        continue
    end
%     if verbose
%         clf
%         stem(x); pause(0.1)
%     end
    %%%%%% BACKWARD STEP  %%%%%%
    if backward
        do_multiple = 1;
        while do_multiple
            resid = (Phi*x-y);
            gradfb =  Phi'*resid;
            
            [~, idx_rem] = min(-coeft.*(At'*gradfb)+0.5*sum(history.^2,1)');
            if idx_rem == size(At,2)
                drops = drops + 1;
            end
            
            % remove an atom and find coefficients
            A_new = At;
            A_new(:,idx_rem) = [];
            coeft_new = coeft;
            coeft_new(idx_rem) = [];
            Hnew = history;
            Hnew(:,idx_rem) = [];
%             if verbose
%                 fprintf('Entering Backward GP: iter = %d, obj = %f \n', iter, obj_fwd);
%                 
%             end
            [coeft_new,~]=grad_proj(Phi*A_new,y,coeft_new,tau,gpiter,gptol);
            
            if iter==2
                obj(iter) = obj_fwd;
                time(iter) = toc;
                do_multiple = 0;
                continue
            end
            
            x_rem = A_new*coeft_new;
            obj_back = 0.5* norm(y - Phi*x_rem)^2;
            
            if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
                At = A_new;
                x = x_rem;
                coeft = coeft_new;
                history = Hnew;
                obj(iter) = obj_back;
                time(iter) = toc;            
                back_count = back_count +1;
                if verbose
                    fprintf('#backward steps = %d \n', back_count);
                end
                if size(At,2)<=1
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

if debias
    if verbose
        fprintf('\n Entering Debiasing Phase. %d atoms seleted \n',cols(At));
    end
    coeft = (Phi*At)\y;
    x = At*coeft;
end

end

