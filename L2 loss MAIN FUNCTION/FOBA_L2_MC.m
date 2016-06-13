function [x, At, iter, obj, time, back_count] = FOBA_L2_MC(y, Phi, tau, Ainit, matsize ,svt, varargin)
% usage [x, iters, obj, time] = FOBA_L2(y, Phi, tau, x0, selfun , varargin)
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - Phi x ||^2  s.t ||x||_A \leq tau
% INPUTS:
%
% Ainit       = first atom selected
% matsize     = Size of the target matrix [n1 X n2]
% svt         = singular value threshold
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
% 'do_fw'     = default = 0 if set, do a frank-wolfe update in forward step
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
% Last Update: Sep 2 2013
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
do_fw     = 0;
gp_forward= 1;
sparsify  = 0;

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
            case 'do_fw'
                do_fw = varargin{i+1};
            case 'gp_forward'
                gp_forward = varargin{i+1};
            case 'sparsify'
                sparsify = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% check if atom selection is working OK
n1 = matsize(1);
n2 = matsize(2);
selfun =@(gradf) find_next_atom_nucnorm(gradf,n1,n2);
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
At = Ainit;
coeft(1) = tau;
history = Phi*Ainit*coeft;
obj = 9999*ones(1,maxiter);
obj(1) = 0.5 * norm(y - Phi*x)^2;
t0 = cputime;
time = zeros(1,maxiter);
time(1) = cputime - t0;
back_count = 0;
drops = 0;
while iter <= maxiter
    iter  = 1+iter;
    
    %%%%% FORWARD STEP %%%%%%
    
    gradf = Phi.'*(Phi*x - y); % gradient
    
    % pick next atom
    next_atom = selfun(gradf);
    
    At = [At next_atom]; % update the atomic set
    
    % gradient projection
    if verbose
        fprintf('Entering Forward GP: iter = %d, obj = %f \n', iter, obj(end));
    end
    
    if do_fw==1
        resid = Phi*x-y;
        gradf = Phi.'*resid;
        phinext=Phi*(tau*next_atom-x);
        gamma=min(-(resid'*phinext)/(phinext'*phinext),1);
        x=x+gamma*(tau*next_atom-x);
        coeft=(1-gamma)*coeft;
        coeft=[coeft; gamma*tau];
    else
        coeft=[coeft;0];
    end
    
    if gp_forward==1
        [coeft,err]=grad_proj(Phi*At,y,coeft,tau,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x = At*coeft;
    history = [history Phi*next_atom*coeft(end)];
    
    obj_fwd = 0.5 * norm(y - Phi*x)^2;
    
    if (~backward||(iter<=3))
        % update output variables
        time(iter) = cputime - t0;
        obj(iter)  = obj_fwd;
        continue
    end
    if verbose
        clf
        stem(x); pause(0.1)
    end
    %%%%%% BACKWARD STEP  %%%%%%
    if backward
        
        
        Tx = reshape(x,n1,n2);
        [U,S,V] = svd(Tx,'econ');
        S = diag(S);
        cutoff = find(diff(S) == max(diff(S)));
        S = S(1:cutoff);
        %         remove = find(S<=svt); % SVT thresholding step
        U = U(:,1:cutoff);
        V = V';
        V = V(1:cutoff,:);
        %         S(remove) = [];
        numats = length(S);
        S = diag(S);
        A_new = [];
        for ats = 1:numats
            a = U(:,ats)*V(ats,:);
            a = a(:);
            A_new = [A_new a];
        end
        [coeft_new,~]=grad_proj(Phi*A_new,y,zeros(numats,1),tau,gpiter,gptol);
        x_rem = A_new*coeft_new;
        obj_back = 0.5* norm(y - Phi*x_rem)^2;
        if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
            At = A_new;
            x = x_rem;
            coeft = coeft_new;
            %                 history = Hnew;
            obj(iter) = obj_back;
            time(iter) = cputime - t0;
            back_count = back_count +1;
            if verbose
                fprintf('#backward steps = %d \n', back_count);
            end
            if size(At,2)<=1
            end
        else
            obj(iter) = obj_fwd;
            time(iter) = cputime - t0;
        end
        
    end
    
    
    if drops >= dropcount
        if verbose
            fprintf('\n Dropcount achieved. Final objective value = %f \n',obj(iter));
        end
        break
    end
    %     if (obj(iter)) < tol  % convergence check
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

    coeft = (Phi*At)\y;
    x = At*coeft;

end

