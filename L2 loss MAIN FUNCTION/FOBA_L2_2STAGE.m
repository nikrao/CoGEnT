function [x1,x2,At1,At2,iter, obj, time, back_count] = FOBA_L2_2STAGE...
    (y, Phi, tau1, tau2, Ainit1, Ainit2, selfun1, selfun2, varargin)

% usage [x, iters, obj, time] = FOBA_L2_2STAGE(y, Phis, taus, Ainits, selfuns , varargin)
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - Phi x1 - Phi x2 ||^2  s.t ||x1,2||_A1,2 \leq tau1,2
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
% x1,2        = final results
% iters       = iterations performed
% obj         = objective function value
% time        = time taken since start of the iterations
% At1,2       = final atomic sets
% back_count  = number of backward steps taken in total
%
% %%%%%%%%%%%%
% Nikhil Rao : Oct 7 2013
% Last Update: Oct 7 2013
% %%%%%%%%%%%%

% set optional parameter defaults
maxiter   = 5000;
tol       = 1e-10;
backward  = 1;
mbackward = 0;
eta       = 0.5;
gptol     = 1e-3;
gpiter    = 50;
verbose   = 0;
dropcount = 9999;
do_fw     = 1;
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
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

% check if atom selection is working OK
try
    dummy = selfun1(Phi.'*y);
catch
    error('Something is wrong with function handle for atom selection 1')
end
try
    dummy = selfun2(Phi.'*y);
catch
    error('Something is wrong with function handle for atom selection 2')
end

% check if there is some update rule for forward step
if ((do_fw||gp_forward) == 0)
    error('There is no update rule specified for forward step \n');
end

% initialize iterate
iter = 1;
x1 = Ainit1*tau1;
At1 = Ainit1;
coeft1(1) = tau1;
history1 = Phi*Ainit1*coeft1;
x2 = Ainit2*tau2;
At2 = Ainit2;
coeft2(1) = tau2;
history2 = Phi*Ainit2*coeft2;
obj = 9999*ones(1,maxiter);
obj(1) = 0.5 * norm(y - Phi*x2 - Phi*x1)^2;
t0 = cputime;
time = zeros(1,maxiter);
time(1) = cputime - t0;
back_count = 0;

if verbose
    fprintf('\n Beginning FOBA iterations \n')
    figure(1)
end

drops = 0;

while iter <= maxiter
    iter  = 1+iter;
    
    %%%%% FORWARD STEP FOR SET 1 %%%%%%
    
    gradf1 = Phi.'*(Phi*x1 + Phi*x2 - y);  % gradient
    
    % pick next atom
    next_atom1 = selfun1(gradf1);
    
    At1 = [At1 next_atom1]; % update the atomic set
    
    % gradient projection
    if verbose
        fprintf('Entering Forward GP: iter = %d, obj = %f \n', iter, obj(iter));
    end
    
    if do_fw==1
        resid = Phi*x1-y+Phi*x2;
        gradf = Phi.'*resid;
        phinext=Phi*(tau1*next_atom1-x1);
        gamma=-(resid'*phinext)/(phinext'*phinext);
        x1=x1+gamma*(tau1*next_atom1-x1);
        coeft1=(1-gamma)*coeft1;
        coeft1=[coeft1; gamma*tau1];
    else
        coeft1=[coeft1;0];
    end
    
    if gp_forward==1
        [coeft1,err]=grad_proj(Phi*At1,y,coeft1,tau1,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x1 = At1*coeft1;
    history1 = [history1 Phi*next_atom1*coeft1(end)];
    
    %%%%% FORWARD STEP FOR SET 2 %%%%%%
    
    gradf2 = Phi.'*(Phi*x1 + Phi*x2 - y);  % gradient
    
    % pick next atom
    next_atom2 = selfun2(gradf1);
    
    At2 = [At2 next_atom2]; % update the atomic set
    
    % gradient projection
    if verbose
        fprintf('Entering Forward GP: iter = %d, obj = %f \n', iter, obj(iter));
    end
    
    if do_fw==1
        resid = Phi*x1-y+Phi*x2;
        gradf = Phi.'*resid;
        phinext=Phi*(tau2*next_atom2-x2);
        gamma=-(resid'*phinext)/(phinext'*phinext);
        x2=x2+gamma*(tau2*next_atom2-x2);
        coeft2=(1-gamma)*coeft2;
        coeft2=[coeft2; gamma*tau2];
    else
        coeft2=[coeft2;0];
    end
    
    if gp_forward==1
        [coeft2,err]=grad_proj(Phi*At2,y,coeft2,tau2,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x2 = At2*coeft2;
    history2 = [history2 Phi*next_atom2*coeft2(end)];
    
    
    obj_fwd = 0.5 * norm(y - Phi*x1 - Phi*x2)^2;
    
    if ~backward
        % update output variables
        time(iter) = cputime - t0;
        obj(iter)  = obj_fwd;
        continue
    end
    
    
    %%%%%% BACKWARD STEP FOR SET 1  %%%%%%
    if backward
        do_multiple = 1;
        while do_multiple
            gradfb =  Phi'*(Phi*x1-y+Phi*x2);
            
            [~, idx_rem] = min(-coeft1.*(At1'*gradfb)+0.5*sum(history1.^2,1)');
            if idx_rem == size(At1,2)
                drops = drops + 1;
            end
            
            % remove an atom and find coefficients
            A_new1 = At1;
            A_new1(:,idx_rem) = [];
            coeft_new1 = coeft1;
            coeft_new1(idx_rem) = [];
            Hnew1 = history1;
            Hnew1(:,idx_rem) = [];
            if verbose
                fprintf('Entering Backward GP: iter = %d, obj = %f \n', iter, obj_fwd);
                
            end
            [coeft_new1,~]=grad_proj(Phi*A_new1,y,coeft_new1,tau1,gpiter,gptol);
            
            if i==2
                obj(iter) = obj_fwd;
                time(iter) = cputime - t0;
                continue
            end
            
            x_rem1 = A_new1*coeft_new1;
            
            
            %%%%%% BACKWARD STEP FOR SET 2  %%%%%%
            
            
            gradfb =  Phi'*(Phi*x1-y+Phi*x2);
            
            [~, idx_rem] = min(-coeft2.*(At2'*gradfb)+0.5*sum(history2.^2,1)');
            if idx_rem == size(At2,2)
                drops = drops + 1;
            end
            
            % remove an atom and find coefficients
            A_new2 = At2;
            A_new2(:,idx_rem) = [];
            coeft_new2 = coeft2;
            coeft_new2(idx_rem) = [];
            Hnew2 = history2;
            Hnew2(:,idx_rem) = [];
            if verbose
                fprintf('Entering Backward GP: iter = %d, obj = %f \n', iter, obj_fwd);
                
            end
            [coeft_new2,~]=grad_proj(Phi*A_new2,y,coeft_new2,tau2,gpiter,gptol);
            
            if i==2
                obj(iter) = obj_fwd;
                time(iter) = cputime - t0;
                continue
            end
            
            x_rem2 = A_new2*coeft_new2;
            
            obj_back = 0.5* norm(y - Phi*x_rem1 - Phi*x_rem2)^2;
            
            if ~mbackward
                do_multiple = 0;  % take multiple backward steps only if specified
            end
            
            if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
                At1 = A_new1;
                x1 = x_rem1;
                coeft1 = coeft_new1;
                history1 = Hnew1;
                At2 = A_new2;
                x2 = x_rem2;
                coeft2 = coeft_new2;
                history2 = Hnew2;
                
                obj(iter) = obj_back;
                time(iter) = cputime - t0;
                back_count = back_count +1;
                if verbose
                    fprintf('#backward steps = %d \n', back_count);
                end
                if (size(At1,2)<=1 || size(At2,2)<=1 )
                    do_multiple = 0;
                end
            else
                obj(iter) = obj_fwd;
                time(iter) = cputime - t0;
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
    if (obj(iter-1) - obj(iter)) < tol  % convergence check
        %     if (obj(iter)) < tol  % convergence check
        if verbose
            fprintf('\n Convergence achieved. Final objective value = %f \n',obj(iter));
        end
        break
    end
    
end

% if sparsify
%     if verbose
%         fprintf('\n Entering final sparsification step \n');
%     end
%     % remove the zeros
%     nonzers = find(coeft ~= 0);
%     At = At(:,nonzers);
%     coeft = coeft(nonzers);
%     x = At*coeft;
%
% end
end
