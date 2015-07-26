function [x, At, obj, time,mses] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun , varargin)
% usage [x, At, obj, time] = CoGEnT_sparse(y, Phi, tau, Ainit, selfun , varargin)
% FUNCTION PERFORMS A MINIMIZATION OF THE FORM:
% min_x 0.5|| y - Phi x ||^2  s.t ||x||_A \leq tau
% this code is good for cases where the atoms are sparse (L1 norm, group
% sparse etc)
%
% INPUTS:
% Ainit       = first atom selected
% selfun      = function handle to pick an atom. one of its inputs must be the gradient.
% y, Phi      = observarion and data matrix
% tau         = constraint on the atomic norm
%
% OPTIONAL INPUTS
% 'maxiter'   = maximum number of iterations allowed default = 500
% 'tol'       = tolerance parameter
% 'backward'  = default = 1. If set, perform multiple backward steps per iteration
% 'eta'       = default = 0.5. backward parameter.
% 'gptol'     = default = 1e-3. tolerance for gradient projection
% 'gpiter'    = default = 10. maximum iterations of gradient projection
% 'gp_forward'= default = 1 perform gradient projection in forward step
% 'sparsify'  = default = 0. If set, performs a final sparsification step
%               by merging atoms and deleting zeros
% 'debias'    = default = 0; perform a final debiasing step
% 'stopcrit'  = default = 0; stopping criterion. 0 = requires the true vector. 1 = relative
%               change in obj < tol
%
% OUTPUTS :
% x           = final result
% obj         = objective function value
% time        = time taken since start of the iterations
% At          = final atomic set
%
% %%%%%%%%%%%%
% Nikhil Rao, Parikshit Shah and Stephen Wright
% Last Update: Dec 19 2013
% %%%%%%%%%%%%

% set optional parameter defaults
maxiter   = 5000;
tol       = 1e-8;
backward  = 1;
eta       = 0.5;
gptol     = 1e-3;
gpiter    = 10;
gp_forward= 1i;
sparsify  = 0;
debias    = 0;
stopcrit  = 1;

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
            case 'gp_forward'
                gp_forward = varargin{i+1};
            case 'sparsify'
                sparsify = varargin{i+1};
            case 'debias'
                debias = varargin{i+1};
            case 'stopcrit'
                stopcrit = varargin{i+1};
            case 'true_x'
                xtrue = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end

if stopcrit == 0
   if ~exist('xtrue','var')
       error('xtrue not provided \n');
   end
end

% initialize iterate
tic;
iter = 1;
x = Ainit*tau;
index = 1;
p = length(Ainit);
At = spalloc(length(Ainit),maxiter,floor(length(Ainit)*maxiter*0.25));
At(:,index) = Ainit;

coeft = zeros(maxiter,1);
coeft(index) = tau;

history = spalloc(size(Phi,1),maxiter,floor(size(Phi,1)*maxiter*0.25));
history(:,index) = Phi*Ainit;

obj = 9999*ones(1,maxiter);
obj(index) = 0.5 * norm(y - Phi*x)^2;
mses = obj;
if exist('xtrue','var')
mses(index) = norm(x - xtrue)^2/p;
end



time = zeros(1,maxiter);
time(index) = toc;

while iter < maxiter
    iter  = 1+iter;i
    index = 1+index;
    
    %%%%% FORWARD STEP %%%%%%
    
    resid = Phi*x-y; % residual
    gradf = Phi.'*resid;  % gradient
    
    % pick next atom
    next_atom = selfun(gradf);
    
    At(:,index) = next_atom; % update the atomic set
    
    % gradient projection
    if gp_forward==1
        [coeft(1:index),~]=grad_proj(Phi*At(:,1:index),y,coeft(1:index),tau,gpiter,gptol);
    end
    
    % update iterate, and store the atoms and coefficients
    x = At*coeft;
    history(:,index) = Phi*next_atom;
    
    obj_fwd = 0.5 * norm(Phi*x-y)^2;
    
    if ~backward
        % update output variables
        time(iter) = toc;
        obj(iter)  = obj_fwd;
		if exist('xtrue','var')
		
		mses(iter) = norm(x - xtrue)^2/p;
		end
        continue
    end
    
    
    %%%%%% BACKWARD STEP  %%%%%%
    if backward
        if iter==2
            obj(iter) = obj_fwd;
            time(iter) = toc;
			if exist('xtrue','var')
			
			mses(iter) = norm(x - xtrue)^2/p;
			end
            do_multiple = false;
        else
            do_multiple = true;
        end
        while do_multiple
%             resid = (Phi*x-y);
            gradfb =  Phi'*(Phi*x-y);
            
            [~, idx_rem] = min(-coeft(1:index).*(At(:,1:index)'*gradfb)+...
                0.5*sum(history(:,1:index).^2,1)'.*coeft(1:index).^2);
            
            % remove an atom and find coefficients
            A_new = At(:,1:index);
            A_new(:,idx_rem) = [];
            coeft_new = coeft(1:index);
            coeft_new(idx_rem) = [];

            [coeft_new,~]=grad_proj(Phi*A_new,y,coeft_new,tau,gpiter,gptol);
            
            x_rem = A_new*coeft_new;
            obj_back = 0.5* norm(y - Phi*x_rem)^2;
            
            if ((obj_back - obj_fwd < (obj(iter-1)-obj_fwd)*eta))
                Hnew = history(:,1:index);
                Hnew(:,idx_rem) = [];
                
                index = index - 1;
                At(:,1:index) = A_new;
                x = x_rem;
                coeft(1:index) = coeft_new;
                history(:,1:index) = Hnew;
                obj(iter) = obj_back;
                time(iter) = toc;  
				if exist('xtrue','var')
				
				mses(iter) = norm(x - xtrue)^2/p;
				end

                if index<=1
                    do_multiple = 0;
                end
            else
                obj(iter) = obj_fwd;
                time(iter) = toc;
				if exist('xtrue','var')
				mses(iter) = norm(x - xtrue)^2/p;
				end
                do_multiple = 0;
            end
            
            
        end
    end

    % convergence check
    if stopcrit == 0
        if (norm(x - xtrue)^2 < tol * norm(xtrue)^2);
            break
        end
    else
        if obj(iter-1)-obj(iter) < tol*obj(iter-1)
            break
        end
    end
    
end

At = At(:,1:index); coeft = coeft(1:index);
time = time(1:iter);
obj  = obj(1:iter);
if exist('xtrue','var')

mses = mses(1:iter);
else
	mses = 0;
end
	

if sparsify
    % remove the zeros
    nonzers = find(coeft ~= 0);
    At = At(:,nonzers);
    coeft = coeft(nonzers);
    x = At*coeft;
    
end

if debias
    coeft = (Phi*At)\y;
    x = At*coeft;
end

end

