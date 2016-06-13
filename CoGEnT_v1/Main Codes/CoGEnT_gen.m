function [x, obj, time] = CoGEnT_gen(y, Phi, tau, Ainit,lossfun,gradfun, selfun , varargin)

%FUNCTION PERFORMS A FRANK WOLFE OPTIMIZATION OF THE FORM
%min_x f(y,Phi x) s.t. ||x||_A \leq tau

% usage [x, At, obj, time, back_count] = CoGEnT(y, Phi, tau, Ainit, selfun , varargin)
%
% THIS CODE HAS NOT BEEN OPTIMZED FOR SPEED
%
% INPUTS:
% Ainit       = first atom selected
% lossfun     = loss. takes x as the input, user provides y and Phi
% gradfun     = function handle to provide gradient
% selfun      = function handle to pick an atom. one of its inputs must be the gradient.
% y, Phi      = observarion and data matrix
% tau         = constraint on the atomic norm
%
% OPTIONAL INPUTS
% 'maxiter'   = maximum number of iterations allowed default = 500
% 'tol'       = tolerance parameter
% 'verbose'   = default = 0. If set, gives text feedback
% 'sparsify'  = default = 0. If set, performs a final sparsification step
%               by merging atoms and deleting zeros
% 'debias'    = default = 0; perform a final debiasing step
%
% OUTPUTS :
% x           = final result
% obj         = objective function value
% time        = time taken since start of the iterations
% At          = final atomic set
% back_count  = number of backward steps taken in total
%
% %%%%%%%%%%%%
% Nikhil Rao, Parikshit Shah and Stephen Wright
% Last Update: Dec 14 2013
% %%%%%%%%%%%%

% set optional parameter defaults
maxiter   = 5000;
tol       = 1e-8;
verbose   = 0;
G 		  = 2;
lambda    = 0;
cls       = 0.5;
tls       = 0.5;

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
            case 'verbose'
                verbose = varargin{i+1};
			case 'g'
				G = varargin{i+1};
			case 'lambda'
				lambda = varargin{i+1};
            otherwise
                % Hmmm, something wrong with the parameter string
                error(['Unrecognized option: ''' varargin{i} '''']);
        end
    end
end


% initialize iterate
iter = 1;
x = Ainit*tau;
obj = 9999*ones(1,maxiter);
tic;
time = zeros(1,maxiter);
time(1) = toc;
obj(1) = lossfun(x);

while iter < maxiter
    iter  = 1+iter;
    
    %%%%% FORWARD STEP %%%%%%
    
    gradf = gradfun(x);  % gradient
	gradf = gradf + lambda*x;
    
    % pick next atom
    next_atom = selfun(gradf);
    
    % "line search" 
	gamma = G;
	j = 1; tt = cls*norm(gradf)^2;
	while (lossfun(x) - lossfun(x - gamma*gradf) < gamma*tt)
		j = j+1;
		gamma = gamma*tls;
		if j >= 10000
			break
		end
	end
	
    % update
    x = (1 - gamma)*x + gamma*tau*next_atom;
   
    % update objective etc
    obj(iter) = lossfun(x);
    time(iter) = toc;
    
    if verbose
       fprintf('iter : %d \t, loss = %f \n',iter,obj(iter)); 
    end
    
    if (obj(iter)) < tol  % convergence check
        if verbose
            fprintf('\n Convergence achieved. Final objective value = %f \n',obj(iter));
        end
        break
    end
    
end

time = time(1:iter);
obj  = obj(1:iter);


end

