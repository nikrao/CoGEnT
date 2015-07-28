function [x,err]=grad_proj_matrix(A,y,x0,tau,maxits,epsilon)
% gradient projection for min ||Ax-y||_2^2 over the simplex defined by
% x>=0, e'x<=tau.

% some internal algorithmic constants
alpha_min=1.e-6;   % min allowable value of steplength alpha
alpha_red=.5;      % reduction factor for alpha
c1=1.e-3;          % used in sufficient decrease test

nrows = size(x0,1); ncols = size(x0,2);

if tau<=0
    fprintf(1,' *** error: gradient projection: tau<=0\n');
    err=1;
    return;
end

e=ones(nrows,ncols);  % matrix of 1s to define simplex
% project given initial point onto the feasible simplex
[x,mu,err_p]=projection_simplex_matrix(x0,e,tau);


% take gradient projection steps with initial steplenth 1 (which
% should be appropriate given that the matrix A should have
% normalized columns and nice RIP.

% get initial residual and function value
r=A*x-y;
f=0.5*(trace(r'*r));

for iter=1:maxits
    g=A'*r;
    alpha=1.0;
    stop_flag=0;
    
    while (alpha>alpha_min) & (stop_flag==0)
        [x_new,mu_new,err_p]=projection_simplex_matrix(x-alpha*g,e,tau);
        
        % check convergence tolerance
        tol_check=norm(x-x_new,'fro');
        if (alpha==1) & (iter > 1)
            
            if tol_check <= epsilon
                stop_flag=1;
            end
        end
        
        % otherwise calculate new residual and function values and
        % check for sufficient decrease
        r_new=A*x_new-y; f_new=0.5*(trace(r_new'*r_new));
        if f_new <= f + c1*trace((x_new-x)'*g)
            break;
        end
        alpha=alpha_red*alpha;
    end
    if stop_flag==1
        % time to return, successfully
        err=0;
        return;
    end
    if alpha<alpha_min
        % unsuccessful exit from step length reduction loop
        
        err=1;
        return;
    end
    % take step
    r=r_new; x=x_new; f=f_new;
end
if iter==maxits
    
    err=1;
    return;
end



