function [x,err]=gp_reduced_complex(A,y,x0,tau,maxits,epsilon)
% gradient projection for min ||Ax-y||_2^2 over the simplex defined by
% x>=0, e'x<=tau.
  
% some internal algorithmic constants
  alpha_min=1.e-6;   % min allowable value of steplength alpha
  alpha_red=.5;      % reduction factor for alpha
  c1=1.e-3;          % used in sufficient decrease test
  
  n=length(x0); m=size(A,1);
  
  % do some sanity checking on the inputs
  if (size(A,2) ~= n) | (length(y)~=m)
      size(A,2)
      length(y)
     fprintf(1,' *** error: gp_reduced: dimensions do not match\n');
    err=1;
    return;
  end
  if tau<=0
%     fprintf(1,' *** error: gp_reduced: tau<=0\n');
    err=1;
    return;
  end
  
  % evaluate and print function at given input x0 (before projection)
  rr=real(A)*x0-real(y);ri=imag(A)*x0-imag(y);
  f=0.5*(rr'*rr) + 0.5*(ri'*ri);
%   fprintf(1,'\n gp: before projection: f=%12.5e, ||x0||_1=%12.5e\n',...
%            f, norm(x0,1));
       
  e=ones(n,1);  % vector of 1s to define simplex  
  % project given initial point onto the feasible simplex
  [x,mu,err_p]=projection_simplex(x0,e,tau);
  
  % evaluate and print function at projected initial point
  rr=real(A)*x-real(y);ri=imag(A)*x-imag(y);
  f=0.5*(rr'*rr) + 0.5*(ri'*ri);
%   fprintf(1,'\n gp: after projection:  f=%12.5e, ||x||_1=%12.5e\n',...
%            f, norm(x,1));
  
  % take gradient projection steps with initial steplenth 1 (which
  % should be appropriate given that the matrix A should have
  % normalized columns and nice RIP.
  
  % get initial residual and function value
  rr=real(A)*x-real(y);ri=imag(A)*x-imag(y);
  f= ((rr'*rr) + (ri'*ri))*0.5;
  
  for iter=1:maxits
      g = real(A)'*rr + imag(A)'*ri;
%     g=A'*r;
    alpha=1.0;
    stop_flag=0;
    
    while (alpha>alpha_min) & (stop_flag==0)
      [x_new,mu_new,err_p]=projection_simplex(x-alpha*g,e,tau);
      
      % check convergence tolerance
      tol_check=norm(x-x_new);
      if (alpha==1) & (iter > 1)
        
        if tol_check <= epsilon
          stop_flag=1;
        end
      end
      
      % otherwise calculate new residual and function values and
      % check for sufficient decrease
        rrnew=real(A)*x_new-real(y);rinew=imag(A)*x_new-imag(y);
        f_new=((rrnew'*rrnew) + (rinew'*rinew))*0.5;
%       r_new=A*x_new-y; f_new=0.5*(r_new'*r_new);
      if f_new <= f + c1*(x_new-x)'*g
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
%       fprintf(1,[' *** in gp_reduced: no suitable steplength found\' ...
%                  'n']);
      err=1;
      return;
    end
    % take step
    rr = rrnew; ri = rinew;
    %     r=r_new;
    x=x_new; f=f_new;
    % print some diagnostics
%    fprintf(1,' gp: %4d: f=%12.5e, alpha=%8.3e, tol=%12.5e, ||x||_1=%12.5e\n',...
%            iter, f, alpha, tol_check, norm(x,1));
  end
  if iter==maxits
%    fprintf(1,[' *** error in gp_reduced: max iterations exceeded\' ...
%               'n']);
    err=1;
    return;
  end  