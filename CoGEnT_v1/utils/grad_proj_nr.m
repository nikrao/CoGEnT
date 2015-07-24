function [x,err]=grad_proj_nr(A,y,x0,tau)
A=[A -A];
[nr nc]=size(A);
maxits=500;
epsilon=1e-5;
% gradient projection for min ||Ax-y||_2^2 over the (scaled) simplex defined by
% x>=0, e'x<=tau.
  
% some internal algorithmic constants
  alpha_min=1.e-6;   % min allowable value of steplength alpha
  alpha_red=.5;      % reduction factor for alpha
  c1=1.e-3;          % used in sufficient decrease test
  
  n=length(x0); m=size(A,1);
  
  % do some sanity checking on the inputs
  if (size(A,2) ~= n) || (length(y)~=m)
    fprintf(1,' *** error: gradient projection: dimensions do not match\n');
    err=1;
    return;
  end
  if tau<=0
    fprintf(1,' *** error: gradient projection: tau<=0\n');
    err=1;
    return;
  end
  
  % evaluate and print function at given input x0 (before projection)
       
%    e=ones(n,1);  % vector of 1s to define simplex  
  % project given initial point onto the feasible simplex
%     [x,~,~]=projection_simplex(x0,e,tau);
    x = projsplx(x0,tau);
%       x= proj_simplex(x0,tau);
  
  % evaluate and print function at projected initial point
  
  % take gradient projection steps with initial steplenth 1 (which
  % should be appropriate given that the matrix A should have
  % normalized columns and nice RIP.
  
  % get initial residual and function value
  r=A*x-y;
  f=0.5*(r'*r);
  
  for iter=1:maxits
    g=A'*r;
    alpha=1.0;
    stop_flag=0;
    
    while (alpha>alpha_min) && (stop_flag==0)
%       [x_new,~,~]=projection_simplex(x-alpha*g,e,tau);
       x_new = projsplx(x - alpha*g,tau);
      
      % check convergence tolerance
      tol_check=norm(x-x_new);
      if (alpha==1) && (iter > 1)
        
        if tol_check <= epsilon
          stop_flag=1;
        end
      end
      
      % otherwise calculate new residual and function values and
      % check for sufficient decrease
      r_new=A*x_new-y; f_new=0.5*(r_new'*r_new);
      if f_new <= f + c1*(x_new-x)'*g
        break;
      end
      alpha=alpha_red*alpha;
    end
    if stop_flag==1
      % time to return, successfully
      err=0;
       xa=x(1:floor(nc/2));
         xb=x(floor(nc/2)+1:end);
        clear x;
        x=xa-xb;
      return;
    end
    if alpha<alpha_min
      % unsuccessful exit from step length reduction loop
%       fprintf(1,[' *** in gradient projection: no suitable steplength found\' ...
%                  'n']);
      err=1;
       xa=x(1:floor(nc/2));
        xb=x(floor(nc/2)+1:end);
        clear x;
        x=xa-xb;
      return;
    end
    % take step
    r=r_new; x=x_new; f=f_new;
  end
  if iter==maxits

    err=1;
     xa=x(1:floor(nc/2));
       xb=x(floor(nc/2)+1:end);
        clear x;
        x=xa-xb;
    return;
     xa=x(1:floor(nc/2));
  xb=x(floor(nc/2)+1:end);
  clear x;
  x=xa-xb;
  end
  
 
   
  