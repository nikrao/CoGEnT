function [t,mu,err]=projection_simplex_complex(z,p,b)

% inputs:
% z (R^n): point to be projected
% p (R^n): all positive numbers, b (positive scalar), define
% the simplex <p,t><=b, t>=0.
%
% output:
% t (R^n) Euclidean projection of z onto  { t : <p,t><=b, t>=0 }
% mu: value of the Lagrange multiplier for which t(mu) is optimal.
% err: 0 if projection was found, 1 if error was detected
  
% KKT conditions are that there exists scalar mu such that
% 0 <= t--mu*p (perp) t >= 0
% 0 <= b-p'*t (perp) mu >= 0.

  n=length(p);
  % check that p has all strictly positive elements
  if ~isempty(find(p<=0))
    fprintf(1,' *** error in projection_simplex: need p>0\n');
    mu=0;
    err=1;
    return;
  end
  
  % first check to see if mu=0 gives the solution
  t=max(z,0);
%   if true
%       t = t.*sign(z);
%   end
  if (p'*t<=b) 
    mu=0;
    err=0;
    return;
  end

  % otherwise seek a positive value of mu using method outlined in
  % notes.
  
  % calculate ratio z/p,
  zp=z./p;
  % sort in ascending order
  [zpsort,isort]=sort(zp);
  % initialize pp and pz partial sums that keep track of partial sums of
  % p_i^2 and p_i*z_i
  pp_sum=p'*p; pz_sum=p'*z; found=0;
  for i=1:n
    % does the root come before the next breakpoint?
    mu_temp = (pz_sum-b) / pp_sum;
    if mu_temp <= zpsort(i);
      mu=mu_temp; found=1; break;
    end
    % No. Adjust the partial sums to remove the next terms from the
    % partial sums
    pp_sum=pp_sum-p(isort(i))*p(isort(i));
    pz_sum=pz_sum-p(isort(i))*z(isort(i));
  end
  if found==0
    fprintf(1,' ERROR: root not found\n');
    t=zeros(n,1); mu=-inf; err=1;
    return;
  end
  % calculate solution based on mu
  t=max(z-mu*p,0);
%   if true
%   %% rao modification
%   t = sign(z).*t;
%   end
  
  err=0;
  return;
  
  % We could make the code above slightly more efficient by starting the
  % search from mu=0, initializing the variables pp_sum and pz_sum to
  % take account of negative elements of z.
  
  
  
    