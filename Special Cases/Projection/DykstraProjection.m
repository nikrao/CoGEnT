function [x, lambda, mu, nu, iter] = DykstraProjection(z,B,u)
% Returns the vecter x that solves
%   min   || x-z ||^2
%   s.t.  sum x_i <= B   (lambda)     ---- (1)
%         0 <= x_i <= u  (mu, nu)     ---- (2)
%
% The x is the projection of z onto the polygon defined by (1) and (2)
%
% This is Dykstra's Algorithm. (NOT Dijkstra's Algorithm)
%                               ~~~~~~~~~~~~~~~~~~~~~~~~

ndim = length(z);
x = z; y=zeros(ndim,1); p = zeros(ndim,1); q = zeros(ndim,1);
quit = false;

iter = 0;
while(~quit)
  iter = iter +1;
  y_ = max(0, min(u, x+p));
  p_ = (x - y_) + p;
  [x_, lambda] = ScaledSimplexProjection(y_+q, B);
  q_ = (y_ - x_) + q;

  quit = ((norm(x_-x) + norm(y_-y) + norm(p_-p) + norm(q_-q)) < 1e-12);
  x = x_; y = y_; p = p_; q = q_;
  % fprintf('\r Iter %d..', iter);
end

if nargout > 1
  idx_ld = find(x>1e-9 & u-x>1e-9);
  if length(idx_ld) == 0
    lambda = 0;
  else
    lambda = mean(z(idx_ld)-x(idx_ld));
  end
  temp = x-z+lambda;
  idx_mu = find(temp>0); idx_nu = find(temp<0);
  mu = zeros(ndim,1); nu = zeros(ndim,1);
  mu(idx_mu) = temp(idx_mu); nu(idx_nu) = -temp(idx_nu);
end
% if nargout > 1
%   if ~(sum(x) < B)
%     [~, lambda] = ScaledSimplexProjection(z, B);
%   else
%     lambda = 0;
%   end
% end
% fprintf('\n');
