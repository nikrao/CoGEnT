clear

ndim = 100; ntests=100; failures=0;
for i=1:ntests
  z = randn(ndim,1);
  B = ndim*rand(1,1);
  u = rand(1,1);
  [x, lambda, mu, nu, iter] = DykstraProjection(z,B,u);
  % check for optimality
  opt_simplex = abs(min(B-sum(x), lambda));
  opt_lb = norm(min(mu, x));
  opt_ub = norm(min(nu, u-x));
  opttest = opt_simplex + opt_lb + opt_ub;
  if opttest>1.e-8*ndim
    failures=failures+1;
    fprintf('*');
  end
  fprintf(' test %4d: opt_simplex=%.2e, opt_lb=%.2e, opt_ub=%.2e opttest=%.2e - # Iter : %d\n', ...
             i, opt_simplex, opt_lb, opt_ub, opttest, iter);
end

if failures==0
  fprintf(' No Failures\n');
else
  fprintf(' %d Failures\n', failures);
end


