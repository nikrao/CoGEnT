clear

ndim = 100; ntests=100; failures=0;
for i=1:ntests
  z = randn(ndim,1);
  B = ndim*rand(1,1);
  [x, lambda] = ScaledFaceProjection(z,B);
  % check for optimality
  opttest = max( norm(min(x-z+lambda,x)), ...
                 abs(B-sum(x)));
  fprintf(' test %4d: B-sum(x)=% .2e, lambda=% .2e, opttest=%.2e\n', ...
             i, B-sum(x), lambda, opttest);
  % fprintf(' test %4d: |x|_1=%.2e, lambda=%.2e, opttest=%6.2e\n', i, norm(x,1)-B, lambda, opttest);
  % fprintf(' test %4d: opttest=%6.2e\n', i, opttest);
  if opttest>1.e-8*ndim
    failures=failures+1;
  end
end

if failures==0
  fprintf(' No Failures\n');
else
  fprintf(' %d Failures\n', failures);
end


