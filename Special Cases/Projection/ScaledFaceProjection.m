function [x,lambdaopt] = ScaledSimplexProjection(z,B)
% returns the vector x in R^n lying in the scaled probability simplex
% (e'x=B, x>=0 for a given B>0) that's closest to a given vector z,
% together with the Lagrange multipler lambda.
%
% SJW 10/1/08; modified 9/6/11
%
% first, order the vector

w=z/B; nw=length(w);
w = sort(w,'descend');
wcum = cumsum(w);
lambdavec = (wcum-1.0)./[1:nw]';

wtest=[w(2:nw);-inf];
lambdadiff = lambdavec - wtest;
lambdadiff = sign(lambdadiff);
if lambdadiff(1)==1
  inx=0;
else
  lambdadiff = diff(lambdadiff);
  inx = find(lambdadiff ~= 0,1);
  % if length(inx)~=1
  %   fprintf(1,' error: inx does not have exactly one element\n');
  %   return;
  % end
end
lambdaopt = lambdavec(inx+1);

lambdaopt = B*lambdaopt;
x = max(0,z-lambdaopt);

