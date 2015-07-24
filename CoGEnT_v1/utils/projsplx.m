function x = projsplx(y,tau)
% project an n-dim vector y to the scaled simplex Dn
% Dn = { x : x n-dim, 1 >= x >= 0, sum(x) = tau}

% (c) Xiaojing Ye
% xyex19@gmail.com
%
% Algorithm is explained as in the linked document
% http://arxiv.org/abs/1101.6081
%
% Jan. 14, 2011.
% Modified by Nikhil Rao. Dec 20 2013

m = length(y); bget = false;

s = sort(y,'descend'); tmpsum = 0;

for ii = 1:m-1
    tmpsum = tmpsum + s(ii);
    tmax = (tmpsum - tau)/ii;
    if tmax >= s(ii+1)
        bget = true;
        break;
    end
end
    
if ~bget
	tmax = (tmpsum + s(m) -tau)/m; 
end

x = max(y-tmax,0);

return;