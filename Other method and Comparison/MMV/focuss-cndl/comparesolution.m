function [dist, nummatch] = comparesolution(xorig, xlearn, index, tolerance)
% comparesolution - Compares each column of the matricies and looks
%       for matches within the tolerance.  The x rows of the x matrix
%       are resorted according to index (if not empty).
%
% [dist, nummatch] = comparesolution(xorig, xlearn, index, tolerance)
%
%
% xlearn    - Learned solutions (each column is an x_k)
% xorig     - Correct solutions
% index     - Sorts x
% tolerance - Consider a match if norm is less than tolerance
% 
% Returns:
% dist      - Distance between normalized vectors
%
% Called from: plotresults.m
%
% JFM  10/27/2000
% Rev: 11/15/2000


% Resort the x matrix
if(~isempty(index))
    x = xlearn(index, :);
else
    x = xlearn;
end

[rows cols] = size(xorig);

nummatch = 0;

for i = 1:cols
    x1 = xorig(:,i);
    nx1 = norm(x1);
    if(nx1 == 0) 
        x1 = x1 * 0;
    else
        x1 = x1 / nx1;
    end
    x1 = abs(x1);
    
    x2 = x(:,i);
    nx2 = norm(x2);
    if(nx2 == 0)
        x2 = x2 * 0;
    else
        x2 = x2 / nx2;
    end
    x2 = abs(x2);
    
    
   % if(i == 55) 
   %     [x1 x2]
   % end
    
    %dist(i) = norm(x2-x1);
    dist(i) = abs( 1 - abs(x2'*x1) ) ;
    
    if(dist(i) < tolerance)
        nummatch = nummatch + 1;
    end
end

%dist(55)