function [num, mindist, normf1, matchedframe, index] = compareframe(f1, f2, tolerance)
% compareframe - Compares each column of the matricies and looks
%       for matches within the tolerance
%
% [num, mindist, normf1, matchedframe, index] = compareframe(f1, f2, tolerance)
%
%
% f1        - Frame 1 (usually Aorig)
% f2        - Frame 2 (usually learned A)
% tolerance - Column dot product must be within tolerance of 1
%
%
% Returns
% num       - Number of columns within tolerance
% mindist   - Distance of closest matching vector for each column (vector)
% normf1    - 
% matchedframe - Frame with columns of f2 rearranged to match f1
% index     - Index of columns used to match
%
%
% JFM   7/25/2000
% Rev: 11/10/2000

[rows cols] = size(f1);

if(size(f1) ~= size(f2))
    disp('Input matricies must be the same size');
    return;
end

% Normalize
for i = 1:cols
    c = f1(:,i);
    c = c / norm(c);
    f1(:, i) = c;
end

for i = 1:cols
    c = f2(:,i);
    c = c / norm(c);
    f2(:, i) = c;
end

match = zeros(cols, 1);
for i = 1:cols
    c1 = f1(:, i);
    
    dist = zeros(cols, 1);
    for j = 1:cols
        c2 = f2(:, j);
        x = 1 - abs(c2' * c1);
        dist(j) = abs(x);
    end
    
    % Don't match a column that has already been matched
    for j = 1:cols
        if(match(j) == 1)
            dist(j) = 1;
        end
    end
    
    [t, index1] = min(dist);
    matchedframe(:, i) = f2(:, index1);
    index(i) = index1;
    
    mindist(i) = t;
    
    if(abs(t) < tolerance)        
        match(index1) = 1;
    end
    
end

normf1 = f1;
num = sum(match);
