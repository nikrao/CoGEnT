function [aout] = sortcolumns(ain)
% sortcolumns - Rearranges the columns of the matrix ain so that each column is
%               near a closes matching column (in 2-norm sense) after normalizing
%               to unit norm.  1st column stays in the same place
%
% [aout,index] = sortcolumns(ain)
%
%
% ain       - Input matrix
%
% Returns
% aout      - Output rearranged matrix (columwise arrangement) (normalized)
%
% JFM  2/5/2002
% Rev: 2/5/2002

[rows cols] = size(ain);

% Normalize
for i = 1:cols
    c = ain(:,i);
    c = c / norm(c);
    ain(:, i) = c;
end

index = zeros(cols, 1);
aout(:,1) = ain(:,1);

for i = 2:cols
    c1 = ain(:, i-1);
    mindist = 2;
    minindex = i;
    
    for j = i+1:cols
        c2 = ain(:, j);
        dist = norm(c1 - c2);
        if(dist < mindist)
            mindist = dist;
            minindex = j;
        end
    end
    
    % Swap columns to match
    aout(:,i) = ain(:,minindex);
    c2 = ain(:,minindex);
    ain(:,minindex) = ain(:,i);
    ain(:,i) = c2;
    
end
