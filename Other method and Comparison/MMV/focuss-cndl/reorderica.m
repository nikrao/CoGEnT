function [xreorder, index] = reorderica(x, WA)
%  reorderica   - Reorders the rows of the data matrix learned by ica using the orignal
%                 mixing matrix and the estimated independent components.
%
% function [snr,snrk] = [xreorder, index] = reorderica(x, WA)
%               x - input signal (assumed that each signal is in one column)
%               WA - (or A * Aorig) reordering matrix 
% Returns:
%               xreorder - input data reordered to correspond to generating s
%               index - Matrix used to rearrange the data
%
%
% JFM    1/21/2002
% Rev:   1/21/2002

[n, N] = size(x);

index = zeros(n, n);

for i = 1:n
    [m, j] = max(abs(WA(i, :)));
    index(i,j) = WA(i,j);
end

index = index';
xreorder = index * x;
