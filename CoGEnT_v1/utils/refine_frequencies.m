function [fnew,cnew] = refine_frequencies(fold,cold,eps)

% function takes old set of frequencies fold, and forms a new set fnew by
% averaging 2 frequencies that are within epsilon of each other

% sort frequencies
fold = sort(fold);
fnew = fold;
cnew = cold;

tomerge = find(diff(fold)<=eps); %indices to merge
go_on = ~isempty(tomerge);
while go_on
    
    % do this if there are indices to merge
    for inds = 1:length(tomerge)
        
        k = tomerge(inds);
        fnew(k) = mean([fold(k), fold(k+1)]);
        cnew(k) = cold(k) + cold(k+1);
        fnew(k+1) = 0;
        cnew(k+1) = 0;
        
    end
    % remove zeros
    inds = find(fnew ==0);
    fnew(inds) = [];
    cnew(inds) = [];
    
    fold = fnew;
    cold = cnew;
    tomerge = find(diff(fold)<=eps); %indices to merge
    go_on = ~isempty(tomerge);
    
end

end
        
