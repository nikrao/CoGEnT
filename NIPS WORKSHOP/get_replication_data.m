function [repindex, ginds, groupR] = get_replication_data(G)

    if isrow(G);
        G = G';
    end
    K = length(G);
    MaxGroupSize = uint32(max(cellfun(@length,G)));
    RepSpaceSize = uint32(sum(cellfun(@length,G)));
% Create RepIndex
    repindex = cell2mat(G');

% Create ginds
    temp = G;
    for i = 1:K
        temp{i}(:) = i;
    end
    ginds = cell2mat(temp');

% Create groupR
    groupR = repmat(RepSpaceSize+1,K,MaxGroupSize);
    a = 1;
    for i = 1:K
        glen = length(G{i});
        b = a+glen-1;
        groups(i,1:glen) = a:b;
        a = b+1;
    end
end
