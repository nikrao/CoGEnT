function G = FormGroups(type,dim,n)

% usage G = FormGroups(type,dim)
% INPUTS:
% type = type of groups. 'pc' =>parent child pairs
%                        'pcin' => parent child pairs and in scale dependencies
%                        'path' => path from root to leaf of the tree
% dim  = dimension of the input signal
% n    = length of 1d signal. dim=2, image must be nXn
% OUTPUT
% G    = cell array of groups, each row is a group

power = log2(n);
G = cell(0);

if dim == 1 %1D signals
    switch lower(type)

        case 'pc'
            fprintf('parent-child pairs for 1D DWT \n');

            for i = 2:2^(power-1)

                group1 = {[i (2*i - 1)]};
                group2 = {[i 2*i]};
                G      = [G ; group1 ;group2];

            end
            G = [G;{1}];

        case 'pcin'
            fprintf('parent child dependencies for 1D DWT, and in scale dependencies also');

            for i = 2:2^(power-1)

                group1 = {[i (2*i - 1)]};
                group2 = {[i 2*i]};
                group3 = {[2*i 2*i-1]};
                G      = [G ; group1 ;group2;group3];

            end
            G = [G;{1}];

        case 'path'
            fprintf('paths on the 1D DWT from root to leaf \n')

            for i = n:-1:(2^(power-1)) + 1
                g = i;
                r = ceil(i/2);
                while r>=2
                    g = [g r];
                    r = ceil(r/2);
                end
                g = wrev(g);
                G = [G;{[g]}];
            end
            G = [G;{[1]}];

        otherwise
            error('this case not implemented \n')
    end
    
else
    switch lower(type)
        case 'pc'
            fprintf('parent-child pairs for 2D DWT \n');
            J = [];
            for k = 1:2^(power-1)
                for i = 1:n/2
                    if ((i+n*(k-1))==1)
                        continue
                    end               
                    J = [J i+n*(k-1)];

                end
            end            

            for i = J
                G1 = {[i 2*i]};
                G2 = {[i 2*i-1]};
                G3 = {[i 2*i+n]};
                G4 = {[i 2*i+n-1]};
                G  = [G;G1;G2;G3;G4];

            end
            
            
            
            G = [G;{1}];

        case 'pcin'
            
            fprintf('parent-child pairs for 2D DWT, and in scale dependencies \n');
            J = [];
            for k = 1:2^(power-1)
                for i = 1:n/2
                    if ((i+n*(k-1))==1)
                        continue
                    end               
                    J = [J i+n*(k-1)];

                end
            end            

            for i = J
                G1 = {[i 2*i]};
                G2 = {[i 2*i-1]};
                G3 = {[i 2*i+n]};
                G4 = {[i 2*i+n-1]};
                G5 = {[2*i 2*i-1]};
                G6 = {[2*i 2*i+n]};
                G7 = {[2*i 2*i+n-1]};
                G8 = {[2*i-1 2*i+n]};
                G9 = {[2*i-1 2*i+n-1]};
                G10 = {[2*i+n 2*i+n-1]};
                G  = [G;G1;G2;G3;G4;G5;G6;G7;G8;G9;G10];

            end
            G = [G; {[1]}];
            
            
            
        case 'path'
            fprintf('paths on the 2D DWT from root to leaf \n')
            
            leaves = [];
            leaves = [leaves n^2:-1:(n^2/2) - (n/2) + 1];
            k = leaves(end);
            seed = k-1+(n/2) : -n : k-1 +(n/2) - (n*(n/2-1));
            L = zeros(n/2);
            L(n/2,:) = seed;
            for j = n/2 - 1 : -1 : 1;
                L(j,:) = L(j+1,:) - 1;
            end
            L = L(:);
            leaves = leaves';
            leaves = [leaves ;L];
            leaves = unique(leaves); % these are the leaves of the tree
            roots = [2;2+n-1;2+n];   % these are the roots of the tree
            
            for ii = 1:3
                go = 1;
                temp = roots(ii);
                while go
                    
                    Qchildren = getChildren(temp(:,end),n);
                    tempi = [];
                    for jj = 1:length(temp(:,end))
                        tempi = [tempi; repmat(temp(jj,:),4,1)];
                    end
                    temp = tempi;
                    temp = [temp Qchildren];
                    
                    if ~isempty(intersect(Qchildren,leaves))
                        go  = 0;
                    end
                    
                    
                end
                
                for jj = 1:length(temp(:,end))
                    
                    G = [G; {[temp(jj,:)]}];
                    
                end
                
            end
            G = [G;{1}];
            
        otherwise
            error('this case not implemented \n')
    end
end
end
            
    