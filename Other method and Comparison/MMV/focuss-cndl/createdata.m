function [Aorig, btrain, x] = createdata(numtrain, numnonzero, rows, cols, SNRdB, divrand)
%  createdata - Gaussian data for dictionary learning. 
%               See Engan 1999 section 4
%
% numtrain - Number of training vectors (btrain) to create
% numnonzero - Number of nonzero elements in each x
% rows     - Rows in A matrix
% cols     - Cols in A matrix
% SNRdB    - Signal-to-noise ratio (10log10 var_x/var_n) , 0 for no noise
% divrand  - Add randomness to diversity if > 0, div = [numnonzero - divrand, numnonzero]
%
% Returns
% Aorig    - (rows x cols)
% btrain   - Training vectors
% x        - Sparse vectors used to generate btrain
%
%
%  JFM   7/20/2000
%  Rev:  5/4/2002


independent = 0;

% Create Aorig
Aorig = randn(rows, cols);

% Normalize Frobenius
%Aorig = Aorig / norm(Aorig, 'fro');

% Normalize columns (Frobenius norm = 1)
for i = 1:cols
    nrm = norm(Aorig(:, i) );
    Aorig(:, i) = Aorig(:, i) / (sqrt(cols) * nrm);
end

% Create sparse original x's
x = zeros(cols, numtrain);

if(numnonzero >= cols)
    disp('Numnonzero must be less than cols');
    return;
end

for i = 1:numtrain
       
    if(divrand ~= 0)
        num = numnonzero - fix(divrand * rand);
    else
        num = numnonzero;
    end

   %   Old way of creating vectors, non-independent
   if(independent == 0)
     for j = 1:num
        r = fix(cols * rand) + 1;
        while x(r, i) ~= 0
            r = fix(cols * rand) + 1;
        end
        
        x(r, i) = randn;
       %x(r,i) = rand;
        % Limit how small the components can be
        if(abs(x(r, i)) < 0.1)
            x(r, i) = 0.1;
        end                  
     end
   else
      nnz = 0;
      % New independent way of creating vectors
      for j = 1:cols
        % Is this column non-zero?
        if(rand < (num/cols))
            x(j, i) = randn;
            % Limit how small the components can be
            if(abs(x(j, i)) < 0.1)
                x(j, i) = 0.1;
            end 
            nnz = nnz + 1;
        end
      end    
    
      if(nnz == 0)
        r = fix(cols * rand) + 1;
        x(r, i) = randn;
        if(abs(x(r, i)) < 0.1)
            x(r, i) = 0.1;
        end 
      end
    end
end


% Create training vectors (b) with noise
btrain = Aorig * x;

if(SNRdB ~= 0) 
    varb = var(btrain);
    varn = 10^(-SNRdB/10);
    
    for i = 1:numtrain
        btrain(:, i) = btrain(:, i) + varb(i) * varn * randn(rows, 1);
    end
    
end



    
