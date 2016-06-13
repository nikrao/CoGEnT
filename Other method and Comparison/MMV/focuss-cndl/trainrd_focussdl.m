function [rmse, A, x, diversity, x100norm] = trainrdnocolnorm(y, x, params, Ainit, targetdiv, Aorig, xorig)
% traidrdnocolnorm - Learn a random dictionary and sparse solutions.  See
%           Section 3.2 of Kreutz 2000 (July).  Changed from trainrd to use only
%           a subset (chosen randomly) of x's to use in the A update during the
%           iteration.  No column normalizaiton
%
%   y       - Training vectors Ax = y
%   x       - Input x's (if [], a pseudoinverse solution is used to init)
%   params  - Struct of parameters
%   Ainit   - Initial choice of A matrix (can be random)
%   targetdiv - Target (original) diversity
%
% Returns:
%   rmse    - Vector of rmse at each iteration
%   A       - Learned dictionary
%   x       - Learned sparse represntations
%   diversity - Diversity of learned x vectors at each iteration
%   x100norm - p-norm of vector x_100 at each iteration
%
%  Copyright 2005 Joseph F. Murray
%
% JFM   4/5/2001
% Rev:  8/19/2002

% Save the results to a .mat file if a filename was given
outputfilename = params.outputfilename;
if(isempty(outputfilename) )
   saveresults = 0;
   disp('Not saving output');
else
   saveresults = 1;
end

saveatiter = params.saveatiter;
maxiter = params.maxiter;
Atolerance = params.Atolerance;
Xtolerance = params.Xtolerance;

% Choice of diversity measure
p = params.p;

% Learning rate, gamma for A matrix, alpha for x's, mu for balance
alpha = params.alpha;
gammaA = params.gammaA;
mu = params.mu;
lambdamax = params.lambdamax;

% Reinitialize x_k's if the diversity is too high after many iterations
reinit = params.reinit;
reinititer = params.reinititer;

display = params.display;

% Number of times to update x vectors after each update of A matrix
% For algorithm in Kreutz 2000, this should be 1.
xiter = params.xiter;

sparsifyiter = params.sparsifyiter;

% Whiten input data if necessary
if(params.whiten == 1)
    [y, whiten, dewhiten] = whitendata(y, []);
    
    if(strcmp(params.Ainit,'btrain') == 1)
        Ainit = y(:, 1:size(Aorig,2));  % First columns of data vector   
        Ainit = Ainit / norm(Ainit, 'fro');
    end
    
end


% N = number of training vector
A = Ainit;
[m N] = size(y);
[m n] = size(Ainit);

ysigma = sqrt(var(reshape(y, m*N, 1)));

%deltax = zeros(maxiter, N);
deltax_all = zeros(N, 1) + 1;
diversity = zeros(maxiter, N);

% Create original x's
if(isempty(x))
    disp('Using pseudoinverse to initialize x');
    x = pinv(A) * y;
end

Ahistorycount = 1;
flopscount = flops;

for iter = 1:maxiter
    % Lambda as a function of the iteration number
  %  lambda = 1.0e-4;
  %   lambda = (lambdamax/2)*(tanh(0.005*(iter-350))+1);
  %   lambda = (lambdamax)*(tanh(0.0025*(iter-700))+1);    % 10/20/2000 
  %  lambda = (lambdamax/350)*(iter);    % 10/24/2000 
    if(reinit == 0)
        lambdaiter = iter;
        lambda = (lambdamax)*(tanh(0.001*( lambdaiter - 1500))+1);    % 10/25/2000 

    else
        %lambdaiter = mod(iter, reinititer);
        lambdaiter = iter;
        lambda = (lambdamax)*(tanh(0.001*( lambdaiter - 1500))+1);    % 10/25/2000 
        if(iter >= reinititer)
        %    lambda = lambdamax;
        end
        
    end
    
    lambda = lambdamax;     % 2/13/2001
  %  lambda = (lambdamax)*(tanh(0.001*( lambdaiter - 1500))+1);    % 10/25/2000 

  %  lambda = (lambdamax/2)*(tanh(0.005*(-iter+500))+1);   % Flipped right-left, didn't work well
    lambdaseries(iter) = lambda;
    
    % Calculate MSE
    resid = A*x - y;
    for k = 1:N
        err(k)=norm(resid(:,k))^2;
        normyk = norm(y(:, k));
        
        if(normyk == 0) 
            epsilon(k) = 0;
        else 
            epsilon(k) = sqrt(err(k)) / normyk;
        end
        
    end
    
    mse(iter)=sum(err)/(m*N);
    mse(1) = 1;
    rmse(iter) = sqrt(mse(iter))/ysigma;
    
    
    x100resid(iter) = err(100);
    if(p ~= 0)
        % x100norm(iter) = norm(x(:,100));        % 2-norm
        x100norm(iter) = sum(abs(x(:,100)).^p);   % p-norm that we're interested in
    else
        x100norm(iter) = sum( abs(x(:,100)) > 1e-4);   % p = 0 norm with threshold
    end
           
    % Calculate numerosity/diversity
    diversity(iter, :) = numerosity(x);
    avgdiversity(iter) = mean(diversity(iter, :));
    
    % Match A's and x solutions
    if(params.compareA == 1 & mod(iter, 5) == 1)
        [nummatch(iter), mindist, normf1, matchedframe, indexx] = compareframe(Aorig, A, Atolerance);
        [dist, nummatchx(iter)] = comparesolution(xorig, x, indexx, Xtolerance);
    else
        if(iter > 1)
            nummatch(iter) = nummatch(iter-1);
            nummatchx(iter) = nummatchx(iter-1);
        else
            nummatch(iter) = 0;
            nummatchx(iter) = 0;
        end      
    end
    
    % Perturb A's at iteration
  %  if(iter == 5000)
  %      for k = 1:n
  %          if(mindist(k) > 0.05)
  %              A(:, k) = A(:, k) + 0.05 * randn(m, 1);
  %          end
  %      end
  %      
  %      % Normalize A Frobenius
  %      A = A ./ norm(A, 'fro');
  %  end
    
    if(params.Asavehistory ~= 0)
        if(mod(iter, params.Asavehistory) == 0)
            if(params.whiten == 1) 
                Ahistory{Ahistorycount} = dewhiten * A;
            else
                Ahistory{Ahistorycount} = A;
            end
                            
            Ahistorycount = Ahistorycount + 1;
        end
    end
    
  
    % -- Plot results --
    if(display == 1) 
        figure(2);
        subplot(5, 1, 1);
        semilogy(mse); 
        %plot(mse);
        plot(2:iter, rmse(2:iter));
        title(sprintf('Iteration %d |A|_F=%5.3f, A = (%d x %d), p = %5.3f, RMSE = %f', ...
            iter, norm(A,'fro'), m, n , p, rmse(iter) ));  
       % xlabel('Iteration');
        ylabel('MSE');
    
        subplot(5, 1, 2);
        
        plot(avgdiversity);
        title(sprintf('Average sparsity = %f', avgdiversity(iter)) );
       % xlabel('Iteration');
        ylabel('Average Diversity');
        %if(iter > 1)
        %    semilogy(deltax(iter-1, :), '.');
        %end
        
        % title('Normalized x updates');
       % semilogy(err, '.');
       % title('Residules');
    
        subplot(5, 1, 3);
        plot(nummatch);
        title(sprintf('Matching columns in A, alpha = %f, gamma = %f, mu = %f', ...
            alpha, gammaA, mu) );
    
        subplot(5, 1, 4);
        plot(mindist, '.');
        title('Minimum distance');
        axis([1 n 0 1]);
    
        subplot(5, 1, 5);
        plot(diversity(iter, :), '.');
        title(sprintf('Sparsity %d > %d, lambda = %f, lambdamax = %f', ...
            sum((diversity(iter,:) > targetdiv)), targetdiv, lambda, lambdamax ));
    else
        disp(sprintf('Iter %d RMSE = %f  nummatch = %d  Num > = %d  Avg Div = %f', ...
            iter, rmse(iter), nummatch(iter), sum((diversity(iter,:) > targetdiv)),...
            avgdiversity(iter) ) );
    end
    
    % -- Save results --
    if(saveresults == 1 & (mod(iter, saveatiter) == 0) )
      %  [fid, message] = fopen('block', 'r');
      %  if(fid ~= -1)
      %      fclose(fid);
      %      disp('Block file found, not writing output');
      %  else            
            disp(sprintf('Saving %s...', outputfilename) );
            mser = mse;          % Matlab bug, won't load properly
            blocking = 'true';
            % save block blocking -ASCII
      %      !ls >block
            save(outputfilename); 
      %      !rm block   
      %  !cp trainrdresults.mat in.mat
            disp('Done saving');
      %  end
  
    end
     
 
 % -- Algorithm begins here --
 
    % Reinitialize x_k's if the diversity hasn't fallen fast enough
    if(reinit == 1 & (mod(iter,reinititer) == 0) )
    %    x = pinv(A) * y;    
    %    gamma = 0;
    
        [s, index] = sort(diversity(iter, :));
        index = fliplr(index);
        s = fliplr(index);
        for k = 1:(.15 * N)
            if(s(k) > targetdiv)
                x(:, index(k)) = randn(n,1);
            end
        end
    
    %    for k = 1:N
    %        if((deltax_all(k,1) < params.focussconverge) & (diversity(iter,k) > targetdiv + 5) )
    %            x(:, k) = 10 * randn(n,1);
    %            deltax_all(k,1) = 1;
    %        end
    %    end        
    end
    
    index_xsubset = 1;    
    clear xsubset;
    clear ysubset;

    % New estimate of x_k's    
    for indexk = 1:N
        if(params.timeticks == 1 & mod(indexk, 10) == 0)
            disp(sprintf('...Iteration %d, x_%d, flops %d...', iter, indexk, flops - flopscount));
            flopscount = flops;
        end
        
        %k = indexk;           % Use vectors in order        
        k = ceil(N * rand);    % Use to randomize the order of x vectors
        
        
        if(gammaA ~= 0 | deltax_all(k, 1) > params.focussconverge)      
                                   % This is the stopping condition for FOCUSSS
            for l = 1:xiter
                xk = x(:,k);
                % Pi is the scaling matrix, see Kreutz 1997 Table 1
                invPi = diag( abs(xk).^(2-p) );
                
                if(diversity(iter, k) < 2)  % < targetdiv )
                    lambdax = 1e-6;
                %    lambdax = lambda;
                else
                    lambdax = (1-min(epsilon(k),1)) * lambda;   % 4/16/2001  New Lambda rule
                  %  lambdax = lambda;   % 9/26/2001
                end
                               
                temp = inv(lambdax * alpha * eye(m) + A * invPi * A' + 1e-10 * eye(m));
                xk1 = invPi * A' * temp * y(:, k);
            
                x(:, k) = mu*xk1 + (1-mu)*xk;
                normxk = norm(xk);
                if(normxk == 0)
                    % A zero length vector was learned
                else
                    deltaxk = norm(xk1-xk)/normxk;
                end
 
                deltax_all(k, 1) = abs(deltaxk);
                
                % Stop if update to xk small enough
                if(deltax_all(k, 1) < params.focussconverge)
                %    disp(sprintf('Vector x(%d) converged', k));
                    break;
                end
                
                % Add to the matrix of x's to use in the A update step
                if(params.Aintupdate == 1)
                    xsubset(:, index_xsubset) = x(:,k);
                    ysubset(:, index_xsubset) = y(:,k);
                    index_xsubset = index_xsubset + 1;
                end
                
            %   disp(sprintf('xiter = %d delta = %f', l, norm(xk1-xk)/norm(xk) ));
            end
        end        
        
        % Update the A matrix 
        if(iter > 1 & params.Aintupdate == 1 & mod(indexk, params.Aupdateiter) == 0 )
            %disp('Begin update A');
            xp = xsubset;
            
            if(size(xp,2) ~= params.Aupdateiter)
                fprintf('!! cols of xp %d, params.Aupdateiter %d\n', size(xp,2), params.Aupdateiter);
            end
            
            % Sparsify x
            if(mod(iter, sparsifyiter) == 0) 
                for k = 1:size(xp,2)
                    xp(:, k) = (pickhighest(xp(:, k), targetdiv))';
                end
            end
                
            % New estimate of dictionary A   
            corrxx = (1/size(xp,2)) * xp * xp';
            corryx = (1/size(xp,2)) * ysubset * xp';
            
            deltaA = A * corrxx - corryx;
            
            updateA = (gammaA)*(deltaA - trace(A'*deltaA)*A);
             
            normupdateA(iter) = norm(updateA, 'fro');
            
            A = A - updateA;
            
            % Normalize columns (Frobenius norm = 1)
            %for i = 1:n
            %    nrm = norm(A(:, i) );
            %    A(:, i) = A(:, i) / (sqrt(n) * nrm);
            %end
            
            % Normalize A Frobenius
            A = A ./ norm(A, 'fro');
            
            % Reset the collection of update x's
            index_xsubset = 1;
        end
        
    end
         
end

