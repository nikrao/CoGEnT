function [rmse, A, x, diversity, x100norm] = trainrd_fastica(y, x, params, Ainit, targetdiv, Aorig, xorig)
% traidrd_extica - Calls the FastICA algorithm (Hyvarinen
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
%
% JFM   1/22/2002
% Rev:  1/22/2002

addpath('\ica\fastica\fastica_21');

% Whiten input data if necessary
if(params.whiten == 1)
    [y, whiten, dewhiten] = whitendata(y, []);
end


% N = number of training vectors
A = Ainit;
[m N] = size(y);
[m n] = size(Ainit);
maxiter = params.maxiter;

% Call the FastICA learning algorithm 
% (automatically does prewhitening)
%[x,A] = ext_ica(y);
[x, Afastica, A] = fastica(y, 'displayMode', 'off', 'g', 'tanh');


% Performance (Fake the iteration dependent results)
% Calculate MSE
%ysphered = spheredata(y);
resid = A*x - y;  
for k = 1:N
    err(k)=norm(resid(:,k))^2;
end

ysigma = sqrt(var(reshape(y, m*N, 1)));

mse(1:maxiter) = sum(err)/(m*N);
rmse(1:maxiter) = sqrt(mse(maxiter))/ysigma;

% Calculate numerosity/diversity
diversity = zeros(maxiter, N);
diversity(maxiter, :) = numerosity(x);
avgdiversity(maxiter) = mean(diversity(maxiter, :));

% 2 norm (be careful when comparing to p-norm x100norm in FOCUSS algorithms)
x100norm(1:maxiter) = norm(x(:,100)); 


% Save the results to a .mat file if a filename was given
outputfilename = params.outputfilename;
if(isempty(outputfilename) )
    saveresults = 0;
    disp('Not saving output');
else
    saveresults = 1;
    save(outputfilename); 
end

