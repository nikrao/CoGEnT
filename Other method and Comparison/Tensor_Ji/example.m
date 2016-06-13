%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Name: LRTC (Low Rank Tensor Complement) algorithm 
% Version: 2.0
% Time: 12/29/2009
% Reference: "Tensor Completion for Estimating Missing Values in Visual Data", ICCV 2009.
% URL: http://www.public.asu.edu/~jliu87/papers/2009-ICCV-Ji-TensorCompletion.pdf
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
clc;
clear;
inData = double(imread('testImg.png','png'));
mark = inData(:,:,:)>254;         % the mark matrix indicating the position 
                                        % of the missing value: 1 means
                                        % mising, 0 means observed

maxIter = 10;                  % maximun iteration number
initial = [];                           % '[]' means no initial setting for this problem; 
                                           % otherwise the function will use the 'initial' value as the starting point

%% parameters setting %%
alpha = [1, 1, 1];                  % relaxation parameter: a larger value 
                                    % means that the solution in eqn.(12) in our paper is 
                                    % closer to eqn. (8) in our paper. In the matrix
                                    % case, eqn.(8) is exactly equivalent to
                                    % eqn.(12). You can choose any value
                                    % for \alpha
beta = [1, 1, 1];                   % you can always set its value as 1
gamma = [100, 100, 0];              % the weights of trace norm terms
  

%% tensor completion
% rImg:     recover Result 
% errList:  the function value
% R:        the real rank value of each mode
tic;
[rImg, errList, R] = LRTC(inData, alpha, beta, gamma, mark, maxIter, initial);
toc;

figure; plot(errList);                       % plot the convergence curve
figure; imshow(inData/256);         % show the original image
figure; imshow(rImg/256);            % show the result