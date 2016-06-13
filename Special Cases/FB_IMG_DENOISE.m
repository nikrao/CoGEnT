%%% script to denoise an image using SOSLASSO and FOBA

clear;
clc;
close all;

% image and patch parameters
imsize = 128;
psize = 32;
noise = 0.1;

% clustering parameters
k = 500;

% group parameters
gsize = 10;
shift = 6;

% read data
I = imread('cameraman.tif');
I = I(21:148, 71:198);
I = im2double(I);
I_noi = I + noise*randn(size(I));

% extract patches
[m n] = size(I_noi); overlap = psize-2;
Ind_x = 1:psize - overlap : m-psize+1;
if (mod(m-1,psize - overlap)~=0)
    Ind_x = [Ind_x, m-psize+1];
end
Ind_y = 1:psize - overlap : n-psize+1;
if (mod(n-1,psize - overlap)~=0)
    Ind_y = [Ind_y, n-psize+1];
end
patches = extract_patches(I_noi, Ind_y, Ind_x, psize);

% cluster data
[Labels,~] = kmeanspp(patches,k);

% make the wavelet matrix
W = MakeHaarBases(psize^2,1);
W = sparse(W);

% form groups
G = FormGroups('pc',1,psize^2); % parent child pairs

% foba parameters
maxiter = 1000;
eta = 0.5;
gptol = 1e-3;
gpiter = 50;
dropcount = 20;
do_fw = 0;
gp_forward = 1;
tau = 50;
gamma = 5;
selfun =@(gradf) find_next_atom_soslasso(gradf,G,gamma);

% run FOBA for each cluster

recovered_patches = zeros(size(patches));
W = MakeHaarbases(psize^2,1); % this is the INVERSE DWT matrix
for clus = 1:k
    ind = find(Labels==clus);
    P = patches(:,ind); % these are the patches
    
    
    
    

