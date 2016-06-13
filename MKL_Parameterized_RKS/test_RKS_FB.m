clear;close all;

addpath(genpath('/scratch/cluster/nikhilr'));
clc;

load foo
	
algparams.tau = 100;
algparams.T = 10;


[w,Xout] = GaussMKL_RKS(y,X,SIGS,D, algparams);