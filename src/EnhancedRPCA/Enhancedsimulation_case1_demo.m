clear all;
clc;	
%% simulated experiment 1
% ----------------------------load image-----------------------------------
addpath(genpath('TV_operator'))
addpath(genpath('quality assess'))
addpath(genpath('ALM_RPCA'))
load simu_indian
Ohsi       = Ori_H;

noiselevel = 0.1*ones(1,224); 
% ------------------------ Simulation experiment --------------------------
Nhsi       = Ohsi;
[M,N,p]    = size(Ohsi);
tensorsize = [M,N,p];
%% Gaussian noise
for i = 1:p
     Nhsi(:,:,i)=Ohsi(:,:,i)+noiselevel(i)*randn(M,N);
end
%% TV sparsity denoising
% output_image = SGSPCA(Nhsi);
% D = reshape(Nhsi,[M*N,p]);
% [A_hat,E_hat,iter] = inexact_alm_rpca(D);
% output_image = reshape(A_hat,[M,N,p]);
% output_image = SGSPCA(Nhsi,3/sqrt(M*N));
NoiseMatrix = zeros(M*N,p);

for i=1:p 
    bandp  = Nhsi(:,:,i);
    NoiseMatrix(:,i) = bandp(:); 
end
tv_flag = [1,1,0];
[A_hat,E_hat] = tv_rpca(NoiseMatrix,tensorsize,tv_flag);
output_image  = reshape(A_hat,tensorsize);
[mpsnr,mssim,ergas]=msqia(Ohsi,output_image);
