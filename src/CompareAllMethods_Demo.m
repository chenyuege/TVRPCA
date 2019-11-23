%==========================================================================
% This script compares hyperspectral-spectral imagery (HSI) noise removal 
% methodslisted as follows:
%   1. ALM_RPCA
%   2. WNNM_RPCA
%   3. LRMR
%   4. LRTV
%   5. EnhancedRPCA 
% three quality assessment (QA) indices -- PSNR, SSIM, ERGAS
%     -- are calculated for each methods after denoising.
%
% You can:
%       1. Type 'Demo' to to run various methods and see the all methods computed results. 
%       2. Change noise level by modifying ratio and noiselevel in Demo.m
%       3. You can change the order of add noise,when order=1,frist add
%          Gaussian noise,then S&P noise,when order=2,frist add S&P
%          noise,then Gaussian noise;
%       4. In compare methods,you can run the BM4D,ALM_RPCA,WSNM_RPCA
%          alone.
% more detail can be found in [1]
%
% [1] paper
%
% by Jiangjun peng, 2017.
%% ==================================================================
clear all;
clc;
addpath(genpath('compared methods/'))
addpath(genpath('EnhancedRPCA/'))
% initial environmental variable
EN_ALM_RPCA     = 1;  % set to 0 for turning off;
EN_WNNM_RPCA    = 1;
EN_LRMR         = 1;
EN_LRTV         = 1;
EN_EnhancedRPCA = 1;
methodname={'Noisy','ALM_RPCA','WNNM_RPCA','LRMR','LRTV','EnhancedRPCA'};

dataname = 'simu_indian.mat';
Omsi    = cell2mat(struct2cell(load(dataname)));
Nmsi    = Omsi;
[M,N,p] = size(Omsi);
%% ===========================noise simulated=======================
GaussianSwitch     = 1;
SaltPepperSwitch   = 0;
DeadlineSwitch     = 0; 
StripSwitch        = 0;
noiselevel   = 0.1*ones(p,1);%0.075*ones(p,1);0.05+0.1*rand(p,1);
ratio        = 0.15*ones(p,1);%+0.1*rand(p,1);

% Gaussian noise
if GaussianSwitch ==1
    for i =1:p
        Nmsi(:,:,i)=Omsi(:,:,i)  + noiselevel(i)*randn(M,N);
    end
end
% S&P noise
if SaltPepperSwitch ==1
    for i =1:p
        Nmsi(:,:,i)=imnoise(Nmsi(:,:,i),'salt & pepper',ratio(i));
    end
end
% dead line ---选择80-130加dead line列上加dead line
if DeadlineSwitch==1
    for i=91:130
        indp=randperm(10,1)+2;
        ind=randperm(N-1,indp);
        an=funrand(2,length(ind));
        % searching the location of an which value is 1
        loc1=find(an==1);loc2=find(an==2);loc3=find(an==3);
        Nmsi(:,ind(loc1),i)=0; 
        Nmsi(:,ind(loc2):ind(loc2)+1,i)=0;
        Nmsi(:,ind(loc3)-1:ind(loc3)+1,i)=0;
    end 
end
% dead line ---选择141-160加strip line
if StripSwitch ==1
    for band=141:160
        num = 19+randperm(21,1);
        loc = ceil(N*rand(1,num));
        t = rand(1,length(loc))*0.5-0.25;
        Nmsi(:,loc,band) = bsxfun(@minus,Nmsi(:,loc,band),t);
    end
end
% the noisy hsis quality assess value;
index  = 1;
R_image{index}  =Nmsi;
[mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
Time(index) = 0;

%% ALM_RPCA
index  = index +1;
if EN_ALM_RPCA
    lambda = 1/sqrt(M*N);
    tol=1e-6;
    maxIter=100;
    fprintf('\n');
    disp(['performing ',methodname{index}, ' ... ']);
    tic;
    D = zeros(M*N,p);
    for i=1:p 
        bandp  = Nmsi(:,:,i);
        D(:,i) = bandp(:); 
    end
    [A,~,~] = inexact_alm_rpca(D, lambda, tol, maxIter);
    R_image{index} = reshape(A,[M,N,p]);
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
end

%% WNNM_RPCA
index  = index +1;
if EN_WNNM_RPCA
    C = 0.01;
    tol=1e-6;
    maxIter=100;
    fprintf('\n');
    disp(['performing ',methodname{index}, ' ... ']);
    tic;
    D = zeros(M*N,p);
    for i=1:p 
        bandp  = Nmsi(:,:,i);
        D(:,i) = bandp(:); 
    end
    [A_hat,E_hat,iter] = inexact_alm_WNNMrpca(D,C,tol, maxIter);
    R_image{index}     = reshape(A_hat,[M,N,p]);
    Time(index)        = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
end

%% LRMR
index  = index +1;
if EN_LRMR
    r = 8;% dc_mall:4,simu_indian
    slide =20;
    s = 0.175;
    stepsize = 8;
    fprintf('\n');
    disp(['performing ',methodname{index}, ' ... ']);
    tic;
    R_image{index}    = LRMR_HSI_denoise( Nmsi,r,slide,s,stepsize );
    Time(index)       = toc;
    sizeN             = size(R_image{index});
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi(1:sizeN(1),1:sizeN(2),:),R_image{index});
end
%% LRTV
index  = index +1;
if EN_LRTV
    tau = 0.01;
    lambda =20/sqrt(M*N);%40/sqrt(M*N) for dc_pure,20 for simu_indian
    rank = 10;%10,10,13,13,dc_mall:6
    fprintf('\n');
    disp(['performing ',methodname{index}, ' ... ']);
    tic
    [R_image{index}, ~] = LRTV(Nmsi, tau,lambda, rank);
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
end

% %% LRTDTV
% index   =  index +1;
% if EN_LRTDTV
%     tau    = 0.005;
%     lambda = 40/sqrt(M*N);%1450/sqrt(M*N);
%     rank   = [160,160,6];
%     fprintf('\n');
%     disp(['performing ',methodname{index}, ' ... ']);
%     tic
%     [R_image{index},~,~,~] = LRTDTV(Nmsi, tau,lambda,rank);
%     Time(index) = toc;
%     [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
% end

%% EN_EnhancedRPCA
index   =  index +1;
if EN_EnhancedRPCA
    fprintf('\n');
    disp(['performing ',methodname{index}, ' ... ']);
    tic
    R_image{index} = SGSPCA(Nmsi);
    Time(index) = toc;
    [mpsnr(index),mssim(index),ergas(index)]=msqia(Omsi,R_image{index});
end
