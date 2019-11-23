clear all;clc
addpath(genpath('../compared methods/'))
addpath(genpath('../EnhancedRPCA/'))
row = 1200;
col = 100;
rank_num= col*0.05;
k_num = 0.2*row*col;
piece_size = rank_num;
tensorsize = [60,20,100];
mode = 'row_mode';
real_lowrank = generate_lowrank_matrix(row,col,rank_num);
real_lowrank = InverseTV(real_lowrank);
subplot(1,2,1);imshow(real_lowrank(1:50,1:50),[])
% real_lowrank = lowrank_piece_matrix(row,col,rank_num,tensorsize);
% real_lowrank = generate_lowrank_piece_matrix(row,col,rank_num,tensorsize,mode);
real_sparse  = generate_sparse_matrix(row,col,k_num);
Noise        = real_lowrank + real_sparse;
lambda       = 1/sqrt(max(row,col));
tol=1e-6;
maxIter=100;
fprintf('\n');
disp(['performing RPCA', ' ... ']);
tv_flag = [1,0,0];
[estimate_lowrank,~] =tv_rpca(Noise,tensorsize,tv_flag);
estimate_lowrank     = real(estimate_lowrank);
estimate_sparse  = Noise - estimate_lowrank;
outer_tensor = reshape(estimate_sparse,tensorsize);
subplot(1,2,2);imshow(estimate_lowrank(1:50,1:50),[])
% subplot(1,4,4);imshow(outer_tensor(:,:,1),[])
disp('============= print the result of RPCA: =============')
lowrank_estimation = norm(real_lowrank-estimate_lowrank,'fro')/norm(real_lowrank,'fro');
sparse_estimation  = norm(real_sparse-estimate_sparse,'fro')/norm(real_sparse,'fro');
disp(['    lowrank estimation error :    ',num2str(lowrank_estimation)]);
disp(['    sparse  estimation error :    ',num2str(sparse_estimation)]);
[mpsnr,mssim,ergas]=msqia(reshape(real_lowrank,tensorsize),reshape(estimate_lowrank,tensorsize));
disp(['  mpsnr :  ',num2str(mpsnr),'  mssim :  ',num2str(mssim),'  ergas :  ',num2str(ergas)]);
