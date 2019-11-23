
function [real_lowrank,real_sparse,estimate_lowrank,estimate_sparse]=theory_simulation(row,col,rank,k_num)
real_lowrank = generate_lowrank_matrix(row,col,rank);
real_sparse  = generate_sparse_matrix(row,col,k_num);
Noise = real_lowrank + real_sparse;
estimate_lowrank = MatrixSGSPCA(Noise);
estimate_sparse  = Noise - estimate_lowrank;
