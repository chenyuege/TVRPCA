function [outer_matrix] = generate_lowrank_piece_matrix(row,col,rank,tensorsize,mode)

% outer_matrix1 = normrnd(0,1/row,row,rank);
% outer_matrix2 = normrnd(0,1/col,col,rank);
outer_matrix1 = normrnd(0,1,row,rank);
outer_matrix2 = normrnd(0,1,col,rank);
outer_matrix  = outer_matrix1*outer_matrix2';
[m,n]         = size(outer_matrix);
% outer_tensor  = generate_tv_matrix(outer_matrix,tensorsize,mode);
outer_matrix  = reshape(outer_tensor,[m,n]);
% piece_size = min(rank,piece_size);
% outer_matrix1 = generate_piece_matrix(outer_matrix1,piece_size,mode);
% outer_matrix  = outer_matrix1*outer_matrix2';        
% outer_matrix = generate_piece_matrix(outer_matrix,piece_size,mode);