function [outer_matrix] =generate_sparse_matrix(row,col,k_num)
outer_matrix = zeros(row,col);
rand_index   = randperm(row*col);
choose_index = rand_index(1:k_num);
k_half = round(k_num/2);
outer_matrix(choose_index(1:k_half))=1;
outer_matrix(choose_index(k_half+1:k_num))=-1;
    
