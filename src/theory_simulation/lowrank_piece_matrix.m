function [outer_matrix] = lowrank_piece_matrix(row,col,rank,tensorsize)
% outer_matrix1 = normrnd(0,1/row,row,rank);
% outer_matrix2 = normrnd(0,1/col,col,rank);
outer_matrix1 = normrnd(0,1,row,rank);
outer_matrix2 = normrnd(0,1,col,rank);
outer_matrix  = outer_matrix1*outer_matrix2';
outer_tensor  = reshape(outer_matrix,tensorsize);
C_output      = outer_tensor;
figure;subplot(1,4,1);imshow(C_output(:,:,1),[])
for i =1:tensorsize/5
    for j = (i-1)*5+1:i*5
        C_output(j,:,:) = outer_tensor((i-1)*5+1,:,:);
    end
end
subplot(1,4,2);imshow(C_output(:,:,1),[])
outer_matrix = reshape(C_output,[tensorsize(1)*tensorsize(2),tensorsize(3)]);
out_tensor = reshape(outer_matrix,[tensorsize(1),tensorsize(2),tensorsize(3)]);
subplot(1,4,3);imshow(out_tensor(:,:,1),[])
