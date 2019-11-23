function matrix=generate_piece_matrix(matrix,piece_size,mode)
    [m,n] = size(matrix);
    M_inner = ceil(m/piece_size);
    N_inner = ceil(n/piece_size);
    if strcmp(mode,'all_mode')
        for i = 1:M_inner
            for j = 1:N_inner
                row_bound = min(m,i*piece_size);
                col_bound = min(n,j*piece_size);
                matrix((i-1)*piece_size+1:row_bound,(j-1)*piece_size+1:col_bound)=matrix((i-1)*piece_size+1,(j-1)*piece_size+1);
            end
        end 
    elseif strcmp(mode,'row_mode')
        for i = 1:M_inner
            row_bound = min(m,i*piece_size);
            for j = (i-1)*piece_size+1:row_bound
                matrix(j,:)=matrix((i-1)*piece_size+1,:);
            end
        end 
    else strcmp(mode,'col_mode')
        for i = 1:N_inner
            col_bound = min(n,i*piece_size);
            for j = (i-1)*piece_size+1:col_bound
                matrix(:,j)=matrix(:,(i-1)*piece_size+1);
            end
        end 
    end
end