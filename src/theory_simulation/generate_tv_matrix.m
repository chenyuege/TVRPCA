function outmatrix=generate_tv_matrix(innermatrix,tensorsize,mode)
    matrix    = reshape(innermatrix,tensorsize);
    outmatrix = reshape(innermatrix,tensorsize);
    [m,n,~] = size(matrix);
    if strcmp(mode,'all_mode')
        for i = 1:m-1
            outmatrix(i,1:n-1,:)=outmatrix(i+1,1:n-1,:)+matrix(i,1:n-1,:);
        end
        for i = 1:n-1
            outmatrix(1:m-1,i,:)=outmatrix(1:m-1,i+1,:)+matrix(1:m-1,i,:);
        end
    elseif strcmp(mode,'row_mode')
        for i = 1:m-1
            outmatrix(i,:,:)=outmatrix(i+1,:,:)+matrix(i,:,:);
        end 
    else strcmp(mode,'col_mode')
        for i = 1:n-1 
            outmatrix(:,i,:)=outmatrix(:,i+1,:)+matrix(:,i,:);
        end 
    end
end