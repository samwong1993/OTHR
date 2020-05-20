function Omega = covariance(sigma,M)

N = M*(M-1)/2;
Omega = zeros(N,N);

tmp = 0;
for i = 1:M-1
    Omega(tmp+1:tmp+M-i,tmp+1:tmp+M-i) = eye(M-i) + ones(M-i);
    tmp = tmp + M-i;
end

row_tmp = M-1; 
for j = 2:M-1
    col_tmp = 0; neg_col = j-2;
    for i = 1:j-1
        Omega(row_tmp+1:row_tmp+M-j,col_tmp+neg_col+1) = -ones(M-j,1);
        Omega(col_tmp+neg_col+1,row_tmp+1:row_tmp+M-j) = -ones(1,M-j);
        Omega(row_tmp+1:row_tmp+M-j,col_tmp+neg_col+2:col_tmp+neg_col+1+M-j) = eye(M-j);
        Omega(col_tmp+neg_col+2:col_tmp+neg_col+1+M-j,row_tmp+1:row_tmp+M-j) = eye(M-j);
        col_tmp = col_tmp + M-i;
        neg_col = neg_col - 1;
    end
    row_tmp = row_tmp + M-j;
end

Omegad = sigma^2*Omega;

end

