function n_mat=mulMatCross(mat1,mat2)
%usage
 [ma,na] = size(mat1);
 [mb,nb] = size(mat2);
 
 mc = max([ma+mb-1,ma,mb]);
 nc = max([na+nb-1,na,nb]);
 
 m=min(mc,nc);
 n_mat = nan(m,m);
 
i_size = size(mat2,1);
j_size = size(mat2,2);

[work_mat,npad_i,npad_j] = pad_edges(mat1,mat2,1);

for i = 1:2*size(mat1,1)-1
    for j = 1:2*size(mat2,2)-1
        
        sub_mat = work_mat(npad_i+i-floor(i_size):npad_i+i-1, ...
            npad_j+j-floor(j_size):npad_j+j-1  );
        nan_sub_mat=sub_mat .* mat2;
        notnan_inds = find(~isnan(nan_sub_mat));  %normalized to the number of nontnan components (average)
        
        n_mat(i,j)=length(notnan_inds);
    end
end