function [signals,PC,V] = pca1(data,flag) 
% data - MxN matrix of input data 
% (M dimensions, N trials) 
% signals - MxN matrix of projected data 
% PC - each column is a PC 
% V - Mx1 matrix of variances 
[M,N] = size(data); 
% subtract off the mean for each dimension 
mn = mean(data,2); 
if nargin <2 %didn't pass flag of zero mean
    flag = 0;
end
if flag==0
        data = data - repmat(mn,1,N); 
end
% calculate the covariance matrix 
covariance =  data * data'; 
% find the eigenvectors and eigenvalues 
[PC, V] = eig(covariance); 
% extract diagonal of matrix as vector 
V = diag(V); 
% sort the variances in decreasing order 
[junk, rindices] = sort(-1*V); 
V = V(rindices); 
PC = PC(:,rindices); 
% project the original data set 
signals = PC' * data;
end