function Y = zeromean(X, dim)
% Remove empirical mean.
%
% INPUTS
% X: data matrix
% dim: direction of data (1=samples are rows, 2=samples are columns)
%
% Copyright 2013 Christian Sigg (christian@sigg-iten.ch)
% See license file for details.
%
% This code is no longer actively maintained. See my 'nsprcomp' R package
% for a better documented and more feature complete implementation.


if nargin < 2 || dim == 1
    Y = X - repmat(mean(X),size(X,1),1);
else
    Y = X - repmat(mean(X,2),1,size(X,2));
end