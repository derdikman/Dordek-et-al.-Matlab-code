function [Variance] = VarFunc(Conv,V)
%return the Variance of a matrix X
%Assumption: Mean is zero

% [m,n] = size(X);
% if (m>n)
%     X = X';
% end
% 
% Conv = X*X';
Variance = V*Conv*V';
end